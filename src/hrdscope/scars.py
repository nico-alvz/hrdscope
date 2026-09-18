"""Genomic scar scores from allele-specific copy-number segments.

Implements the three classical homologous-recombination-deficiency (HRD) scar
metrics used by clinical assays and by the scarHRD R package:

* HRD-LOH  (Abkevich et al., Br J Cancer 2012): number of loss-of-heterozygosity
  regions longer than 15 Mb that do not span the whole chromosome.
* NtAI     (Birkbak et al., Cancer Discov 2012): number of regions of allelic
  imbalance that extend to the telomere but do not cross the centromere.
* LST      (Popova et al., Cancer Res 2012): number of chromosomal breaks
  between adjacent regions of at least 10 Mb, after filtering out segments
  shorter than 3 Mb.

HRD-sum = HRD-LOH + NtAI + LST is the quantity thresholded by Myriad
myChoice-like assays (commonly at 42) and by TCGA studies (33, 42 or 63).

Input is a list of ``Segment`` objects with major/minor allele copy numbers,
e.g. loaded from GDC open-access ASCAT3 files with ``load_gdc_ascat``.
Only autosomes are scored.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

AUTOSOMES = tuple(f"chr{i}" for i in range(1, 23))


@dataclass(frozen=True)
class Segment:
    chrom: str
    start: int
    end: int
    major: int
    minor: int

    @property
    def length(self) -> int:
        return self.end - self.start

    @property
    def total(self) -> int:
        return self.major + self.minor

    @property
    def is_loh(self) -> bool:
        return self.minor == 0 and self.major > 0

    @property
    def is_ai(self) -> bool:
        return self.major != self.minor


@dataclass(frozen=True)
class Arm:
    length: int
    cen_start: int
    cen_end: int


@dataclass(frozen=True)
class ScarScores:
    loh: int
    tai: int
    lst: int

    @property
    def hrd_sum(self) -> int:
        return self.loh + self.tai + self.lst


def load_arms(path: str | Path) -> dict[str, Arm]:
    arms: dict[str, Arm] = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            arms[row["chrom"]] = Arm(int(row["length"]), int(row["cen_start"]), int(row["cen_end"]))
    return arms


def load_gdc_ascat(path: str | Path) -> list[Segment]:
    """Load a GDC ``*.ascat3.allelic_specific.seg.txt`` (or ASCAT2/AscatNGS) file."""
    segs: list[Segment] = []
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            chrom = row["Chromosome"]
            if not chrom.startswith("chr"):
                chrom = "chr" + chrom
            if chrom not in AUTOSOMES:
                continue
            segs.append(
                Segment(chrom, int(row["Start"]), int(row["End"]),
                        int(row["Major_Copy_Number"]), int(row["Minor_Copy_Number"]))
            )
    segs.sort(key=lambda s: (AUTOSOMES.index(s.chrom), s.start))
    return segs


def _by_chrom(segs: list[Segment]) -> dict[str, list[Segment]]:
    out: dict[str, list[Segment]] = {}
    for s in segs:
        out.setdefault(s.chrom, []).append(s)
    for v in out.values():
        v.sort(key=lambda s: s.start)
    return out


def merge_adjacent(segs: list[Segment]) -> list[Segment]:
    """Join consecutive segments on the same chromosome with identical allele state."""
    out: list[Segment] = []
    for s in segs:
        if out and out[-1].chrom == s.chrom and (out[-1].major, out[-1].minor) == (s.major, s.minor):
            p = out[-1]
            out[-1] = Segment(p.chrom, p.start, s.end, p.major, p.minor)
        else:
            out.append(s)
    return out


def hrd_loh(segs: list[Segment], min_size: int = 15_000_000) -> int:
    """Abkevich 2012: LOH regions > ``min_size`` not spanning the whole chromosome."""
    n = 0
    for chrom_segs in _by_chrom(merge_adjacent(segs)).values():
        if all(s.is_loh for s in chrom_segs):
            continue  # whole-chromosome LOH is excluded
        run_start = run_end = None
        runs: list[int] = []
        for s in chrom_segs:
            if s.is_loh:
                if run_start is None:
                    run_start = s.start
                run_end = s.end
            elif run_start is not None:
                runs.append(run_end - run_start)
                run_start = run_end = None
        if run_start is not None:
            runs.append(run_end - run_start)
        n += sum(1 for r in runs if r > min_size)
    return n


def ntai(segs: list[Segment], arms: dict[str, Arm], min_size: int = 1_000_000) -> int:
    """Birkbak 2012: telomeric allelic imbalance not crossing the centromere."""
    n = 0
    for chrom, chrom_segs in _by_chrom(merge_adjacent(segs)).items():
        arm = arms[chrom]
        if len(chrom_segs) == 1:
            continue  # whole-chromosome AI is excluded
        first, last = chrom_segs[0], chrom_segs[-1]
        if first.is_ai and first.length >= min_size and first.end < arm.cen_start:
            n += 1
        if last.is_ai and last.length >= min_size and last.start > arm.cen_end:
            n += 1
    return n


def _split_at_centromere(chrom_segs: list[Segment], arm: Arm) -> tuple[list[Segment], list[Segment]]:
    p: list[Segment] = []
    q: list[Segment] = []
    for s in chrom_segs:
        if s.end <= arm.cen_start:
            p.append(s)
        elif s.start >= arm.cen_end:
            q.append(s)
        else:
            if s.start < arm.cen_start:
                p.append(Segment(s.chrom, s.start, arm.cen_start, s.major, s.minor))
            if s.end > arm.cen_end:
                q.append(Segment(s.chrom, arm.cen_end, s.end, s.major, s.minor))
    return p, q


def _drop_small(arm_segs: list[Segment], min_size: int) -> list[Segment]:
    """Iteratively remove segments shorter than ``min_size`` and re-merge (Popova smoothing)."""
    segs = merge_adjacent(arm_segs)
    while True:
        small = [i for i, s in enumerate(segs) if s.length < min_size]
        if not small or len(segs) == 1:
            return segs
        i = small[0]
        s = segs[i]
        if i == 0:
            nxt = segs[1]
            segs[1] = Segment(nxt.chrom, s.start, nxt.end, nxt.major, nxt.minor)
        elif i == len(segs) - 1:
            prv = segs[i - 1]
            segs[i - 1] = Segment(prv.chrom, prv.start, s.end, prv.major, prv.minor)
        else:
            prv, nxt = segs[i - 1], segs[i + 1]
            mid = (s.start + s.end) // 2
            segs[i - 1] = Segment(prv.chrom, prv.start, mid, prv.major, prv.minor)
            segs[i + 1] = Segment(nxt.chrom, mid, nxt.end, nxt.major, nxt.minor)
        del segs[i]
        segs = merge_adjacent(segs)


def lst(segs: list[Segment], arms: dict[str, Arm], min_size: int = 10_000_000,
        smooth_size: int = 3_000_000) -> int:
    """Popova 2012: breaks between adjacent segments >= ``min_size`` per chromosome arm."""
    n = 0
    for chrom, chrom_segs in _by_chrom(segs).items():
        for arm_segs in _split_at_centromere(chrom_segs, arms[chrom]):
            if not arm_segs:
                continue
            clean = _drop_small(arm_segs, smooth_size)
            for a, b in zip(clean, clean[1:]):
                if a.length >= min_size and b.length >= min_size and (b.start - a.end) < smooth_size:
                    n += 1
    return n


def score(segs: list[Segment], arms: dict[str, Arm]) -> ScarScores:
    return ScarScores(hrd_loh(segs), ntai(segs, arms), lst(segs, arms))
