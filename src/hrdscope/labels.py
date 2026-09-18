"""Build the open TCGA-OV HRD label table from GDC allele-specific copy number.

Every input is open access (no dbGaP): GDC ASCAT3 segments for TCGA-OV, the
Knijnenburg et al. 2018 (Cell Rep) HRD scores used as an external check, and
the GDC slide manifest. The output is one row per patient with scar components,
HRD-sum, binary labels at the three thresholds used in the literature, and
which slide types are available.
"""

from __future__ import annotations

import csv
import statistics
import urllib.request
from collections import defaultdict
from pathlib import Path

from .scars import load_arms, load_gdc_ascat, score

GDC_DATA = "https://api.gdc.cancer.gov/data/"
THRESHOLDS = (33, 42, 63)


def download_gdc_files(manifest: Path, out_dir: Path, workers: int = 8) -> int:
    """Download files listed in a manifest (columns ``file_id``, ``file_name``)."""
    import concurrent.futures as cf

    out_dir.mkdir(parents=True, exist_ok=True)
    rows = list(csv.DictReader(open(manifest, newline=""), delimiter="\t"))

    def fetch(row: dict) -> int:
        dest = out_dir / row["file_name"]
        if dest.exists() and dest.stat().st_size > 0:
            return 0
        data = urllib.request.urlopen(GDC_DATA + row["file_id"], timeout=300).read()
        dest.write_bytes(data)
        return 1

    with cf.ThreadPoolExecutor(workers) as ex:
        return sum(ex.map(fetch, rows))


def build_label_table(ascat_manifest: Path, ascat_dir: Path, arms_path: Path,
                      slides_manifest: Path, knijnenburg_path: Path | None,
                      out_path: Path) -> list[dict]:
    arms = load_arms(arms_path)
    per_patient: dict[str, list] = defaultdict(list)
    for row in csv.DictReader(open(ascat_manifest, newline=""), delimiter="\t"):
        segs = load_gdc_ascat(ascat_dir / row["file_name"])
        per_patient[row["patient"]].append(score(segs, arms))

    slides: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))
    for row in csv.DictReader(open(slides_manifest, newline=""), delimiter="\t"):
        key = "dx" if row["strategy"].startswith("Diagnostic") else "ts"
        slides[row["patient"]][key] += 1

    knij: dict[str, dict] = {}
    if knijnenburg_path and knijnenburg_path.exists():
        knij = {r["patient"]: r for r in csv.DictReader(open(knijnenburg_path, newline=""), delimiter="\t")}

    rows: list[dict] = []
    for patient in sorted(per_patient):
        scores = per_patient[patient]
        loh = statistics.mean(s.loh for s in scores)
        tai = statistics.mean(s.tai for s in scores)
        lst = statistics.mean(s.lst for s in scores)
        hrd = loh + tai + lst
        row = {
            "patient": patient,
            "n_ascat_profiles": len(scores),
            "hrd_loh": round(loh, 2), "ntai": round(tai, 2), "lst": round(lst, 2),
            "hrd_sum": round(hrd, 2),
        }
        for t in THRESHOLDS:
            row[f"hrd_ge{t}"] = int(hrd >= t)
        row["knijnenburg2018_hrd"] = knij.get(patient, {}).get("HRD_Score", "")
        row["n_dx_slides"] = slides[patient]["dx"]
        row["n_ts_slides"] = slides[patient]["ts"]
        rows.append(row)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    return rows


def concordance_report(rows: list[dict]) -> str:
    """Compare our scores with the published Knijnenburg 2018 values."""
    from scipy.stats import pearsonr, spearmanr

    paired = [(float(r["hrd_sum"]), float(r["knijnenburg2018_hrd"])) for r in rows if r["knijnenburg2018_hrd"]]
    if len(paired) < 3:
        return "no overlap with Knijnenburg 2018"
    ours, theirs = zip(*paired)
    lines = [f"patients scored: {len(rows)}; overlap with Knijnenburg 2018: {len(paired)}",
             f"HRD-sum Pearson r = {pearsonr(ours, theirs)[0]:.3f}, Spearman rho = {spearmanr(ours, theirs)[0]:.3f}",
             f"HRD-sum median ours = {statistics.median(ours):.1f}, published = {statistics.median(theirs):.1f}"]
    for t in THRESHOLDS:
        agree = sum((o >= t) == (p >= t) for o, p in paired) / len(paired)
        lines.append(f"binary agreement at HRD-sum >= {t}: {agree:.3f}")
    for t in THRESHOLDS:
        n = sum(int(r[f"hrd_ge{t}"]) for r in rows)
        lines.append(f"all patients: HRD-sum >= {t}: {n}/{len(rows)} ({100 * n / len(rows):.1f}%)")
    return "\n".join(lines)
