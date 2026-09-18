"""Command-line interface: ``hrdscope <command>``."""

from __future__ import annotations

import argparse
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data"


def cmd_labels(args: argparse.Namespace) -> None:
    from .labels import build_label_table, concordance_report, download_gdc_files

    ascat_dir = Path(args.ascat_dir)
    n = download_gdc_files(DATA / "manifests" / "tcga_ov_ascat3.tsv", ascat_dir)
    print(f"downloaded {n} new ASCAT3 files into {ascat_dir}")
    rows = build_label_table(
        DATA / "manifests" / "tcga_ov_ascat3.tsv", ascat_dir,
        DATA / "reference" / "hg38_chrom_arms.tsv",
        DATA / "manifests" / "tcga_ov_slides.tsv",
        DATA / "manifests" / "knijnenburg2018_ov_hrd.tsv",
        Path(args.out),
    )
    print(f"wrote {args.out}")
    print(concordance_report(rows))


def cmd_score(args: argparse.Namespace) -> None:
    from .scars import load_arms, load_gdc_ascat, score

    arms = load_arms(DATA / "reference" / "hg38_chrom_arms.tsv")
    for path in args.files:
        s = score(load_gdc_ascat(path), arms)
        print(f"{path}\tLOH={s.loh}\tTAI={s.tai}\tLST={s.lst}\tHRD_sum={s.hrd_sum}")


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(prog="hrdscope", description=__doc__)
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("labels", help="download open GDC ASCAT3 profiles and build the TCGA-OV HRD label table")
    s.add_argument("--ascat-dir", default=str(DATA / "raw" / "ascat3"))
    s.add_argument("--out", default=str(DATA / "labels" / "tcga_ov_hrd.tsv"))
    s.set_defaults(func=cmd_labels)

    s = sub.add_parser("score", help="compute LOH/TAI/LST/HRD-sum for allele-specific segment files")
    s.add_argument("files", nargs="+")
    s.set_defaults(func=cmd_score)

    args = p.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
