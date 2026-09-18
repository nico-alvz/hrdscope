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


def cmd_embed(args: argparse.Namespace) -> None:
    import json

    from .embed import Backbone, embed_slide

    bb = Backbone(args.backbone, device=args.device)
    out_dir = Path(args.out_dir)
    for path in args.slides:
        out = out_dir / bb.name / (Path(path).stem + ".h5")
        if out.exists() and not args.overwrite:
            print(f"skip {out.name} (exists)")
            continue
        info = embed_slide(path, out, bb, tile_px=args.tile_px, tile_mpp=args.tile_mpp,
                           batch_size=args.batch_size, max_tiles=args.max_tiles, mpp=args.mpp)
        print(json.dumps(info))


def cmd_stream(args: argparse.Namespace) -> None:
    """Download slides from a GDC manifest, embed them, delete the slide.

    The next slide is downloaded in a background thread while the current one
    is being embedded, so the GPU rarely waits for the network."""
    import csv
    import json
    import threading
    import urllib.request

    from .embed import Backbone, embed_slide

    bb = Backbone(args.backbone, device=args.device)
    rows = list(csv.DictReader(open(args.manifest, newline=""), delimiter="\t"))
    if args.strategy:
        rows = [r for r in rows if r.get("strategy", "").startswith(args.strategy)]
    if args.patients:
        keep = set(Path(args.patients).read_text().split())
        rows = [r for r in rows if r["patient"] in keep]
    rows.sort(key=lambda r: int(r.get("file_size", 0)))
    out_dir, tmp_dir = Path(args.out_dir) / bb.name, Path(args.tmp_dir)
    tmp_dir.mkdir(parents=True, exist_ok=True)
    out_dir.mkdir(parents=True, exist_ok=True)
    rows = [r for r in rows if not (out_dir / (Path(r["file_name"]).stem + ".h5")).exists()]
    log = open(Path(args.out_dir) / f"stream_{bb.name}.jsonl", "a")

    def download(r: dict) -> Path:
        slide = tmp_dir / r["file_name"]
        if slide.exists() and slide.stat().st_size == int(r.get("file_size", -1)):
            return slide
        part = slide.with_suffix(".part")
        with urllib.request.urlopen("https://api.gdc.cancer.gov/data/" + r["file_id"], timeout=600) as resp, open(part, "wb") as fh:
            while chunk := resp.read(1 << 22):
                fh.write(chunk)
        part.rename(slide)
        return slide

    class Prefetch(threading.Thread):
        def __init__(self, r):
            super().__init__(daemon=True); self.r, self.path, self.error = r, None, None
        def run(self):
            try:
                self.path = download(self.r)
            except Exception as exc:  # reported by the consumer
                self.error = exc

    nxt = Prefetch(rows[0]) if rows else None
    if nxt:
        nxt.start()
    for i, r in enumerate(rows):
        cur, nxt = nxt, (Prefetch(rows[i + 1]) if i + 1 < len(rows) else None)
        cur.join()
        if nxt:
            nxt.start()
        out = out_dir / (Path(r["file_name"]).stem + ".h5")
        slide = cur.path
        try:
            if cur.error:
                raise cur.error
            info = embed_slide(slide, out, bb, tile_px=args.tile_px, tile_mpp=args.tile_mpp,
                               batch_size=args.batch_size, max_tiles=args.max_tiles)
            info.update(patient=r.get("patient"), file_id=r["file_id"], i=i, n=len(rows))
        except Exception as exc:  # keep going, record the failure
            info = {"slide": r["file_name"], "error": repr(exc), "i": i, "n": len(rows)}
        finally:
            if slide and slide.exists() and not args.keep:
                slide.unlink()
        print(json.dumps(info), flush=True)
        log.write(json.dumps(info) + "\n")
        log.flush()


def _add_embed_args(s: argparse.ArgumentParser) -> None:
    s.add_argument("--backbone", default="midnight", choices=["midnight", "hibou-b", "h-optimus"])
    s.add_argument("--device", default=None, help="cuda or cpu (default: cuda when available)")
    s.add_argument("--out-dir", default=str(DATA / "features"))
    s.add_argument("--tile-px", type=int, default=224)
    s.add_argument("--tile-mpp", type=float, default=0.5)
    s.add_argument("--batch-size", type=int, default=32)
    s.add_argument("--max-tiles", type=int, default=None)


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

    s = sub.add_parser("embed", help="tile local slides and write foundation-model embeddings to HDF5")
    s.add_argument("slides", nargs="+")
    s.add_argument("--mpp", type=float, default=None, help="override microns per pixel when metadata is missing")
    s.add_argument("--overwrite", action="store_true")
    _add_embed_args(s)
    s.set_defaults(func=cmd_embed)

    s = sub.add_parser("stream", help="download GDC slides one by one, embed, delete")
    s.add_argument("--manifest", default=str(DATA / "manifests" / "tcga_ov_slides.tsv"))
    s.add_argument("--strategy", default="Diagnostic", help="prefix filter on the manifest strategy column")
    s.add_argument("--patients", default=None, help="file with patient ids to keep, one per line")
    s.add_argument("--tmp-dir", default=str(DATA / "raw" / "slides"))
    s.add_argument("--keep", action="store_true", help="do not delete slides after embedding")
    _add_embed_args(s)
    s.set_defaults(func=cmd_stream)

    s = sub.add_parser("splits", help="write fixed patient-level stratified folds for the benchmark")
    s.add_argument("--labels", default=str(DATA / "labels" / "tcga_ov_hrd.tsv"))
    s.add_argument("--folds", type=int, default=5)
    s.add_argument("--seed", type=int, default=20260918)
    s.add_argument("--threshold", type=int, default=42)
    s.set_defaults(func=cmd_splits)

    s = sub.add_parser("benchmark", help="cross-validated attention-MIL benchmark on the fixed splits")
    s.add_argument("--features", default=str(DATA / "features" / "midnight"))
    s.add_argument("--labels", default=str(DATA / "labels" / "tcga_ov_hrd.tsv"))
    s.add_argument("--splits", default=str(DATA / "splits" / "tcga_ov_dx_folds.tsv"))
    s.add_argument("--out-dir", default="runs/dx_midnight")
    s.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2])
    s.add_argument("--max-tiles", type=int, default=4000)
    s.add_argument("--epochs", type=int, default=40)
    s.add_argument("--device", default=None)
    s.add_argument("--params", default=None, help="best_params.json from hrdscope tune")
    s.set_defaults(func=cmd_benchmark)

    s = sub.add_parser("tune", help="nested Optuna search of aggregator hyper-parameters on the inner folds")
    s.add_argument("--features", default=str(DATA / "features" / "midnight"))
    s.add_argument("--labels", default=str(DATA / "labels" / "tcga_ov_hrd.tsv"))
    s.add_argument("--splits", default=str(DATA / "splits" / "tcga_ov_dx_folds.tsv"))
    s.add_argument("--out-dir", default="runs/tune_dx_midnight")
    s.add_argument("--trials", type=int, default=30)
    s.add_argument("--epochs", type=int, default=40)
    s.add_argument("--device", default=None)
    s.set_defaults(func=cmd_tune)

    s = sub.add_parser("predict", help="embed one slide, run the trained ensemble, write JSON and attention heatmap")
    s.add_argument("slide")
    s.add_argument("--run-dir", default="runs/dx_midnight")
    s.add_argument("--out-dir", default="predictions")
    s.add_argument("--backbone", default="midnight", choices=["midnight", "hibou-b", "h-optimus"])
    s.add_argument("--device", default=None)
    s.add_argument("--batch-size", type=int, default=16)
    s.add_argument("--mpp", type=float, default=None)
    s.set_defaults(func=cmd_predict)

    args = p.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
