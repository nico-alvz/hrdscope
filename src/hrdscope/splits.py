"""Fixed, patient-level, stratified cross-validation splits for the open benchmark.

The split files are committed so that every result in the paper and every
external contribution is computed on exactly the same patients.
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np


def make_folds(label_path: Path, out_path: Path, domain: str = "dx", n_folds: int = 5,
               seed: int = 20260918, threshold: int = 42) -> list[dict]:
    rows = list(csv.DictReader(open(label_path, newline=""), delimiter="\t"))
    key = "n_dx_slides" if domain == "dx" else "n_ts_slides"
    rows = [r for r in rows if int(r[key]) > 0]
    rng = np.random.default_rng(seed)
    out: list[dict] = []
    for label in (0, 1):
        group = [r for r in rows if int(r[f"hrd_ge{threshold}"]) == label]
        rng.shuffle(group)
        for i, r in enumerate(group):
            out.append({"patient": r["patient"], "fold": i % n_folds, "hrd_sum": r["hrd_sum"],
                        f"hrd_ge{threshold}": label, "domain": domain})
    out.sort(key=lambda r: r["patient"])
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0].keys()), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out)
    return out
