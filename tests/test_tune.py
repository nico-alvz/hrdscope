import csv
import json

import h5py
import numpy as np

from hrdscope.tune import tune


def test_tune_smoke(tmp_path):
    feat = tmp_path / "feat"
    feat.mkdir()
    rng = np.random.default_rng(0)
    with open(tmp_path / "labels.tsv", "w", newline="") as fh, open(tmp_path / "folds.tsv", "w", newline="") as fs:
        wl = csv.writer(fh, delimiter="\t")
        ws = csv.writer(fs, delimiter="\t")
        wl.writerow(["patient", "hrd_sum", "hrd_ge33", "hrd_ge42", "hrd_ge63"])
        ws.writerow(["patient", "fold"])
        for i in range(16):
            p = f"TCGA-00-{i:04d}"
            hrd = 20 + 60 * (i % 2)
            with h5py.File(feat / f"{p}-01Z-00-DX1.h5", "w") as h5:
                h5.create_dataset("features", data=rng.normal(loc=(i % 2), size=(20, 8)).astype(np.float16))
            wl.writerow([p, hrd, int(hrd >= 33), int(hrd >= 42), int(hrd >= 63)])
            ws.writerow([p, i % 3])
    best = tune(feat, tmp_path / "labels.tsv", tmp_path / "folds.tsv", tmp_path / "tune", n_trials=2, epochs=3, device="cpu")
    assert set(best) == {"0", "1", "2"} and "lr" in best["0"]["params"]
    assert json.load(open(tmp_path / "tune" / "best_params.json"))
