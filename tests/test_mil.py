import csv
import numpy as np
import h5py
import torch

from hrdscope.benchmark import run_cv, triage_points
from hrdscope.mil import GatedAttentionMIL, predict


def test_model_shapes():
    m = GatedAttentionMIL(16)
    pr, pb, a = m(torch.randn(50, 16))
    assert pr.shape == () and pb.shape == (3,) and a.shape == (50,)
    assert abs(float(a.sum()) - 1) < 1e-4
    out = predict(m, torch.randn(10, 16))
    assert set(out["prob"]) == {"ge33", "ge42", "ge63"}


def test_triage_points_perfect_classifier():
    y = np.array([0] * 50 + [1] * 50)
    p = np.linspace(0, 1, 100)
    t = triage_points(y, p)
    assert t["rule_out"]["sensitivity"] >= 0.95 and t["rule_out"]["npv"] >= 0.95
    assert t["rule_in"]["ppv"] >= 0.9


def test_run_cv_smoke(tmp_path):
    feat = tmp_path / "feat"; feat.mkdir()
    rng = np.random.default_rng(0)
    pats = [f"TCGA-00-{i:04d}" for i in range(20)]
    with open(tmp_path / "labels.tsv", "w", newline="") as fh, open(tmp_path / "folds.tsv", "w", newline="") as fs:
        wl = csv.writer(fh, delimiter="\t"); ws = csv.writer(fs, delimiter="\t")
        wl.writerow(["patient", "hrd_sum", "hrd_ge33", "hrd_ge42", "hrd_ge63"]); ws.writerow(["patient", "fold"])
        for i, p in enumerate(pats):
            hrd = 20 + 60 * (i % 2)
            x = rng.normal(loc=(i % 2), size=(30, 8)).astype(np.float16)
            with h5py.File(feat / f"{p}-01Z-00-DX1.h5", "w") as h5:
                h5.create_dataset("features", data=x)
            wl.writerow([p, hrd, int(hrd >= 33), int(hrd >= 42), int(hrd >= 63)]); ws.writerow([p, i % 4])
    rep = run_cv(feat, tmp_path / "labels.tsv", tmp_path / "folds.tsv", tmp_path / "runs", seeds=(0,), epochs=5, device="cpu")
    assert rep["n_patients"] == 20 and "ge42" in rep and (tmp_path / "runs" / "oof_predictions.tsv").exists()
