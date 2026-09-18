"""Cross-validated benchmark on the fixed splits: train, predict out-of-fold, report metrics."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
import torch

from .mil import GatedAttentionMIL, THRESHOLDS, index_features, load_bag, predict, read_labels, train_one


def bootstrap_auc(y: np.ndarray, p: np.ndarray, n: int = 1000, seed: int = 0) -> tuple[float, float, float]:
    from sklearn.metrics import roc_auc_score

    rng = np.random.default_rng(seed)
    auc = roc_auc_score(y, p)
    vals = []
    for _ in range(n):
        idx = rng.integers(0, len(y), len(y))
        if len(set(y[idx])) == 2:
            vals.append(roc_auc_score(y[idx], p[idx]))
    return auc, float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


def triage_points(y: np.ndarray, p: np.ndarray, sens_target: float = 0.95, ppv_target: float = 0.90) -> dict:
    """Rule-out threshold (sensitivity >= target) and rule-in threshold (PPV >= target)."""
    order = np.argsort(-p)
    ys, ps = y[order], p[order]
    out = {}
    tp = np.cumsum(ys); n_pos = ys.sum()
    sens = tp / max(1, n_pos)
    k = np.searchsorted(sens, sens_target)
    if k < len(ps):
        thr = ps[k]
        below = p < thr
        out["rule_out"] = {"threshold": float(thr), "sensitivity": float(sens[k]),
                           "fraction_ruled_out": float(below.mean()),
                           "npv": float(((y == 0) & below).sum() / max(1, below.sum()))}
    ppv = tp / np.arange(1, len(ys) + 1)
    ok = np.where(ppv >= ppv_target)[0]
    if len(ok):
        k = ok[-1]
        out["rule_in"] = {"threshold": float(ps[k]), "ppv": float(ppv[k]), "fraction_ruled_in": float((k + 1) / len(ps))}
    return out


def run_cv(feature_dir: Path, label_path: Path, split_path: Path, out_dir: Path, seeds: tuple[int, ...] = (0, 1, 2),
           max_tiles: int = 4000, device: str | None = None, epochs: int = 40) -> dict:
    device = device or ("cuda" if torch.cuda.is_available() else "cpu")
    feats = index_features(feature_dir)
    labels = read_labels(label_path)
    splits = [r for r in csv.DictReader(open(split_path, newline=""), delimiter="\t") if r["patient"] in feats]
    if not splits:
        raise SystemExit("no patients with both features and split assignment")
    in_dim = load_bag(feats[splits[0]["patient"]][:1], max_tiles=1).shape[1]
    rng = np.random.default_rng(0)
    bags = {r["patient"]: (load_bag(feats[r["patient"]], max_tiles, rng), float(labels[r["patient"]]["hrd_sum"]),
                            np.array([int(labels[r["patient"]][f"hrd_ge{t}"]) for t in THRESHOLDS], dtype=np.float32))
            for r in splits}
    folds = sorted({int(r["fold"]) for r in splits})
    oof = {p: {"hrd_sum_pred": [], **{f"ge{t}": [] for t in THRESHOLDS}} for p in bags}
    out_dir.mkdir(parents=True, exist_ok=True)
    for seed in seeds:
        for f in folds:
            test = [r["patient"] for r in splits if int(r["fold"]) == f]
            rest = [r["patient"] for r in splits if int(r["fold"]) != f]
            val_fold = (f + 1) % len(folds)
            val = [r["patient"] for r in splits if int(r["fold"]) == val_fold]
            train = [p for p in rest if p not in set(val)]
            model = GatedAttentionMIL(in_dim)
            model, info = train_one(model, [bags[p] for p in train], [bags[p] for p in val], epochs=epochs,
                                    device=device, seed=seed * 100 + f)
            torch.save(model.state_dict(), out_dir / f"abmil_seed{seed}_fold{f}.pt")
            for p in test:
                pr = predict(model, bags[p][0], device)
                oof[p]["hrd_sum_pred"].append(pr["hrd_sum_pred"])
                for t in THRESHOLDS:
                    oof[p][f"ge{t}"].append(pr["prob"][f"ge{t}"])
            print(json.dumps({"seed": seed, "fold": f, **info, "n_train": len(train), "n_test": len(test)}), flush=True)
    pats = sorted(bags)
    y_sum = np.array([bags[p][1] for p in pats])
    p_sum = np.array([np.mean(oof[p]["hrd_sum_pred"]) for p in pats])
    from scipy.stats import pearsonr, spearmanr

    report = {"n_patients": len(pats), "in_dim": in_dim, "seeds": list(seeds),
              "hrd_sum_pearson": float(pearsonr(y_sum, p_sum)[0]), "hrd_sum_spearman": float(spearmanr(y_sum, p_sum)[0])}
    for i, t in enumerate(THRESHOLDS):
        y = np.array([bags[p][2][i] for p in pats]).astype(int)
        pb = np.array([np.mean(oof[p][f"ge{t}"]) for p in pats])
        auc, lo, hi = bootstrap_auc(y, pb)
        auc_r, lo_r, hi_r = bootstrap_auc(y, p_sum)
        report[f"ge{t}"] = {"prevalence": float(y.mean()), "auroc_binary_head": [auc, lo, hi],
                            "auroc_regression_head": [auc_r, lo_r, hi_r], "triage_binary_head": triage_points(y, pb)}
    with open(out_dir / "oof_predictions.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["patient", "hrd_sum", "hrd_sum_pred"] + [f"prob_ge{t}" for t in THRESHOLDS])
        for p, ys, ps in zip(pats, y_sum, p_sum):
            w.writerow([p, ys, round(float(ps), 2)] + [round(float(np.mean(oof[p][f"ge{t}"])), 4) for t in THRESHOLDS])
    json.dump(report, open(out_dir / "report.json", "w"), indent=2)
    return report
