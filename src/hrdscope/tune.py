"""Nested hyper-parameter search with Optuna.

For each outer fold of the fixed split, Optuna searches the aggregator's
hyper-parameters using only the inner folds (train / inner-validation), so the
outer test fold and the external cohort never influence the choice. The best
configuration per outer fold is stored and reused by the benchmark.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
import torch

from .mil import (
    THRESHOLDS,
    GatedAttentionMIL,
    index_features,
    load_bag,
    predict,
    read_labels,
    train_one,
)

SPACE = {
    "hidden": ("cat", [128, 256, 512]),
    "attn": ("cat", [64, 128, 256]),
    "dropout": ("float", 0.0, 0.5),
    "lr": ("log", 1e-5, 1e-3),
    "weight_decay": ("log", 1e-4, 1e-1),
    "reg_scale": ("cat", [25.0, 50.0, 100.0]),
    "max_tiles": ("cat", [500, 1000, 2000, 4000]),
}


def _suggest(trial):
    cfg = {}
    for k, spec in SPACE.items():
        if spec[0] == "cat":
            cfg[k] = trial.suggest_categorical(k, spec[1])
        elif spec[0] == "float":
            cfg[k] = trial.suggest_float(k, spec[1], spec[2])
        else:
            cfg[k] = trial.suggest_float(k, spec[1], spec[2], log=True)
    return cfg


def _auroc(model, bags, device, reg_scale):
    from sklearn.metrics import roc_auc_score

    y = np.array([b[2][1] for b in bags])  # threshold 42 is the primary endpoint
    p = np.array([predict(model, b[0], device, reg_scale)["prob"]["ge42"] for b in bags])
    return roc_auc_score(y, p) if len(set(y)) == 2 else 0.5


def tune(feature_dir: Path, label_path: Path, split_path: Path, out_dir: Path, n_trials: int = 30,
         device: str | None = None, epochs: int = 40, seed: int = 0) -> dict:
    import optuna

    optuna.logging.set_verbosity(optuna.logging.WARNING)
    device = device or ("cuda" if torch.cuda.is_available() else "cpu")
    feats, labels = index_features(feature_dir), read_labels(label_path)
    splits = [r for r in csv.DictReader(open(split_path, newline=""), delimiter="\t") if r["patient"] in feats]
    folds = sorted({int(r["fold"]) for r in splits})
    rng = np.random.default_rng(seed)
    full = {r["patient"]: (load_bag(feats[r["patient"]], None, rng), float(labels[r["patient"]]["hrd_sum"]),
                           np.array([int(labels[r["patient"]][f"hrd_ge{t}"]) for t in THRESHOLDS], dtype=np.float32))
            for r in splits}
    in_dim = next(iter(full.values()))[0].shape[1]
    out_dir.mkdir(parents=True, exist_ok=True)
    best_per_fold = {}

    def subsample(p, n):
        x, yr, yb = full[p]
        if n and len(x) > n:
            x = x[torch.from_numpy(np.sort(rng.choice(len(x), n, replace=False)))]
        return x, yr, yb

    for outer in folds:
        inner_folds = [f for f in folds if f != outer]

        def objective(trial):
            cfg = _suggest(trial)
            scores = []
            for v in inner_folds:
                train = [r["patient"] for r in splits if int(r["fold"]) not in (outer, v)]
                val = [r["patient"] for r in splits if int(r["fold"]) == v]
                model = GatedAttentionMIL(in_dim, cfg["hidden"], cfg["attn"], cfg["dropout"])
                model, _ = train_one(model, [subsample(p, cfg["max_tiles"]) for p in train],
                                     [subsample(p, cfg["max_tiles"]) for p in val], epochs=epochs, lr=cfg["lr"],
                                     weight_decay=cfg["weight_decay"], device=device, reg_scale=cfg["reg_scale"],
                                     seed=seed)
                scores.append(_auroc(model, [full[p] for p in val], device, cfg["reg_scale"]))
                trial.report(float(np.mean(scores)), len(scores))
                if trial.should_prune():
                    raise optuna.TrialPruned()
            return float(np.mean(scores))

        study = optuna.create_study(direction="maximize", sampler=optuna.samplers.TPESampler(seed=seed),
                                    pruner=optuna.pruners.MedianPruner(n_warmup_steps=1))
        study.optimize(objective, n_trials=n_trials)
        best_per_fold[str(outer)] = {"params": study.best_params, "inner_auroc_ge42": study.best_value,
                                     "n_trials": len(study.trials)}
        print(json.dumps({"outer_fold": outer, **best_per_fold[str(outer)]}), flush=True)
        study.trials_dataframe().to_csv(out_dir / f"trials_outer{outer}.csv", index=False)
    json.dump(best_per_fold, open(out_dir / "best_params.json", "w"), indent=2)
    return best_per_fold
