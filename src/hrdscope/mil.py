"""Attention-based multiple-instance learning on tile embeddings.

One bag = all embedded tiles of one patient (all their slides). The model
outputs a continuous HRD-sum estimate and a logit per binary threshold.
"""

from __future__ import annotations

import csv
from pathlib import Path

import h5py
import numpy as np
import torch
from torch import nn

THRESHOLDS = (33, 42, 63)


class GatedAttentionMIL(nn.Module):
    def __init__(self, in_dim: int, hidden: int = 256, attn: int = 128, dropout: float = 0.25, n_bin: int = 3):
        super().__init__()
        self.proj = nn.Sequential(nn.Linear(in_dim, hidden), nn.GELU(), nn.Dropout(dropout))
        self.att_v = nn.Sequential(nn.Linear(hidden, attn), nn.Tanh())
        self.att_u = nn.Sequential(nn.Linear(hidden, attn), nn.Sigmoid())
        self.att_w = nn.Linear(attn, 1)
        self.head_reg = nn.Linear(hidden, 1)
        self.head_bin = nn.Linear(hidden, n_bin)

    def forward(self, x: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        h = self.proj(x)  # (n, hidden)
        a = self.att_w(self.att_v(h) * self.att_u(h)).squeeze(-1)  # (n,)
        a = torch.softmax(a, dim=0)
        z = (a.unsqueeze(-1) * h).sum(0)  # (hidden,)
        return self.head_reg(z).squeeze(-1), self.head_bin(z), a


def load_bag(paths: list[Path], max_tiles: int | None = None, rng: np.random.Generator | None = None) -> torch.Tensor:
    feats = []
    for p in paths:
        with h5py.File(p, "r") as h5:
            feats.append(h5["features"][:].astype(np.float32))
    x = np.concatenate(feats, 0)
    if max_tiles and len(x) > max_tiles:
        idx = (rng or np.random.default_rng(0)).choice(len(x), max_tiles, replace=False)
        x = x[np.sort(idx)]
    return torch.from_numpy(x)


def index_features(feature_dir: Path) -> dict[str, list[Path]]:
    """Map patient id (first 12 chars of a TCGA slide name) to its HDF5 files."""
    out: dict[str, list[Path]] = {}
    for p in sorted(feature_dir.glob("*.h5")):
        out.setdefault(p.name[:12], []).append(p)
    return out


def read_labels(label_path: Path) -> dict[str, dict]:
    return {r["patient"]: r for r in csv.DictReader(open(label_path, newline=""), delimiter="\t")}


def train_one(model: GatedAttentionMIL, bags: list[tuple[torch.Tensor, float, np.ndarray]], val: list, epochs: int = 40,
              lr: float = 2e-4, weight_decay: float = 1e-2, patience: int = 8, device: str = "cpu",
              reg_scale: float = 50.0, seed: int = 0) -> tuple[GatedAttentionMIL, dict]:
    torch.manual_seed(seed)
    model.to(device)
    opt = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=weight_decay)
    bce = nn.BCEWithLogitsLoss()
    best, best_state, bad = float("inf"), None, 0
    rng = np.random.default_rng(seed)
    for ep in range(epochs):
        model.train()
        order = rng.permutation(len(bags))
        for i in order:
            x, y_reg, y_bin = bags[i]
            x = x.to(device)
            pr, pb, _ = model(x)
            loss = ((pr - y_reg / reg_scale) ** 2) + bce(pb, torch.tensor(y_bin, dtype=torch.float32, device=device))
            opt.zero_grad(); loss.backward(); opt.step()
        vl = evaluate_loss(model, val, device, reg_scale)
        if vl < best - 1e-4:
            best, bad = vl, 0
            best_state = {k: v.detach().clone() for k, v in model.state_dict().items()}
        else:
            bad += 1
            if bad >= patience:
                break
    if best_state:
        model.load_state_dict(best_state)
    return model, {"best_val_loss": best, "epochs": ep + 1}


@torch.no_grad()
def evaluate_loss(model, bags, device, reg_scale):
    model.eval()
    bce = nn.BCEWithLogitsLoss()
    tot = 0.0
    for x, y_reg, y_bin in bags:
        pr, pb, _ = model(x.to(device))
        tot += float(((pr - y_reg / reg_scale) ** 2) + bce(pb, torch.tensor(y_bin, dtype=torch.float32, device=device)))
    return tot / max(1, len(bags))


@torch.no_grad()
def predict(model, x: torch.Tensor, device: str = "cpu", reg_scale: float = 50.0) -> dict:
    model.eval()
    pr, pb, a = model(x.to(device))
    return {"hrd_sum_pred": float(pr) * reg_scale,
            "prob": {f"ge{t}": float(torch.sigmoid(pb[i])) for i, t in enumerate(THRESHOLDS)},
            "attention": a.cpu().numpy()}
