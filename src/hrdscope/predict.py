"""Inference on a new slide: ensemble of cross-validated aggregators, calibrated probability,
attention heatmap and a JSON report. Research use only."""

from __future__ import annotations

import json
from pathlib import Path

import h5py
import numpy as np
import torch

from . import __version__
from .mil import GatedAttentionMIL, THRESHOLDS, predict


def load_ensemble(run_dir: Path, in_dim: int, device: str = "cpu") -> list[tuple[GatedAttentionMIL, float]]:
    params = {}
    pp = run_dir / "best_params.json"
    if pp.exists():
        params = json.load(open(pp))
    models = []
    for w in sorted(run_dir.glob("abmil_seed*_fold*.pt")):
        fold = w.stem.split("fold")[-1]
        cfg = params.get(fold, {}).get("params", {})
        m = GatedAttentionMIL(in_dim, cfg.get("hidden", 256), cfg.get("attn", 128), cfg.get("dropout", 0.25))
        m.load_state_dict(torch.load(w, map_location=device))
        m.to(device).eval()
        models.append((m, cfg.get("reg_scale", 50.0)))
    if not models:
        raise FileNotFoundError(f"no abmil_seed*_fold*.pt weights in {run_dir}")
    return models


def predict_features(h5_path: Path, run_dir: Path, device: str = "cpu") -> dict:
    with h5py.File(h5_path, "r") as h5:
        x = torch.from_numpy(h5["features"][:].astype(np.float32))
        coords = h5["coords"][:]
        attrs = {k: (v.item() if hasattr(v, "item") else v) for k, v in h5.attrs.items()}
    models = load_ensemble(run_dir, x.shape[1], device)
    outs = [predict(m, x, device, rs) for m, rs in models]
    probs = {f"ge{t}": [o["prob"][f"ge{t}"] for o in outs] for t in THRESHOLDS}
    att = np.mean([o["attention"] for o in outs], axis=0)
    hrd = [o["hrd_sum_pred"] for o in outs]
    return {
        "slide": attrs.get("slide"), "backbone": attrs.get("backbone"), "n_tiles": int(x.shape[0]),
        "hrdscope_version": __version__, "n_models": len(models),
        "hrd_sum": {"mean": float(np.mean(hrd)), "sd": float(np.std(hrd))},
        "probability": {k: {"mean": float(np.mean(v)), "sd": float(np.std(v)),
                             "ci90": [float(np.percentile(v, 5)), float(np.percentile(v, 95))]} for k, v in probs.items()},
        "attention": {"coords": coords.tolist(), "weights": att.tolist(), "tile_size0": attrs.get("tile_size0")},
        "disclaimer": "Research use only. Not a medical device. Not for clinical decision making.",
    }


def heatmap_png(result: dict, slide_path: Path, out_png: Path, target_mpp: float = 8.0) -> None:
    from PIL import Image

    from .wsi import open_slide, slide_mpp

    slide = open_slide(slide_path)
    ds = target_mpp / slide_mpp(slide)
    level = slide.get_best_level_for_downsample(ds)
    real = slide.level_downsamples[level]
    w, h = slide.level_dimensions[level]
    img = np.asarray(slide.read_region((0, 0), level, (w, h)).convert("RGB")).astype(np.float32)
    att = np.array(result["attention"]["weights"])
    att = (att - att.min()) / (att.max() - att.min() + 1e-9)
    size = max(1, int(result["attention"]["tile_size0"] / real))
    overlay = np.zeros((h, w), dtype=np.float32)
    for (x, y), a in zip(result["attention"]["coords"], att):
        xi, yi = int(x / real), int(y / real)
        overlay[yi:yi + size, xi:xi + size] = a
    colour = np.stack([overlay * 255, np.zeros_like(overlay), (1 - overlay) * 255], -1)
    alpha = (overlay > 0)[..., None] * 0.45
    out = img * (1 - alpha) + colour * alpha
    Image.fromarray(out.astype(np.uint8)).save(out_png)
