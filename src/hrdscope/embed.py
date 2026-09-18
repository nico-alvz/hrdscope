"""Tile embeddings with open-weight pathology foundation models.

Backbones (all permissively licensed):
* ``midnight``  kaiko-ai/midnight (Midnight-12k, MIT, ViT-g/14, 1536-d)
* ``hibou-b``   histai/hibou-b (Apache-2.0, ViT-B/14, 768-d, gated download)
* ``h-optimus`` bioptimus/H-optimus-0 (Apache-2.0, ViT-g/14, 1536-d, gated download)

Embeddings are written per slide to HDF5: ``features`` (n, d) float16, ``coords`` (n, 2) int32
level-0 top-left coordinates, ``tissue_frac`` (n,), plus slide metadata as attributes.
"""

from __future__ import annotations

import time
from pathlib import Path

import h5py
import numpy as np
import torch
from PIL import Image

from . import __version__
from .wsi import grid_tiles, open_slide, read_tile, slide_mpp

BACKBONES = {
    "midnight": {"hf": "kaiko-ai/midnight", "dim": 1536, "mean": (0.5, 0.5, 0.5), "std": (0.5, 0.5, 0.5)},
    "hibou-b": {"hf": "histai/hibou-b", "dim": 768, "mean": (0.7068, 0.5755, 0.7220), "std": (0.1950, 0.2316, 0.1816)},
    "h-optimus": {"hf": "bioptimus/H-optimus-0", "dim": 1536, "mean": (0.7076, 0.5814, 0.7038), "std": (0.2118, 0.2305, 0.1855)},
}


class Backbone:
    def __init__(self, name: str = "midnight", device: str | None = None, dtype: torch.dtype = torch.float16):
        if name not in BACKBONES:
            raise ValueError(f"unknown backbone {name}; choose from {list(BACKBONES)}")
        self.name = name
        self.cfg = BACKBONES[name]
        self.device = device or ("cuda" if torch.cuda.is_available() else "cpu")
        self.dtype = dtype if self.device == "cuda" else torch.float32
        if name == "h-optimus":
            import timm
            self.model = timm.create_model("hf_hub:bioptimus/H-optimus-0", pretrained=True, init_values=1e-5, dynamic_img_size=False)
            self._forward = self._forward_timm
        else:
            from transformers import AutoModel
            from transformers.utils import logging as hf_logging

            hf_logging.disable_progress_bar()
            hf_logging.set_verbosity_error()
            # load straight into the target dtype so that no fp32 copy of the 1.1 B parameters is made
            self.model = AutoModel.from_pretrained(self.cfg["hf"], trust_remote_code=(name == "hibou-b"), dtype=self.dtype)
            self._forward = self._forward_hf
        self.model.eval().to(self.device, self.dtype)
        if self.device == "cpu":
            torch.set_num_threads(max(1, torch.get_num_threads()))
        self.mean = torch.tensor(self.cfg["mean"], device=self.device, dtype=self.dtype).view(1, 3, 1, 1)
        self.std = torch.tensor(self.cfg["std"], device=self.device, dtype=self.dtype).view(1, 3, 1, 1)

    @property
    def dim(self) -> int:
        return self.cfg["dim"]

    @property
    def out_dim(self) -> int:
        return 2 * self.dim if self.name == "midnight" else self.dim

    def _forward_hf(self, x: torch.Tensor) -> torch.Tensor:
        out = self.model(pixel_values=x).last_hidden_state
        if self.name == "midnight":
            # Midnight recommends concatenating CLS with the mean of patch tokens
            return torch.cat([out[:, 0], out[:, 1:].mean(dim=1)], dim=-1)
        return out[:, 0]

    def _forward_timm(self, x: torch.Tensor) -> torch.Tensor:
        return self.model(x)

    @torch.inference_mode()
    def __call__(self, images: list[Image.Image]) -> np.ndarray:
        arr = np.stack([np.asarray(im, dtype=np.uint8) for im in images])  # (b, h, w, 3)
        x = torch.from_numpy(arr).to(self.device).permute(0, 3, 1, 2).to(self.dtype) / 255.0
        x = (x - self.mean) / self.std
        return self._forward(x).float().cpu().numpy()


def embed_slide(slide_path: str | Path, out_path: str | Path, backbone: Backbone, tile_px: int = 224,
                tile_mpp: float = 0.5, batch_size: int = 32, max_tiles: int | None = None,
                mpp: float | None = None, seed: int = 0) -> dict:
    slide_path, out_path = Path(slide_path), Path(out_path)
    t0 = time.time()
    slide = open_slide(slide_path)
    tiles, size0 = grid_tiles(slide, tile_px=tile_px, tile_mpp=tile_mpp, mpp=mpp)
    n_total = len(tiles)
    if max_tiles and n_total > max_tiles:
        rng = np.random.default_rng(seed)
        tiles = [tiles[i] for i in sorted(rng.choice(n_total, max_tiles, replace=False))]
    feats = np.zeros((len(tiles), backbone.out_dim), dtype=np.float16)
    coords = np.array([(t.x, t.y) for t in tiles], dtype=np.int32).reshape(-1, 2)
    fracs = np.array([t.tissue_frac for t in tiles], dtype=np.float32)
    for i in range(0, len(tiles), batch_size):
        batch = [read_tile(slide, t, tile_px) for t in tiles[i : i + batch_size]]
        feats[i : i + len(batch)] = backbone(batch).astype(np.float16)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(out_path, "w") as h5:
        h5.create_dataset("features", data=feats, compression="gzip")
        h5.create_dataset("coords", data=coords)
        h5.create_dataset("tissue_frac", data=fracs)
        h5.attrs.update({
            "slide": slide_path.name, "backbone": backbone.name, "hf_id": backbone.cfg["hf"],
            "tile_px": tile_px, "tile_mpp": tile_mpp, "tile_size0": size0, "mpp0": mpp or slide_mpp(slide),
            "n_tiles_total": n_total, "n_tiles": len(tiles), "hrdscope_version": __version__,
            "seconds": round(time.time() - t0, 1),
        })
    slide.close()
    return {"slide": slide_path.name, "n_tiles_total": n_total, "n_tiles": len(tiles), "seconds": round(time.time() - t0, 1)}
