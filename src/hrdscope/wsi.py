"""Whole-slide image access: tissue detection and tiling at a fixed physical resolution."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from PIL import Image

try:
    import openslide
except ImportError as exc:  # pragma: no cover
    raise ImportError("install hrdscope[slides] for whole-slide image support") from exc


@dataclass(frozen=True)
class Tile:
    x: int  # level-0 coordinates of the top-left corner
    y: int
    size0: int  # side length in level-0 pixels
    tissue_frac: float


def slide_mpp(slide: "openslide.OpenSlide") -> float:
    p = slide.properties
    for key in (openslide.PROPERTY_NAME_MPP_X, "aperio.MPP", "openslide.mpp-x"):
        if key in p:
            try:
                return float(p[key])
            except ValueError:
                pass
    mag = p.get(openslide.PROPERTY_NAME_OBJECTIVE_POWER)
    if mag:
        return 10.0 / float(mag)  # 40x -> 0.25, 20x -> 0.5
    raise ValueError("cannot determine microns per pixel; pass --mpp explicitly")


def tissue_mask(slide: "openslide.OpenSlide", target_mpp: float = 16.0,
                sat_threshold: int = 20, min_value: int = 30, max_value: int = 230) -> tuple[np.ndarray, float]:
    """Boolean tissue mask on a low-resolution thumbnail plus its downsample factor."""
    mpp0 = slide_mpp(slide)
    ds = target_mpp / mpp0
    level = slide.get_best_level_for_downsample(ds)
    w, h = slide.level_dimensions[level]
    img = slide.read_region((0, 0), level, (w, h)).convert("RGB")
    real_ds = slide.level_downsamples[level]
    hsv = np.asarray(img.convert("HSV"))
    rgb = np.asarray(img)
    sat, val = hsv[..., 1], hsv[..., 2]
    mask = (sat > sat_threshold) & (val > min_value) & (val < max_value)
    # drop pen marks and black artefacts: strongly coloured pixels with low red/green balance
    r, g, b = rgb[..., 0].astype(int), rgb[..., 1].astype(int), rgb[..., 2].astype(int)
    pen = ((b - r) > 40) | ((g - r) > 40)
    mask &= ~pen
    return mask, real_ds


def grid_tiles(slide: "openslide.OpenSlide", tile_px: int = 224, tile_mpp: float = 0.5,
               min_tissue: float = 0.5, mpp: float | None = None) -> tuple[list[Tile], int]:
    """Enumerate non-overlapping tiles covering tissue at ``tile_mpp`` microns per pixel."""
    mpp0 = mpp or slide_mpp(slide)
    size0 = int(round(tile_px * tile_mpp / mpp0))
    mask, ds = tissue_mask(slide)
    W, H = slide.dimensions
    tiles: list[Tile] = []
    for y in range(0, H - size0 + 1, size0):
        for x in range(0, W - size0 + 1, size0):
            mx0, my0 = int(x / ds), int(y / ds)
            mx1, my1 = max(mx0 + 1, int((x + size0) / ds)), max(my0 + 1, int((y + size0) / ds))
            frac = float(mask[my0:my1, mx0:mx1].mean()) if my1 <= mask.shape[0] and mx1 <= mask.shape[1] else 0.0
            if frac >= min_tissue:
                tiles.append(Tile(x, y, size0, frac))
    return tiles, size0


def read_tile(slide: "openslide.OpenSlide", tile: Tile, tile_px: int = 224) -> Image.Image:
    """Read a tile at the best pyramid level and resample to ``tile_px`` square."""
    level = slide.get_best_level_for_downsample(tile.size0 / tile_px)
    ds = slide.level_downsamples[level]
    size_l = int(round(tile.size0 / ds))
    img = slide.read_region((tile.x, tile.y), level, (size_l, size_l)).convert("RGB")
    if size_l != tile_px:
        img = img.resize((tile_px, tile_px), Image.BILINEAR)
    return img


def open_slide(path: str | Path) -> "openslide.OpenSlide":
    return openslide.OpenSlide(str(path))
