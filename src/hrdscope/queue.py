"""Offline-friendly producer/consumer queue for slide embedding.

``fetch`` downloads slides into a local cache while the network is up, keeping a
minimum amount of free disk space, and marks each completed slide with a JSON
file. ``work`` embeds whatever is in the cache and deletes each slide after it is
embedded; it never touches the network, so the GPU keeps working overnight
without Wi-Fi as long as the cache holds slides.
"""

from __future__ import annotations

import csv
import json
import shutil
import socket
import time
from pathlib import Path

from .download import IDC_HOST, download_idc, download_url

GDC_HOST = "api.gdc.cancer.gov"


def is_tumour_sample(name: str) -> bool:
    """TCGA sample-type codes 01-09 are tumour, 10-19 normal, 20-29 control; non-TCGA names are kept."""
    if not name.startswith("TCGA-") or len(name) < 15 or not name[13:15].isdigit():
        return True
    return int(name[13:15]) < 10


def online(host: str) -> bool:
    try:
        socket.create_connection((host, 443), timeout=10).close()
        return True
    except OSError:
        return False


def plan_rows(plan: list[dict], backbone: str) -> list[dict]:
    """Expand a plan (list of jobs) into manifest rows annotated with job priority and output directory."""
    rows: list[dict] = []
    for priority, job in enumerate(plan):
        with open(job["manifest"], newline="") as fh:
            items = list(csv.DictReader(fh, delimiter="\t"))
        if job.get("strategy"):
            items = [r for r in items if r.get("strategy", "").startswith(job["strategy"])]
        if job.get("patients"):
            keep = set(Path(job["patients"]).read_text().split())
            items = [r for r in items if r.get("patient") in keep]
        if not job.get("include_normal", False):
            items = [r for r in items if is_tumour_sample(r["file_name"])]
        items.sort(key=lambda r: int(r.get("file_size", 0) or 0))
        for r in items:
            r["_priority"] = priority
            r["_out"] = str(Path(job["out_dir"]) / backbone / (Path(r["file_name"]).stem + ".h5"))
            rows.append(r)
    return rows


def fetch(plan: list[dict], cache: Path, backbone: str = "midnight", min_free_gb: float = 12.0,
          workers: int = 6, log=print) -> None:
    cache.mkdir(parents=True, exist_ok=True)
    (cache / "FETCH_DONE").unlink(missing_ok=True)
    for partial in cache.glob("*.partial"):
        shutil.rmtree(partial, ignore_errors=True)
    for r in plan_rows(plan, backbone):
        stem = Path(r["file_name"]).stem
        if Path(r["_out"]).exists() or (cache / f"{stem}.json").exists() or (cache / f"{stem}.failed.json").exists():
            continue
        size = int(r.get("file_size", 0) or 0)
        while shutil.disk_usage(cache).free - size < min_free_gb * 1e9:
            time.sleep(60)  # wait for the worker to free space
        host = IDC_HOST if r.get("series_uid") else GDC_HOST
        partial = cache / f"{stem}.partial"
        for attempt in range(4):
            while not online(host):
                time.sleep(60)
            try:
                shutil.rmtree(partial, ignore_errors=True)
                if r.get("series_uid"):
                    download_idc(r["series_uid"], partial)
                else:
                    url = r.get("url") or "https://api.gdc.cancer.gov/data/" + r["file_id"]
                    download_url(url, partial / r["file_name"], size or None, r.get("md5") or None, workers=workers)
                break
            except Exception as exc:
                if attempt == 3:
                    log(json.dumps({"fetch_error": r["file_name"], "error": repr(exc)}))
                    shutil.rmtree(partial, ignore_errors=True)
                    partial = None
                    break
                time.sleep(30 * (attempt + 1))
        if partial is None:
            continue
        partial.rename(cache / stem)
        marker = {k: v for k, v in r.items()}
        (cache / f"{stem}.json").write_text(json.dumps(marker))
        log(json.dumps({"fetched": r["file_name"], "gb": round(size / 1e9, 2),
                        "free_gb": round(shutil.disk_usage(cache).free / 1e9, 1)}))
    (cache / "FETCH_DONE").write_text(time.strftime("%Y-%m-%dT%H:%M:%S"))


def _slide_path(item_dir: Path, file_name: str) -> Path:
    direct = item_dir / file_name
    if direct.exists():
        return direct
    dcm = sorted(item_dir.rglob("*.dcm"))
    if dcm:
        return dcm[0]
    raise FileNotFoundError(f"no slide inside {item_dir}")


def _pending(cache: Path) -> list[tuple[dict, Path]]:
    items = [(json.loads(m.read_text()), m) for m in cache.glob("*.json") if not m.name.endswith(".failed.json")]
    items.sort(key=lambda t: (t[0]["_priority"], int(t[0].get("file_size", 0) or 0)))
    return items


def work(cache: Path, staging: Path, backbone_name: str = "midnight", device: str | None = None, batch_size: int = 32,
         tile_px: int = 224, tile_mpp: float = 0.5, max_tiles: int | None = None, log_path: Path | None = None,
         wait: bool = True) -> None:
    """Embed cached slides. Each slide is first copied from the (possibly slow) cache disk to ``staging`` on a fast
    disk; the copy of the next slide overlaps with embedding of the current one."""
    import concurrent.futures as cf

    from .embed import Backbone, embed_slide

    bb = Backbone(backbone_name, device=device)
    staging.mkdir(parents=True, exist_ok=True)
    log = open(log_path, "a") if log_path else None
    pool = cf.ThreadPoolExecutor(1)

    def stage(r: dict) -> Path:
        stem = Path(r["file_name"]).stem
        dst = staging / stem
        if not dst.exists():
            shutil.copytree(cache / stem, staging / f"{stem}.copying")
            (staging / f"{stem}.copying").rename(dst)
        return dst

    staged: dict[str, cf.Future] = {}
    while True:
        pending = _pending(cache)
        if not pending:
            if not wait or (cache / "FETCH_DONE").exists():
                break
            time.sleep(30)
            continue
        r, marker = pending[0]
        stem = Path(r["file_name"]).stem
        if stem not in staged:
            staged[stem] = pool.submit(stage, r)
        if len(pending) > 1:
            nxt = Path(pending[1][0]["file_name"]).stem
            if nxt not in staged:
                staged[nxt] = pool.submit(stage, pending[1][0])
        try:
            local = staged.pop(stem).result()
            info = embed_slide(_slide_path(local, r["file_name"]), Path(r["_out"]), bb, tile_px=tile_px, tile_mpp=tile_mpp,
                               batch_size=batch_size, max_tiles=max_tiles, slide_name=r["file_name"])
            info.update(patient=r.get("patient"), out=r["_out"], queued=len(pending) - 1)
            marker.unlink()
        except Exception as exc:
            info = {"slide": r["file_name"], "error": repr(exc), "queued": len(pending) - 1}
            marker.rename(cache / f"{stem}.failed.json")
        shutil.rmtree(staging / stem, ignore_errors=True)
        shutil.rmtree(cache / stem, ignore_errors=True)
        line = json.dumps(info)
        print(line, flush=True)
        if log:
            log.write(line + "\n")
            log.flush()
    pool.shutdown()
