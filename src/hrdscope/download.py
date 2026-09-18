"""Parallel, resumable download of a GDC open-access file using HTTP range requests."""

from __future__ import annotations

import concurrent.futures as cf
import hashlib
import os
import time
import urllib.request
from pathlib import Path

GDC_DATA = "https://api.gdc.cancer.gov/data/"
IDC_HOST = "idc-open-data.s3.amazonaws.com"
_IDC_CLIENT = None


def download_idc(series_uid: str, dest_dir: Path) -> Path:
    """Download one DICOM whole-slide series from NCI Imaging Data Commons (public S3, no credentials).

    Returns the path of one instance; OpenSlide opens the whole pyramid from the files in its directory."""
    global _IDC_CLIENT
    from idc_index import IDCClient

    if _IDC_CLIENT is None:
        _IDC_CLIENT = IDCClient()
    dest_dir.mkdir(parents=True, exist_ok=True)
    _IDC_CLIENT.download_from_selection(seriesInstanceUID=[series_uid], downloadDir=str(dest_dir), dirTemplate=None,
                                        quiet=True, show_progress_bar=False)
    files = sorted(p for p in dest_dir.rglob("*.dcm"))
    if not files:
        raise OSError(f"IDC returned no DICOM files for {series_uid}")
    return files[0]


def _size(url: str, retries: int = 6) -> int:
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"Range": "bytes=0-0"})
            with urllib.request.urlopen(req, timeout=120) as r:
                return int(r.headers.get("Content-Range", "").split("/")[-1])
        except Exception:
            if attempt == retries - 1:
                raise
            time.sleep(2 ** attempt)
    raise OSError("unreachable")


def _fetch_range(url: str, start: int, end: int, retries: int = 6) -> bytes:
    """GDC drops connections under load; back off exponentially (1, 2, 4 ... s) and retry."""
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"Range": f"bytes={start}-{end}"})
            with urllib.request.urlopen(req, timeout=300) as r:
                data = r.read()
            if len(data) == end - start + 1:
                return data
        except Exception:
            if attempt == retries - 1:
                raise
        time.sleep(2 ** attempt)
    raise OSError(f"short read for bytes {start}-{end}")


def download_gdc(file_id: str, dest: Path, expected_size: int | None = None, md5: str | None = None,
                 workers: int = 8, chunk: int = 16 << 20) -> Path:
    return download_url(GDC_DATA + file_id, dest, expected_size, md5, workers, chunk)


def download_url(url: str, dest: Path, expected_size: int | None = None, md5: str | None = None,
                 workers: int = 8, chunk: int = 16 << 20) -> Path:
    """Parallel range download of any HTTP(S) URL that supports byte ranges (GDC, TCIA PathDB)."""
    file_id = url.rsplit("/", 1)[-1]
    if dest.exists() and expected_size and dest.stat().st_size == expected_size:
        return dest
    size = expected_size or _size(url)
    part = dest.with_suffix(dest.suffix + ".part")
    with open(part, "wb") as fh:
        fh.truncate(size)
    fd = os.open(part, os.O_WRONLY)
    try:
        ranges = [(s, min(s + chunk, size) - 1) for s in range(0, size, chunk)]

        def job(rng):
            data = _fetch_range(url, *rng)
            os.pwrite(fd, data, rng[0])
            return len(data)

        with cf.ThreadPoolExecutor(workers) as ex:
            got = sum(ex.map(job, ranges))
    finally:
        os.close(fd)
    if got != size:
        raise OSError(f"downloaded {got} of {size} bytes for {file_id}")
    if md5:
        h = hashlib.md5()
        with open(part, "rb") as fh:
            for block in iter(lambda: fh.read(1 << 24), b""):
                h.update(block)
        if h.hexdigest() != md5:
            part.unlink()
            raise OSError(f"md5 mismatch for {file_id}")
    part.rename(dest)
    return dest
