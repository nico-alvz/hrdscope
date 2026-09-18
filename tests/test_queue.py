import csv
import json

import hrdscope.queue as q


def _manifest(path, rows):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def test_plan_rows_filters_and_orders(tmp_path):
    _manifest(tmp_path / "m.tsv", [
        {"file_name": "TCGA-01-0001-01A-01-TS1.dcm", "strategy": "Tissue Slide",
         "patient": "TCGA-01-0001", "file_size": 30, "series_uid": "a"},
        {"file_name": "TCGA-01-0001-11A-01-TS1.dcm", "strategy": "Tissue Slide",
         "patient": "TCGA-01-0001", "file_size": 10, "series_uid": "b"},
        {"file_name": "TCGA-01-0002-01Z-00-DX1.svs", "strategy": "Diagnostic Slide",
         "patient": "TCGA-01-0002", "file_size": 20, "series_uid": ""},
        {"file_name": "TCGA-01-0003-01A-01-BS1.dcm", "strategy": "Tissue Slide",
         "patient": "TCGA-01-0003", "file_size": 5, "series_uid": "c"},
    ])
    (tmp_path / "p.txt").write_text("TCGA-01-0001\nTCGA-01-0002\n")
    plan = [{"manifest": str(tmp_path / "m.tsv"), "strategy": "Tissue", "patients": str(tmp_path / "p.txt"), "out_dir": "o"}]
    rows = q.plan_rows(plan, "midnight")
    assert [r["file_name"] for r in rows] == ["TCGA-01-0001-01A-01-TS1.dcm"]
    assert rows[0]["_out"].endswith("o/midnight/TCGA-01-0001-01A-01-TS1.h5")


def test_fetch_then_work_offline(tmp_path, monkeypatch):
    _manifest(tmp_path / "m.tsv", [
        {"file_name": f"TCGA-01-000{i}-01Z-00-DX1.svs", "strategy": "Diagnostic Slide", "patient": f"TCGA-01-000{i}",
         "file_size": 10 * i, "file_id": f"id{i}", "md5": ""} for i in (1, 2)])
    plan = [{"manifest": str(tmp_path / "m.tsv"), "out_dir": str(tmp_path / "feat")}]
    monkeypatch.setattr(q, "online", lambda host: True)

    def fake_download(url, dest, *a, **k):
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(b"x")
        return dest

    monkeypatch.setattr(q, "download_url", fake_download)
    cache = tmp_path / "cache"
    q.fetch(plan, cache, min_free_gb=0, log=lambda line: None)
    assert sorted(p.name for p in cache.glob("*.json")) == ["TCGA-01-0001-01Z-00-DX1.json", "TCGA-01-0002-01Z-00-DX1.json"]
    assert (cache / "FETCH_DONE").exists()

    import hrdscope.embed as embed

    class FakeBackbone:
        name = "midnight"

        def __init__(self, *a, **k):
            pass

    def fake_embed(slide, out, bb, **k):
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text("h5")
        return {"slide": k["slide_name"], "n_tiles": 1}

    monkeypatch.setattr(embed, "Backbone", FakeBackbone)
    monkeypatch.setattr(embed, "embed_slide", fake_embed)
    q.work(cache, tmp_path / "staging", log_path=tmp_path / "work.jsonl")
    assert len(list((tmp_path / "feat" / "midnight").glob("*.h5"))) == 2
    assert not list(cache.glob("*.json")) and not list((tmp_path / "staging").iterdir())
    assert len((tmp_path / "work.jsonl").read_text().splitlines()) == 2
    json.loads((tmp_path / "work.jsonl").read_text().splitlines()[0])
