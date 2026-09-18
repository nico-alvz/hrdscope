"""Regenerate data/manifests/*.tsv from the public GDC API (no token needed)."""
import json
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "manifests"


def query(filters, fields, size=5000):
    body = json.dumps({"filters": filters, "fields": fields, "size": size, "format": "JSON"}).encode()
    req = urllib.request.Request("https://api.gdc.cancer.gov/files", data=body, headers={"Content-Type": "application/json"})
    return json.load(urllib.request.urlopen(req, timeout=120))["data"]["hits"]


def in_(field, values):
    return {"op": "in", "content": {"field": field, "value": values}}


def main():
    ov = in_("cases.project.project_id", ["TCGA-OV"])
    slides = query({"op": "and", "content": [ov, in_("data_type", ["Slide Image"])]},
                   "file_id,file_name,file_size,experimental_strategy,cases.submitter_id,md5sum")
    with open(OUT / "tcga_ov_slides.tsv", "w") as fo:
        fo.write("file_id\tfile_name\tfile_size\tstrategy\tpatient\tmd5\n")
        for h in sorted(slides, key=lambda h: h["file_name"]):
            fo.write(f"{h['file_id']}\t{h['file_name']}\t{h['file_size']}\t{h['experimental_strategy']}\t{h['cases'][0]['submitter_id']}\t{h['md5sum']}\n")
    ascn = query({"op": "and", "content": [ov, in_("data_type", ["Allele-specific Copy Number Segment"]), in_("access", ["open"]),
                                           in_("analysis.workflow_type", ["ASCAT3"])]},
                 "file_id,file_name,cases.submitter_id,cases.samples.sample_type,analysis.workflow_type")
    with open(OUT / "tcga_ov_ascat3.tsv", "w") as fo:
        fo.write("file_id\tfile_name\tpatient\tsample_type\tworkflow\n")
        for h in sorted(ascn, key=lambda h: h["file_name"]):
            fo.write(f"{h['file_id']}\t{h['file_name']}\t{h['cases'][0]['submitter_id']}\t{h['cases'][0]['samples'][0]['sample_type']}\t{h['analysis']['workflow_type']}\n")
    print(len(slides), "slides;", len(ascn), "ASCAT3 profiles")


if __name__ == "__main__":
    main()
