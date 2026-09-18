"""Write data/manifests/tcga_ov_slides_idc.tsv from NCI Imaging Data Commons (DICOM copies of the TCGA-OV slides)."""
import csv
from pathlib import Path

from idc_index import IDCClient

OUT = Path(__file__).resolve().parents[1] / "data" / "manifests" / "tcga_ov_slides_idc.tsv"


def main():
    client = IDCClient()
    client.fetch_index("sm_index")
    version = client.get_idc_version()
    df = client.sql_query("""
        SELECT i.SeriesInstanceUID, i.PatientID, i.series_size_MB, i.series_aws_url, s.ContainerIdentifier,
               s.min_PixelSpacing_2sf, i.license_short_name
        FROM index i JOIN sm_index s USING (SeriesInstanceUID)
        WHERE i.collection_id = 'tcga_ov' ORDER BY s.ContainerIdentifier""")
    with open(OUT, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["series_uid", "file_name", "file_size", "strategy", "patient", "sample_type", "pixel_spacing_mm",
                    "license", "aws_url", "idc_version"])
        for r in df.itertuples():
            barcode = r.ContainerIdentifier
            strategy = "Diagnostic Slide" if barcode.split("-")[-1].startswith("DX") else "Tissue Slide"
            w.writerow([r.SeriesInstanceUID, barcode + ".dcm", int(r.series_size_MB * 1e6), strategy, r.PatientID,
                        barcode[13:15], r.min_PixelSpacing_2sf, r.license_short_name, r.series_aws_url, version])
    print(f"{len(df)} series from IDC {version} -> {OUT}")


if __name__ == "__main__":
    main()
