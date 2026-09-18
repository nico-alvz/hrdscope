# hrdscope

Open, reproducible homologous recombination deficiency (HRD) triage from routine H&E whole-slide images in ovarian cancer. MIT licensed, built only on open data and open-weight models.

About half of high-grade serous ovarian carcinomas are HRD-positive and respond to PARP inhibitors, but the genomic assays that establish HRD status are expensive and unavailable in most of the world. `hrdscope` aims to turn the H&E slide every patient already has into a calibrated triage signal, so that scarce sequencing can be directed to the patients who need it.

Status: early development. See `docs/PLAN.md` for the research plan and `docs/PRIOR_ART.md` for what already exists and why this is different.

## What works today

* `hrdscope labels` downloads the open-access GDC ASCAT3 allele-specific copy-number profiles for TCGA-OV and computes the three genomic scar scores (HRD-LOH, NtAI, LST) and HRD-sum for 555 patients, with binary labels at the 33 / 42 / 63 thresholds. On the 168 patients that overlap with the published Knijnenburg et al. 2018 table the Pearson correlation is 0.96.
* `hrdscope score file.seg.txt` scores any allele-specific segment file in GDC format.

## Install

```bash
pip install -e ".[dev]"          # labels and scoring only
pip install -e ".[slides,dev]"   # plus whole-slide image pipeline
pytest
hrdscope labels                  # writes data/labels/tcga_ov_hrd.tsv
```

## Data sources

* TCGA-OV slides and ASCAT3 segments: NCI Genomic Data Commons, open access.
* Knijnenburg et al. 2018, Cell Reports, TCGA DDR data resource (HRD scores used for validation).
* PTRC-HGSOC slides and clinical data: The Cancer Imaging Archive, CC BY 4.0.
* Chromosome arms: UCSC hg38 cytoBand.

## License

MIT. Research use only; this software is not a medical device and makes no clinical claims.
