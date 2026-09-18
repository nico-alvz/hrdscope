# Data acquisition

Everything below is open access. Nothing requires dbGaP or a data-use agreement.

## TCGA-OV slides and copy-number (NCI Genomic Data Commons)

* Manifests in `data/manifests/` were generated from the GDC API on 2026-09-18 (`scripts/` reproduces them): 107 diagnostic FFPE slides (106 patients, 209 GB), 1,374 frozen tissue slides (589 patients, 260 GB), 555 ASCAT3 allele-specific copy-number profiles.
* `hrdscope labels` downloads the 555 ASCAT3 files (6.5 MB) and writes `data/labels/tcga_ov_hrd.tsv`.
* `hrdscope stream --strategy Diagnostic --patients data/splits/tcga_ov_dx_patients.txt` downloads each slide with 8 parallel range requests (about 15 MB/s from GDC), verifies the md5, embeds it and deletes it. Use `--strategy Tissue` for frozen slides.
* Citation: The Cancer Genome Atlas Research Network, Integrated genomic analyses of ovarian carcinoma, Nature 2011. GDC data are open under the NIH Genomic Data Sharing policy.

## Knijnenburg et al. 2018 HRD scores (validation of our scar scores)

* `TCGA_DDR_Data_Resources.zip` from https://gdc.cancer.gov/about-data/publications/PanCan-DDR-2018 (open). `DDRscores.tsv` rows match `Samples.tsv`; columns match `Scores.tsv` (HRD_TAI, HRD_LST, HRD_LOH, HRD_Score). The ovarian subset is committed as `data/manifests/knijnenburg2018_ov_hrd.tsv`.
* Citation: Knijnenburg TA et al., Genomic and molecular landscape of DNA damage repair deficiency across The Cancer Genome Atlas, Cell Reports 2018.

## PTRC-HGSOC slides and clinical data (The Cancer Imaging Archive), external test set

* 158 patients, 348 FFPE H&E slides (SVS, 20x, 0.5 micron/px), 120 GB, license CC BY 4.0, DOI 10.7937/6RDA-P940.
* Clinical table (platinum refractory vs sensitive, stage, grade, site, neo-adjuvant): `PTRC-HGSOC_List_clincal_data.xlsx` from the collection page, one row per slide with `File Name`, `Patient ID`, `Tumor response`.
* Slides are served directly by TCIA PathDB over HTTPS with byte-range support (about 16 MB/s per stream). `data/manifests/ptrc_hgsoc_slides.tsv` (348 rows, generated with `tcia_utils.pathdb.getImages("PTRC-HGSOC")` and joined with the clinical table) has the URL, patient, platinum response (91 sensitive, 67 refractory patients), site (FHCRC 156, Mayo 100, UAB 92 slides), stage, grade and pixel size (0.4965 micron/px).

```bash
hrdscope stream --manifest data/manifests/ptrc_hgsoc_slides.tsv --out-dir data/features_ptrc
```

* The Aspera Faspex link on the collection page also works with the IBM `ascli` client but is not needed.

* Scar metrics (nTAI, nLST, nLOH, wGII) for the discovery cohort are in Chowdhury et al., Cell 2023, Table S1 (open access on PMC, PMC10414761); download the supplementary xlsx manually from the PMC page because the server blocks scripted downloads.
* Citation: Chowdhury S et al., Proteogenomic analysis of chemo-refractory high-grade serous ovarian cancer, Cell 2023; plus the TCIA citation (Clark K et al., J Digit Imaging 2013).

## Not used

* CPTAC-OV (TCIA): the same patients as TCGA-OV, so it cannot serve as an external cohort.
* UBC-OCEAN (Kaggle): histotype labels only; reserved for a later HGSOC gating module.

## Foundation-model weights

* Midnight-12k: `kaiko-ai/midnight` on Hugging Face, MIT, no gating. Downloaded automatically on first use (4.4 GB).
* Hibou-B (`histai/hibou-b`) and H-optimus-0 (`bioptimus/H-optimus-0`): Apache-2.0, gated; run `huggingface-cli login` after accepting the terms on the model page.
