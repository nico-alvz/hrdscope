# hrdscope: open HRD triage from routine H&E in ovarian cancer

Research plan, version 2026-09-18. All artefacts (code, labels, weights, embeddings, manuscript) are released under MIT / CC BY with a Zenodo DOI.

## 1. Problem

About half of high-grade serous ovarian carcinomas (HGSOC) carry homologous recombination deficiency (HRD) and benefit from PARP inhibitors. HRD status is determined by genomic assays (Myriad myChoice, scarHRD on WGS/WES, SNP arrays) that cost hundreds to thousands of dollars per patient and are unavailable in most of the world. Every ovarian cancer patient, however, already has an H&E slide. Several groups have shown that deep learning can recover part of the HRD signal from H&E (see `PRIOR_ART.md`), but none of that work is usable today: weights are withheld, licenses are academic-only, pipelines are notebooks that need 100+ GPU hours, labels are private, and thresholds differ between papers.

## 2. Goal

Ship a tool that a pathology lab or a researcher can install with `pip` or `docker`, point at an `.svs` file, and get a calibrated HRD probability with an uncertainty estimate, a tile-level attention map and a one-page report, within minutes on an 8 GB consumer GPU (or slower on CPU). Frame the output as **triage**, not diagnosis: which patients can safely skip a confirmatory assay (rule-out) and which should be prioritised for it (rule-in) when sequencing capacity is limited.

## 3. Contributions

C1. **Open label table for TCGA-OV.** A Python reimplementation of the three genomic scar scores (HRD-LOH, NtAI, LST) applied to GDC open-access ASCAT3 allele-specific copy-number segments gives HRD-sum for 555 TCGA-OV patients (versus 173 in the only public table). On the 168 overlapping patients the Pearson correlation with Knijnenburg et al. 2018 is 0.96 (LOH 0.97, LST 0.94, NtAI 0.82). Labels are provided at the three thresholds used in the literature (33, 42, 63) plus the continuous score. Done, see `src/hrdscope/scars.py` and `data/labels/`.

C2. **Open benchmark.** Fixed, patient-level, stratified 5-fold splits on TCGA-OV (FFPE diagnostic slides and frozen tissue slides handled as separate domains), a frozen evaluation protocol (AUROC, AUPRC, calibration, decision curves, triage operating points, bootstrap CIs) and an external test on PTRC-HGSOC (158 patients, CC BY 4.0). Anyone can plug in a new backbone or aggregator and compare.

C3. **Model.** Tile embeddings from an MIT-licensed pathology foundation model (Midnight-12k; Hibou-B and H-optimus-0 as Apache-2.0 alternatives), attention-based multiple-instance learning with a regression head on HRD-sum (regression transfers better than binary labels, El Nahhas 2024) and a binary head per threshold, temperature-scaled calibration, Monte-Carlo dropout or ensemble uncertainty. The whole training run must fit in one day on one 8 GB GPU so that anyone can reproduce it.

C4. **External validation.** PTRC-HGSOC slides from TCIA: primary endpoint platinum refractory vs sensitive (open label), secondary endpoint scar metrics from Chowdhury et al. 2023 Table S1 where obtainable. Reported as pre-registered analysis with no tuning on the external set.

C5. **Deployable tool.** `hrdscope predict slide.svs` producing JSON, an attention heatmap and a PDF report; Docker image; CPU fallback; weights hosted on Zenodo and Hugging Face under MIT. Research-use-only disclaimer, no clinical claims.

C6. **Open science.** MIT code on GitHub with Zenodo DOI, embeddings released so that others can train aggregators without downloading 470 GB of slides, medRxiv preprint, submission to JOSS (software) and to an open-access journal that accepts preprinted manuscripts for the validation study.

## 4. Data

| Set | Use | Patients | Slides | Size | Labels |
|---|---|---|---|---|---|
| TCGA-OV diagnostic (FFPE) | train / CV, FFPE domain | 98 labeled of 106 | 107 | 209 GB | HRD-sum from C1 |
| TCGA-OV tissue (frozen) | train / CV, larger frozen domain | 569 labeled of 589 | 1,374 | 260 GB | HRD-sum from C1 |
| PTRC-HGSOC (FFPE) | external test only | 158 | 348 | 120 GB | platinum response (open), scar metrics (Table S1) |

Slides are streamed: download one file, tissue-mask, tile at 20x (0.5 micron per pixel, 224 px), embed, store embeddings as HDF5 (~10 MB per slide), delete the slide. Disk stays under 50 GB, embeddings for everything are ~15 GB and are released.

CPTAC-OV is the same patients as TCGA-OV and is therefore not used as external data. UBC-OCEAN is reserved for a later subtyping module (histotype), which gates the HRD model to HGSOC.

## 5. Methods

1. Tissue detection: Otsu on saturation at 2x, reject tiles with < 50 % tissue or pen marks.
2. Tiling: 224 px at 0.5 mpp (level chosen per slide from `openslide` metadata, resampled when needed), cap at 4,000 tiles per slide by random sampling for training, all tiles at inference.
3. Embedding: Midnight-12k CLS + mean patch token (1536-d), fp16, batch size fitted to free GPU memory; Hibou-B (768-d) as the light/CPU alternative.
4. Aggregation: attention MIL (Ilse 2018) with gated attention, dropout, multi-task output (HRD-sum regression, three binary logits). Patient-level bags (all slides of a patient) for training; slide-level prediction available.
5. Training: 5-fold stratified CV by patient, repeated with 3 seeds, AdamW, early stopping on validation loss, 8 GB budget.
6. Calibration and uncertainty: temperature scaling on out-of-fold predictions; ensemble of the 15 CV models at inference; conformal prediction sets at 90 % coverage.
7. Evaluation: AUROC, AUPRC with bootstrap CIs; calibration (ECE, reliability curves); decision-curve analysis; triage operating points: threshold giving sensitivity >= 0.95 for HRD (rule-out), report NPV and fraction of assays avoidable; threshold giving PPV >= 0.90 (rule-in). Frozen vs FFPE cross-domain transfer measured explicitly. Comparison across backbones and label thresholds.
8. Explainability: attention heatmaps; top-attended tiles per prediction reviewed against known HRD morphology (tumour-infiltrating lymphocytes, solid/pseudo-endometrioid/transitional patterns, nuclear pleomorphism).

## 6. Compute budget (this machine)

RTX 3070 with 8 GB (about 3.5 GB free while other services run), 20 CPU threads, 92 GB free disk, ~100 Mbit download.

| Step | Estimate |
|---|---|
| Download 470 GB of TCGA slides | measured 2.6 MB/s single stream from GDC: ~50 h; parallel range requests used to shorten it |
| Tiling + embedding, ~2,000-6,000 tiles/slide, ~1,500 slides | ~1-2 s per 100 tiles fp16 -> 1-2 GPU days |
| MIL training, 15 models | < 2 hours |
| PTRC-HGSOC embedding | ~ 1 day |

Total: about one week of unattended machine time.

## 7. Risks and mitigations

* Small FFPE cohort (98 patients): report CIs; use frozen slides for pretraining the aggregator; treat FFPE as the primary domain because that is what labs have.
* Label noise: HRD-sum is itself a proxy; publish continuous scores, evaluate at all three thresholds, and show concordance with Knijnenburg.
* Midnight-12k was pretrained on TCGA slides (self-supervised, no labels): document, and repeat the benchmark with Hibou-B / H-optimus-0 which were not.
* PTRC-HGSOC HRD labels may not be in the open supplement: platinum response remains a valid, clinically meaningful endpoint used by DeepHRD and IHGAMP.
* Regulatory: the tool is research-use-only; the report says so on every page.

## 8. Milestones

| # | Deliverable | Status |
|---|---|---|
| M1 | Open TCGA-OV HRD label table + scar module + tests | done |
| M2 | Streaming WSI pipeline (download, tile, embed), MIL, tuning, benchmark and predict commands with tests | done (code); GPU run pending |
| M3 | Embeddings for TCGA-OV FFPE (107 slides) + first CV results | |
| M4 | Embeddings for TCGA-OV frozen + PTRC-HGSOC; benchmark table | |
| M5 | `hrdscope predict` CLI, report, Docker, CPU path, weights on Zenodo/HF | |
| M6 | Preprint on medRxiv; JOSS submission; Zenodo DOI for v1.0 | |

## 9. Publication plan

* Software paper: Journal of Open Source Software (open, free, reviews the repository itself; requires tests, docs, CI, a statement of need).
* Validation study: medRxiv preprint first, then an open-access journal that explicitly allows preprints, e.g. Journal of Pathology Informatics, npj Digital Medicine, PLOS Digital Health, GigaScience, or Bioinformatics Advances.
* Data records: Zenodo deposits for (a) the label table, (b) tile embeddings per cohort, (c) trained weights, each with its own DOI, all CC BY 4.0 / MIT.
