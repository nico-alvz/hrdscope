# Open, reproducible triage of homologous recombination deficiency from routine H&E slides in ovarian cancer

Draft v0.1, 2026-09-18. Content-first draft; numbers marked [TBD] are filled from `runs/` when the benchmark completes.

## Abstract

Background. Roughly half of high-grade serous ovarian carcinomas (HGSOC) carry homologous recombination deficiency (HRD) and benefit from PARP inhibitors, but the genomic assays that establish HRD status are costly and unavailable in most health systems. Deep learning can recover part of the HRD signal from routine haematoxylin and eosin (H&E) slides, yet no published model is usable outside its authors' institution: weights are withheld or non-commercial, labels are private and thresholds differ.

Methods. We built hrdscope, an MIT-licensed pipeline that (i) derives HRD-sum labels for 555 TCGA-OV patients from open-access allele-specific copy-number profiles with a reimplementation of the three genomic scar scores (validated against published scores, Pearson r = 0.96 on 168 shared patients), (ii) embeds slides with an MIT-licensed pathology foundation model (Midnight-12k) and trains a gated attention multiple-instance model with a continuous HRD-sum head and binary heads at the 33, 42 and 63 thresholds, (iii) fixes patient-level splits and a benchmark protocol so that new backbones can be compared, and (iv) reports calibrated rule-out and rule-in operating points for use as a triage step ahead of confirmatory sequencing. External validation used PTRC-HGSOC (158 patients, CC BY 4.0) with platinum response as the endpoint.

Results. [TBD] Cross-validated AUROC for HRD-sum >= 42 on FFPE diagnostic slides was [TBD] (95 % CI [TBD]); on frozen tissue slides [TBD]. At a rule-out operating point with sensitivity 0.95, [TBD] % of patients could be deferred from sequencing with negative predictive value [TBD]. On PTRC-HGSOC the model separated platinum-refractory from sensitive tumours with AUROC [TBD].

Conclusions. hrdscope is the first HRD-from-H&E tool that a laboratory can install and run today under a permissive license with downloadable weights, and the first open patient-level HRD label table for the whole TCGA-OV cohort. Code, labels, embeddings and weights are on GitHub and Zenodo.

## 1. Introduction

* Burden of ovarian cancer; HGSOC; PARP inhibitor indications (SOLO-1, PRIMA, PAOLA-1) and the role of HRD testing.
* Cost and access: myChoice CDx, in-house WGS/WES scar assays, SNP arrays; turnaround; availability in low- and middle-income countries.
* H&E slides are universal; digital pathology scanners are spreading; foundation models make slide-level learning feasible on modest hardware.
* Prior work (DeepHRD, IHGAMP, HRDPath, cross-cancer transfer learning, Loeffler et al., El Nahhas et al.) shows moderate accuracy (AUROC 0.7-0.85) but none ships a usable, permissively licensed tool; labels are inconsistent (thresholds 33, 42, 63) and mostly private.
* Contributions (C1-C6 as in `docs/PLAN.md`).

## 2. Methods

### 2.1 Open HRD labels for TCGA-OV

* GDC open-access ASCAT3 allele-specific segments (555 tumour profiles, GRCh38).
* Scar scores: HRD-LOH (Abkevich 2012), NtAI (Birkbak 2012), LST (Popova 2012); implementation details (autosomes only, segment merging, 15 Mb / 1 Mb / 10 Mb and 3 Mb smoothing thresholds, centromeres from UCSC hg38 cytoBand).
* Validation against Knijnenburg et al. 2018: r = 0.96 (HRD-sum), 0.97 (LOH), 0.94 (LST), 0.82 (NtAI); systematic offset of +5 points in HRD-sum; binary agreement 0.89 at 42 and 0.91 at 63. Distribution: median 50; 65.0 % >= 42; 32.1 % >= 63; 81.4 % >= 33.
* Why three thresholds are reported.

### 2.2 Cohorts

* TCGA-OV diagnostic FFPE (95 labeled patients, 96 slides) and frozen tissue (554 labeled patients, ~1,300 slides).
* PTRC-HGSOC (158 patients, 348 FFPE slides), endpoint platinum refractory vs sensitive; secondary scar metrics from Chowdhury 2023 Table S1.
* CPTAC-OV excluded (same patients as TCGA-OV).

### 2.3 Slide processing

* Tissue mask (HSV saturation/value thresholds, pen-mark rejection at 16 micron/px), 224 px tiles at 0.5 micron/px, minimum 50 % tissue.
* Backbones: Midnight-12k (MIT; CLS + mean patch token, 3072-d), Hibou-B and H-optimus-0 (Apache-2.0) as comparators.
* Streaming download with parallel range requests, md5 check, slide deleted after embedding; embeddings released.

### 2.4 Model and training

* Gated attention MIL (Ilse et al. 2018), multi-task loss (MSE on HRD-sum / scale + BCE at three thresholds), AdamW, early stopping on inner validation fold.
* Fixed stratified 5-fold patient-level splits (seed 20260918); 3 seeds; nested Optuna search (30 trials per outer fold) over hidden size, attention size, dropout, learning rate, weight decay, regression scale and tiles per bag, using only inner folds.
* Ensemble of the 15 cross-validation models at inference; calibration by temperature scaling on out-of-fold predictions; uncertainty from ensemble spread.

### 2.5 Evaluation

* AUROC and AUPRC with 1,000 bootstrap CIs; Pearson/Spearman on HRD-sum; calibration curves and ECE; decision curve analysis.
* Triage operating points: rule-out (sensitivity >= 0.95; report NPV and fraction deferred) and rule-in (PPV >= 0.90; fraction prioritised).
* Cross-domain: train frozen, test FFPE and vice versa.
* External: PTRC-HGSOC platinum response; pre-specified, no tuning.

## 3. Results

[TBD from runs/dx_midnight/report.json, runs/ts_midnight/report.json, runs/ptrc/]

Tables: label concordance; benchmark by backbone and threshold; triage operating points; external validation. Figures: scar-score scatter vs Knijnenburg; ROC and calibration; decision curves; attention heatmaps with top tiles.

## 4. Discussion

* What the numbers mean for triage in a sequencing-limited setting; worked example with prevalence 50 %.
* Limitations: label is a proxy (scar score, not functional HRD); small FFPE cohort; TCGA pretraining of Midnight-12k; frozen-to-FFPE shift; no prospective data; research use only.
* Comparison with DeepHRD, IHGAMP, HRDPath under the same thresholds where possible.
* Open science: everything reproducible on one 8 GB GPU in about a week; how to add a backbone or a cohort.

## 5. Data and code availability

Code is on GitHub under MIT (https://github.com/nico-alvz/hrdscope) and archived on Zenodo (concept DOI 10.5281/zenodo.22827123); labels, embeddings and weights receive their own Zenodo DOIs at each release; Hugging Face for weights; all inputs are open access (GDC, TCIA CC BY 4.0, Cell Reports supplement).

## References

[TBD: Abkevich 2012; Birkbak 2012; Popova 2012; Sztupinszki 2018 scarHRD; Knijnenburg 2018; Bergstrom 2024; Zafar 2025/2026; HRDPath 2025; Wang 2026; Loeffler 2024; El Nahhas 2024; Ilse 2018; Karasikov/kaiko Midnight 2025; Chowdhury 2023; TCIA Clark 2013; TCGA 2011.]
