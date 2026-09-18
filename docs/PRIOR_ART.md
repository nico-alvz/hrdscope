# Prior art: HRD prediction from H&E in ovarian cancer (checked 2026-09-18)

| Work | Year / venue | Ovarian cohorts (n) | HRD label | Backbone | Ovarian result | Code | Weights | License | Deployable tool |
|---|---|---|---|---|---|---|---|---|---|
| DeepHRD (Bergstrom et al.) | 2024, J Clin Oncol | TCGA-OV frozen + FFPE, external Dijon | HRD score (Myriad-like) | Custom CNN, 2-stage MIL (5x then 20x) | AUC ~0.8 (breast/ovarian), platinum response | Yes (github AlexandrovLab/DeepHRD) | Not downloadable ("not currently supported") | UCSD Academic Software License, non-commercial | No (research scripts, no Docker) |
| IHGAMP (Zafar et al.) | 2025 medRxiv, 2026 APMIS | TCGA pan-cancer (OV inside 8,109 labeled); PTRC-HGSOC (158) for platinum response | scarHRD >= 33 (top-20% of TCGA) | OpenCLIP, UNI, OpenSlideFM + MIL | pan-cancer AUROC 0.77-0.83; PTRC platinum AUROC 0.67 | Yes, MIT (19 notebooks) | joblib heads; UNI variant is non-commercial | MIT code; weights inherit backbone license | No CLI/Docker; ~144 GPU hours to reproduce |
| HRDPath | 2025 bioRxiv | UWOV (152, private), TCGA-OV (81), PTRC-HGSOC (93) | HRDetect (UWOV); LOH+LST+TAI > 42 (TCGA/PTRC) | DSMIL + ViT ensemble | AUC 0.85, sens 0.69, spec 0.94 | Not stated | No | Not stated | No |
| Cross-cancer transfer learning (Wang et al.) | 2026, Sci Rep | TCGA-OV (561) | scarHRD >= 42 | UNI + CLAM-SB | OV gains small ("negative transfer") | Yes | Not released | Not stated | No |
| MSK + Zhongda OV model | 2025-2026, Int J Gynecol Cancer | 205 slides, 205 patients | assay HRD status | tumour segmentation + classifier | reported clinically meaningful AUC | Not stated | No | Not stated | No |
| Loeffler et al. | 2024, BMC Biology | TCGA 9-10 tumour types incl. OV | HRD score thresholds | attention MIL (marugoto/STAMP) | modest in OV | Yes (STAMP, MIT) | No HRD weights | MIT | STAMP is a generic pipeline, no HRD model shipped |
| El Nahhas et al. | 2024, Nature Cancer | TCGA pan-cancer incl. HRD | continuous HRD score | regression MIL | regression beats classification | Yes | No | -- | No |
| HRD-informed digital histology (Springer, Interdiscip Sci) | 2025 | HGSOC | HRD | -- | platinum response, prognosis | Not stated | No | -- | No |

## Labels used in the field

* Knijnenburg et al. 2018 (Cell Rep) HRD scores: open, but only 173 TCGA-OV patients have a score.
* scarHRD (R package) run by each group on their own ASCAT/Sequenza calls: not released as a table, thresholds differ (33, 42, 63).
* Commercial assays (Myriad myChoice) on private cohorts.

## Open pathology foundation models (weights usable in an MIT tool)

| Model | Params | Dim | License | Gated | Note |
|---|---|---|---|---|---|
| Midnight-12k (kaiko.ai) | 1.1 B | 1536 | MIT | No | trained only on TCGA slides |
| OpenMidnight (Sophont) | 1.1 B | 1536 | open replication of Midnight | No | |
| H-optimus-0 (Bioptimus) | 1.1 B | 1536 | Apache-2.0 | Yes (contact info) | 0.5 mpp |
| Hibou-B (HistAI) | 86 M | 768 | Apache-2.0 | Yes (contact info) | light enough for CPU |
| UNI / UNI2, CONCH, Virchow2, Prov-GigaPath | -- | -- | CC BY-NC or custom | Yes | not usable in an MIT tool |

## Public ovarian cohorts with slides

| Cohort | Patients / slides | Access | Labels available | Note |
|---|---|---|---|---|
| TCGA-OV diagnostic (FFPE) | 106 / 107 | GDC open, 209 GB | 98 with open ASCAT3 (this repo) | main FFPE training set |
| TCGA-OV tissue (frozen) | 589 / 1,374 | GDC open, 260 GB | 569 with open ASCAT3 | larger, frozen domain |
| PTRC-HGSOC | 158 / 348 FFPE | TCIA, CC BY 4.0, 120 GB | platinum refractory vs sensitive (open); scar metrics in Chowdhury 2023 Table S1 | independent external validation |
| CPTAC-OV | 102 / 222 | TCIA, CC BY 3.0 | -- | same patients as TCGA-OV, not independent |
| UBC-OCEAN | 527 WSIs | Kaggle | histotype | for subtyping, not HRD |

## Gap this project fills

1. No HRD-from-H&E tool exists that a lab can install and run today under a permissive license with downloadable weights.
2. No open, patient-level HRD label table covers the whole TCGA-OV cohort; groups recompute private labels with different thresholds, so results are not comparable.
3. No published work reports calibrated triage operating points (rule-out / rule-in) with decision-curve analysis, which is the actual use case when sequencing is scarce.
4. No fully open benchmark (open backbone + open labels + fixed splits + external cohort) exists for ovarian HRD.
