# Sex-specific nonlinear immune aging at single cell level
Snakemake pipeline for the analysis of sex-specific, non-linear immune aging trajectories from public human PBMC single-cell RNA-sequencing data. This repository contains all code used to generate the figures and results in our manuscript (see [Citation](#citation)).

## Contents
- [Overview](#overview)
- [Data source](#data-source)
- [Repository structure](#repository-structure)
- [Running the code](#running-the-code)
- [Citation](#citation)

---

## Overview

Aging is increasingly recognized as a non-linear process, but investigation of this non-linearity at single-cell resolution and separately by sex remains limited. This pipeline addresses that gap by leveraging four large-scale, healthy-donor human PBMC scRNA-seq datasets (~3.8 million cells, 1,828 donors, ages 19-97) to:

1. **Detect non-linear immune aging dynamics** at single-cell resolution using a differential expression sliding-window analysis (DE-SWAN), applied to both sexes combined and to each sex separately, to identify sex-aware aging features (SAFs).
2. **Cluster SAFs into temporal kinetic trajectories** to characterize shared and sex-specific aging programs.
3. **Link transcriptional trajectories to DNA methylation** by analyzing an independent whole-blood 450K methylation array cohort for age x sex interaction differentially methylated regions (DMRs).
4. **Build sex-stratified biological age clocks** (Elastic Net, XGBoost, and MLP) trained on sex-combined vs. sex-stratified SAFs, and interpret model predictions using SHAP values.

---

## Data source
All datasets used in this study are publicly available. Raw data are **not** included in this repository and must be downloaded separately into the paths listed below before running the pipeline.

### Internal datasets (scRNA-seq)
* [Onek1k](https://cellxgene.cziscience.com/collections/dde06e0f-ab3b-46be-96a2-a8082383c4a1): data/internal_data_prep/onek1k.h5ad
* [AIDAv2](https://cellxgene.cziscience.com/collections/ced320a1-29f3-47c1-a735-513c7084d508): data/internal_data_prep/aida.h5ad
* [Marina et al.](https://www.synapse.org/Synapse:syn49637038/files/): GEX_HTO_processed/all_pbmcs.tar.gz -> data/internal_data_prep/all_pbmcs/all_pbmcs_rna.h5ad
* [Perez et al.](https://cellxgene.cziscience.com/collections/436154da-bcf1-4130-9c8b-120ff9a888f2): data/internal_data_prep/perez.h5ad

### External datasets (scRNA-seq)
* [Ren et al.](https://cellxgene.cziscience.com/collections/0a839c4b-10d0-4d64-9272-684c49a2c8ba) (COVID-19): data/external_dis_data_prep/ren.h5ad
* [Wellcome](https://cellxgene.cziscience.com/collections/ddfad306-714d-4cc0-9985-d9072820c530) (COVID-19): data/external_dis_data_prep/wellcome.h5ad
* [COMBAT](https://cellxgene.cziscience.com/collections/8f126edf-5405-4731-8374-b5ce11f53e82) (COVID-19): data/external_dis_data_prep/combat.h5ad
* [Heimlich et al.](https://cellxgene.cziscience.com/collections/0aab20b3-c30c-4606-bd2e-d20dae739c45) (Clonal Hematopoiesis): data/external_dis_data_prep/ch.h5ad
* [GSE268936](https://cellxgene.cziscience.com/collections/de2cde16-c8d3-4a6d-80be-1be9e879aaca) (Glaucoma): data/external_dis_data_prep/glaucoma.h5ad
* [Binvignat et al.](https://cellxgene.cziscience.com/collections/e1a9ca56-f2ee-435d-980a-4f49ab7a952b) (Rheumatoid arthritis): data/external_dis_data_prep/ra.h5ad
* [Immunobiology of Aging (ImmAge)](https://cellxgene.cziscience.com/collections/60a2676d-9f37-46cc-9b02-c7370a53be9c): data/external_immage_data_prep/immage.h5ad
* [Sound Life (soundlife)](https://cellxgene.cziscience.com/collections/60a2676d-9f37-46cc-9b02-c7370a53be9c): data/external_soundlife_data_prep/soundlife_{group}.h5ad

### Other datasets (DNA methylation)
* [GSE35069](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE35069): data/dna_methylation/{filenames}.idat, data/dna_methylation/samples.csv

### Other datasets (NHANES 2017-2018 Serum Ferritin)
* [NHANES 2017-2018 demographics data](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/DEMO_J.htm): data/internal_clustering/DEMO_J.xpt
* [NHANES 2017-2018 ferritin data](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/FERTIN_J.htm): data/internal_clustering/FERTIN_J.xpt

---

## Repository structure

```
pbmc-aging-snakemake/
Snakefile              # Main Snakemake workflow: defines rule "all" and includes all sub-workflows
config.yaml            # Pipeline configuration 
rules/                 # Modular Snakemake rule files (.smk), one per analysis
scripts/                # R / Python scripts called by the rules
data/                  # Input data and intermediate processed data
tables/                 # Output summary tables 
plots/                  # Output figures 
resources/               # Reference / auxiliary resource files (e.g., gene sets, annotation references)
env/                    # Conda / R environment specification files
pbs_config/             # HPC (PBS) job submission configuration files
pbs_oe/                 # HPC job output/error logs
README.md
```

### Key rule files (`rules/`)
 
| Rule file | Purpose |
|---|---|
| `internal_data_prep.smk` | QC, filtering, and Azimuth cell type annotation of the four internal scRNA-seq datasets |
| `internal_pseudobulk.smk` | Generates cell type-level pseudobulk expression matrices |
| `internal_deswan.smk` | Runs DE-SWAN (sliding-window differential expression) to identify SAFs |
| `internal_clustering.smk` | Batch correction, LOESS trajectory fitting, and Mfuzz clustering of SAFs; functional/positional enrichment |
| `internal_correlation.smk` | Spearman correlation of SAFs with age and GSEA-based functional annotation |
| `internal_celltype_prop.smk` | Cell type proportion changes with age (propeller analysis) |
| `internal_clock.smk` | Trains and evaluates internal sex-combined and sex-stratified biological age clocks; SHAP analysis |
| `external_dis_data_prep.smk` / `external_sc_data_prep.smk` | QC and preparation of external disease/validation scRNA-seq cohorts |
| `external_immage_data_prep.smk` / `external_immage_analysis.smk` | Preparation and DE-SWAN analysis of the ImmAge external cohort |
| `external_soundlife_data_prep.smk` | Preparation of the Sound Life external cohort |
| `external_pseudobulk.smk` | Pseudobulk generation for external datasets |
| `external_clock.smk` | External validation of biological age clocks on healthy and disease cohorts |
| `dna_methylation.smk` | Processing and DMR analysis of whole-blood 450K methylation array data |
 
---

## Running the code

### Setup
 
```bash
git clone https://github.com/HaileyHryPark/pbmc-aging-snakemake.git
cd pbmc-aging-snakemake
conda env create -f env/initial_env.yaml
conda activate <env_name>
```
`initial_env.yaml` installs Snakemake and the dependencies. 
 
### Run
 
```bash
snakemake all --profile ./pbs_config --use-conda --use-singularity
```
This runs `rule all` in the `Snakefile` the way I ran, generating every table and figure used in the manuscript. 
 
> **HPC/environment notes:** This pipeline was run on HPC, using the PBS scheduler (job configs in `pbs_config/`, logs in `pbs_oe/`). R was run inside a Singularity container rather than via conda for some rules in the case of required R packages not available on conda/Bioconductor, which were installed manually and stored in `resources/r_package/`. If reproducing this pipeline on a different system, you may need to adapt the R environment setup (Singularity) and the scheduler submission commands accordingly.

---

## Citation

If you use this pipeline or code in your work, please cite:
Jacques Behmoaras, Harry Park, Nina Le Bert et al. Sex-specific trajectories of nonlinear immune aging at single cell level, 04 March 2026, PREPRINT (Version 1) available at Research Square [https://doi.org/10.21203/rs.3.rs-8944546/v1]

---
 
## Contact
 
For questions about the code or pipeline, please open an issue on this repository or contact [e0859928@u.nus.edu].
