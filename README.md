# pbmc-aging-snakemake
## Contents
- [Overview](#overview)
- [Installation guide](#installation-guide)
- [Repo contents](#repo-contents)
- [Data source](#data-source)
- [Running the code](#running-the-code)
- [Citation](#citation)

## Overview

## Installation guide

## Repo contents

## Data source
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

### External datasets (DNA methylation)
* [GSE35069](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE35069): data/dna_methylation/{filenames}.idat, data/dna_methylation/samples.csv

### External datasets (NHANES 2017-2018 Serum Ferritin)
* [NHANES 2017-2018 demographics data](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/DEMO_J.htm): data/internal_clustering/DEMO_J.xpt
* [NHANES 2017-2018 ferritin data](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/FERTIN_J.htm): data/internal_clustering/FERTIN_J.xpt

## Running the code

## Citation
