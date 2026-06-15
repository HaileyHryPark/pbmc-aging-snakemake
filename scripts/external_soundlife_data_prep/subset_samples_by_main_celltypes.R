#print(.libPaths())
#.libPaths(c("/scratch/users/nus/e0859928/Snakemake/onek1k-analysis-snakemake/resources/r_package", "/home/users/nus/e0859928/opt/miniforge3/envs/rbase/lib/R/library", .libPaths()))
#print(.libPaths())

library(Seurat)
library(rio)
library(tidyr)
library(dplyr)
library(tidyverse)
library(ggpubr)

main_celltypes <- c("CD4 T", "CD8 T", "Mono", "B", "NK")

## Could not get h5ad -> rds conversion. Manually changed with anndata.
seu <- readRDS(snakemake@input[["data"]])
print(seu)
print(unique(seu$predicted.celltype.l1))

meta <- seu@meta.data

cell_counts <- meta %>%
  count(sample_id, sex, age, predicted.celltype.l1) %>%
  pivot_wider(names_from = "predicted.celltype.l1", values_from = n, values_fill = 0)
print(cell_counts)

## Also removed sex unknown and age 89+
eligible_samples <- cell_counts %>%
  mutate(age = as.integer(as.character(age))) %>% 
  filter(if_all(all_of(main_celltypes), ~ .x >= 3), !is.na(age), sex %in% c("female","male")) %>%
  pull(sample_id)

seu_filtered <- subset(seu, subset = sample_id %in% eligible_samples)

saveRDS(seu_filtered, snakemake@output[["subset"]])
print("rds saved")

excluded_samples <- meta %>%
  filter(!(sample_id %in% eligible_samples)) %>%
  distinct(sample_id, sex, age, disease, self_reported_ethnicity)

export(excluded_samples, snakemake@output[["log"]])
