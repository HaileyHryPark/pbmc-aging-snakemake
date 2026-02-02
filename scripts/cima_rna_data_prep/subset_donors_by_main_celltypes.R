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
  count(donor_id, predicted.celltype.l1) %>%
  pivot_wider(names_from = "predicted.celltype.l1", values_from = n, values_fill = 0)
print(cell_counts)

eligible_donors <- cell_counts %>%
  filter(if_all(all_of(main_celltypes), ~ .x >= 3)) %>%
  pull(donor_id)

seu_filtered <- subset(seu, subset = donor_id %in% eligible_donors)

saveRDS(seu_filtered, snakemake@output[["subset"]])
print("rds saved")

excluded_donors <- meta %>%
  filter(!(donor_id %in% eligible_donors)) %>%
  distinct(donor_id, sex, age, disease)

export(excluded_donors, snakemake@output[["log"]])
