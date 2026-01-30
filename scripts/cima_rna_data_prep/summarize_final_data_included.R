library(Seurat)
library(rio)
library(tidyr)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(purrr)

main_celltypes <- c("CD4 T", "CD8 T", "Mono", "B", "NK")

data <- readRDS(snakemake@input[["data"]])

meta <- sc@meta.data %>% mutate(age = as.character(age))
print(head(meta))
print(summary(meta))

df <- data.frame(dataset = "Supercentenarian",
			 n_donors = n_distinct(meta$donor_id),
			 n_control_donors = n_distinct(meta %>% filter(disease == "CT") %>% pull(donor_id)),
			 n_sc_donors = n_distinct(meta %>% filter(disease == "SC") %>% pull(donor_id)),
			 n_samples = n_distinct(meta$sample_id),
			 n_cells_main = nrow(meta %>% filter(predicted.celltype.l1 %in% main_celltypes)),
			 n_cells_total = nrow(meta),
			 n_female = meta %>% filter(sex == "female") %>% pull(donor_id) %>% unique() %>% length(),
			 n_male = meta %>% filter(sex == "male") %>% pull(donor_id) %>% unique() %>% length(),
			 age_min = "50s",
			 age_max = "110s")


export(df, snakemake@output[["summary"]])
