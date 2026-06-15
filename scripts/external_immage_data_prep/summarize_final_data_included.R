library(Seurat)
library(rio)
library(tidyr)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(purrr)

main_celltypes <- c("CD4 T", "CD8 T", "Mono", "B", "NK")

data <- map(snakemake@input[["data"]], readRDS)

meta_all <- bind_rows(lapply(data, \(x) x@meta.data)) %>% mutate(age = as.integer(as.character(age)))
print(head(meta_all))
print(summary(meta_all))

summary <- data.frame(dataset = "immage",
			  n_donors = n_distinct(meta_all$donor_id),
			  n_healthy_donors = n_distinct(meta_all %>% filter(disease == "normal") %>% pull(donor_id)),
			  n_disease_donors = n_distinct(meta_all %>% filter(disease != "normal") %>% pull(donor_id)),
			  n_samples = n_distinct(meta_all$sample_id),
			  n_cells_main = nrow(meta_all %>% filter(predicted.celltype.l1 %in% main_celltypes)),
			  n_cells_total = nrow(meta_all),
			  n_female = meta_all %>% filter(sex == "female") %>% pull(donor_id) %>% unique() %>% length(),
			  n_male = meta_all %>% filter(sex == "male") %>% pull(donor_id) %>% unique() %>% length(),
			  age_min = min(meta_all$age, na.rm = TRUE),
			  age_max = max(meta_all$age, na.rm = TRUE))


export(summary, snakemake@output[["summary"]])
