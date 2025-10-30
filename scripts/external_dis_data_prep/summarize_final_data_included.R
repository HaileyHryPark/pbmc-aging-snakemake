library(Seurat)
library(rio)
library(tidyr)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(purrr)

main_celltypes <- c("CD4 T", "CD8 T", "Mono", "B", "NK")

ren <- map(snakemake@input[["ren"]], readRDS)
wellcome <- map(snakemake@input[["wellcome"]], readRDS)
combat <- map(snakemake@input[["combat"]], readRDS)
ch <- map(snakemake@input[["ch"]], readRDS)
glaucoma <- map(snakemake@input[["glaucoma"]], readRDS)
ra <- map(snakemake@input[["ra"]], readRDS)
sle <- map(snakemake@input[["sle"]], readRDS)


objs <- list(ren, wellcome, combat, ch, glaucoma, ra, sle)
names(objs) <- c("ren","wellcome","combat","ch","glaucoma","ra", "sle")

summary <- lapply(as.list(names(objs)), function(obj){

	meta_all <- bind_rows(lapply(objs[[obj]], \(x) x@meta.data)) %>% mutate(age = as.integer(as.character(age)))
	print(obj)
	print(head(meta_all))
	print(summary(meta_all))

	return(data.frame(dataset = obj,
			  n_donors = n_distinct(meta_all$donor_id),
			  n_healthy_donors = n_distinct(meta_all %>% filter(disease == "normal") %>% pull(donor_id)),
			  n_disease_donors = n_distinct(meta_all %>% filter(disease != "normal") %>% pull(donor_id)),
			  n_samples = n_distinct(meta_all$sample_id),
			  n_cells_main = nrow(meta_all %>% filter(predicted.celltype.l1 %in% main_celltypes)),
			  n_cells_total = nrow(meta_all),
			  n_female = meta_all %>% filter(sex == "female") %>% pull(donor_id) %>% unique() %>% length(),
			  n_male = meta_all %>% filter(sex == "male") %>% pull(donor_id) %>% unique() %>% length(),
			  age_min = min(meta_all$age, na.rm = TRUE),
			  age_max = max(meta_all$age, na.rm = TRUE)))

})

export(bind_rows(summary), snakemake@output[["summary"]])
