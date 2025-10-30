library(Seurat)
library(rio)
library(tidyr)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(purrr)

main_celltypes <- c("CD4 T", "CD8 T", "Mono", "B", "NK")

onek1k <- map(snakemake@input[["onek1k"]], readRDS)
aida <- map(snakemake@input[["aida"]], readRDS)
perez <- map(snakemake@input[["perez"]], readRDS)
marina <- map(snakemake@input[["marina"]], readRDS)


objs <- list(onek1k, aida, perez, marina)
names(objs) <- c("onek1k","aida","perez","marina")

summary <- lapply(as.list(names(objs)), function(obj){

	meta_all <- bind_rows(lapply(objs[[obj]], \(x) x@meta.data))

	return(data.frame(dataset = obj,
			  n_donors = n_distinct(meta_all$donor_id),
			  n_cells_main = nrow(meta_all %>% filter(predicted.celltype.l1 %in% main_celltypes)),
			  n_cells_total = nrow(meta_all),
			  n_female = meta_all %>% filter(sex == "female") %>% pull(donor_id) %>% unique() %>% length(),
			  n_male = meta_all %>% filter(sex == "male") %>% pull(donor_id) %>% unique() %>% length(),
			  n_asian = meta_all %>% filter(self_reported_ethnicity %in% c("Asian", "Indian", "Japanese", "Korean", "Singaporean Chinese", "Singaporean Indian", "Singaporean Malay", "Thai")) %>% pull(donor_id) %>% unique() %>% length(),
			  age_min = min(meta_all$age, na.rm = TRUE),
			  age_max = max(meta_all$age, na.rm = TRUE)))

})

export(bind_rows(summary), snakemake@output[["summary"]])
