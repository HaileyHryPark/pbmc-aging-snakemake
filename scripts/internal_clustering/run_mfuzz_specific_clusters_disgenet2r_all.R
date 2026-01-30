.libPaths(c("/scratch/users/nus/e0859928/Snakemake/onek1k-analysis-snakemake/resources/r_package", "/home/users/nus/e0859928/opt/miniforge3/envs/rbase/lib/R/library", .libPaths()))

library(dplyr)
library(tidyverse)
library(ggplot2)
library(rio)
library(disgenet2r)


api_key <- "b9ff0e60-2dc1-4c68-8e77-28bc80901780"
Sys.setenv(DISGENET_API_KEY= api_key)

cluster_level = c("Early\nincrease", "Early\ndecrease", "Continuous\ndecrease", "Early\nfluctuation", "Late\nincrease", "Continuous\nincrease")
celltype_level = c("All celltype", "CD4 T", "CD8 T", "NK", "B", "Mono")

## Main
both_df <- import(snakemake@input[["both_df"]]) %>% dplyr::filter(!is.na(final_cluster))
female_df <- import(snakemake@input[["female_df"]]) %>% dplyr::filter(!is.na(final_cluster))
male_df <- import(snakemake@input[["male_df"]]) %>% dplyr::filter(!is.na(final_cluster))

subset1 <- intersect(female_df %>% dplyr::filter(final_cluster == "Continuous\nincrease") %>% pull(feature),
          male_df %>% dplyr::filter(final_cluster == "Early\nincrease") %>% pull(feature))
subset2 <- intersect(female_df %>% dplyr::filter(final_cluster == "Early\nincrease") %>% pull(feature),
          male_df %>% dplyr::filter(final_cluster == "Early\nfluctuation") %>% pull(feature))

## DF with specific subclusters
df <- female_df %>% mutate(subcluster = ifelse(feature %in% subset1, "FCI_MEI",
                                                ifelse(feature %in% subset2, "FEI_MEF", NA))) %>%
                        dplyr::filter(!is.na(subcluster))
print(dim(df))

pdf(snakemake@output[["plot"]], height = 16, width = 8)

all_res <- lapply(as.list(unique(df$subcluster)), function(sub){
ct_res <- lapply(as.list(celltype_level), function(ct){
	if(ct == "All celltype"){
		clustdf <- df %>% dplyr::filter(subcluster == sub)
	}else{
		clustdf <- df %>% dplyr::filter(celltype == ct, subcluster == sub)
	}
	if(nrow(clustdf) < 5) return(data.frame())
	results <- disease_enrichment(entities = unique(clustdf$gene), vocabulary = "HGNC", database = "ALL", common_entities = 3)

	results_df <- results@qresult %>% as.data.frame() %>% mutate(celltype = ct, subcluster = sub)
	if(nrow(results_df) > 0){
		plot( results, type = "Enrichment", count =4,  cutoff= 0.05)
	}
	return(results_df)
})
return(bind_rows(ct_res))
})
export(bind_rows(all_res), snakemake@output[["res"]])
dev.off()

