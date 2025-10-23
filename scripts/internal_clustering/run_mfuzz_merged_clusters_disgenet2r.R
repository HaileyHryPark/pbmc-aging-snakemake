.libPaths(c("/scratch/users/nus/e0859928/Snakemake/onek1k-analysis-snakemake/resources/r_package", "/home/users/nus/e0859928/opt/miniforge3/envs/rbase/lib/R/library", .libPaths()))

library(dplyr)
library(tidyverse)
library(ggpubr)
library(ggplot2)
library(rio)
library(disgenet2r)


api_key <- "b9ff0e60-2dc1-4c68-8e77-28bc80901780"
Sys.setenv(DISGENET_API_KEY= api_key)

cluster_level = c("Early\nincrease", "Early\ndecrease", "Continuous\ndecrease", "Irregular\nchange", "Late\nincrease", "Continuous\nincrease", "Inverted\nU-shape")
celltype_level = c("CD4 T", "CD8 T", "NK", "B", "Mono")

df <- import(snakemake@input[["var_cluster_df"]])

pdf(snakemake@output[["plot1"]], height = 16, width = 8)
all_res <- lapply(as.list(unique(df$final_cluster)), function(clust){

	ct_res <- lapply(as.list(celltype_level), function(ct){
		print(paste(clust, ct))
		clustdf <- df %>% dplyr::filter(celltype == ct, final_cluster == clust)
		if(nrow(clustdf) < 5) return(data.frame())
		results <- disease_enrichment(entities = unique(clustdf$gene), vocabulary = "HGNC")
				#database = "ALL")
		results_df <- results@qresult %>% as.data.frame() %>% mutate(final_cluster = clust, celltype = ct)
		if(nrow(results_df) > 0){
			plot( results, type = "Enrichment", count =4,  cutoff= 0.05)
		}
		return(results_df)
	})
	return(bind_rows(ct_res))
})
dev.off()

export(bind_rows(all_res), snakemake@output[["res1"]])

## Exclude ribosomal genes
df <- df %>% dplyr::filter(!grepl("^(RPS|RPL|MRPS|MRPL|MT-)", gene))

pdf(snakemake@output[["plot2"]], height = 16, width = 8)
all_res <- lapply(as.list(unique(df$final_cluster)), function(clust){

	ct_res <- lapply(as.list(celltype_level), function(ct){
		print(paste(clust, ct))
		clustdf <- df %>% dplyr::filter(celltype == ct, final_cluster == clust)
		if(nrow(clustdf) < 5) return(data.frame())
		results <- disease_enrichment(entities = unique(clustdf$gene), vocabulary = "HGNC",
				#database = "ALL")
				database = "ALL")
		results_df <- results@qresult %>% as.data.frame() %>% mutate(final_cluster = clust, celltype = ct)
		if(nrow(results_df) > 0){
			plot( results, type = "Enrichment", count =4,  cutoff= 0.05)
		}
		return(results_df)
	})
	return(bind_rows(ct_res))
})
dev.off()

export(bind_rows(all_res), snakemake@output[["res2"]])

