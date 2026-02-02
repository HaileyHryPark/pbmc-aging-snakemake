.libPaths("resources/r_package")

library(dplyr)
library(tidyverse)
library(ggpubr)
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
female_df_exc <- female_df %>% filter(!(feature %in% male_df$feature))
male_df_exc <- male_df %>% filter(!(feature %in% female_df$feature))

print(table(female_df_exc$final_cluster))
print(table(male_df_exc$final_cluster))

df_list <- list(both_df, female_df, male_df, female_df_exc, male_df_exc)
names(df_list) <- c("both", "female", "male", "female_only", "male_only")

pdf(snakemake@output[["plot1"]], height = 16, width = 8)

all_res_df <- lapply(as.list(names(df_list)), function(n){

df <- df_list[[n]]

all_res <- lapply(as.list(unique(df$final_cluster)), function(clust){

	ct_res <- lapply(as.list(celltype_level), function(ct){
		print(paste(clust, ct))
		if(ct == "All celltype"){
			clustdf <- df %>% dplyr::filter(final_cluster == clust)
		}else{
			clustdf <- df %>% dplyr::filter(celltype == ct, final_cluster == clust)
		}
		if(nrow(clustdf) < 5) return(data.frame())
		results <- disease_enrichment(entities = unique(clustdf$gene), vocabulary = "HGNC", 
				database = "CURATED")
				#database = "ALL")
		results_df <- results@qresult %>% as.data.frame() %>% mutate(final_cluster = clust, celltype = ct)
		if(nrow(results_df) > 0){
			plot( results, type = "Enrichment", count =4,  cutoff= 0.05)
		}
		return(results_df)
	})
	return(bind_rows(ct_res))
})

return(bind_rows(all_res) %>% mutate(type = n))

})
dev.off()

export(bind_rows(all_res_df), snakemake@output[["res1"]])

## Exclude ribosomal genes

pdf(snakemake@output[["plot2"]], height = 16, width = 8)

all_res_df <- lapply(as.list(names(df_list)), function(n){

df <- df_list[[n]]
df <- df %>% dplyr::filter(!grepl("^(RPS|RPL|MRPS|MRPL|MT-)", gene))

all_res <- lapply(as.list(unique(df$final_cluster)), function(clust){

	ct_res <- lapply(as.list(celltype_level), function(ct){
		print(paste(clust, ct))
		if(ct == "All celltype"){
			clustdf <- df %>% dplyr::filter(final_cluster == clust)
		}else{
			clustdf <- df %>% dplyr::filter(celltype == ct, final_cluster == clust)
		}
		if(nrow(clustdf) < 5) return(data.frame())
		results <- disease_enrichment(entities = unique(clustdf$gene), vocabulary = "HGNC",
				#database = "ALL")
				database = "CURATED")
		results_df <- results@qresult %>% as.data.frame() %>% mutate(final_cluster = clust, celltype = ct)
		if(nrow(results_df) > 0){
			plot( results, type = "Enrichment", count =4,  cutoff= 0.05)
		}
		return(results_df)
	})
	return(bind_rows(ct_res))
})

return(bind_rows(all_res) %>% mutate(type = n))

})
dev.off()

export(bind_rows(all_res_df), snakemake@output[["res2"]])

