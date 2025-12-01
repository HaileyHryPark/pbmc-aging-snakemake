library(dplyr)
library(tidyverse)
library(ggpubr)
library(ggplot2)
library(colorspace)
library(ComplexHeatmap)
library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)
library(ensembldb)
library(circlize)
library(rio)

cluster_level = c("Early\nincrease", "Early\ndecrease", "Continuous\ndecrease", "Irregular\nchange", "Late\nincrease", "Continuous\nincrease")
celltype_level = c("CD4 T", "CD8 T", "NK", "B", "Mono")

### Functions
runGSEA <- function(df){
  
  ## Prepare list for gsea
  df <- df %>% dplyr::filter(!is.na(ENTREZID))
  df <- df[!duplicated(df$ENTREZID),]
  print(head(df))

  gsea_list <- df$rho
  names(gsea_list) <- as.integer(df$ENTREZID)
  gsea_list <- sort(gsea_list, decreasing = TRUE)

  print(gsea_list)
  print(names(gsea_list))
  
  # Run GSEA
  #set.seed(123)
  #kegg_res <- gseKEGG(geneList = gsea_list,
  #      	  pvalueCutoff = 1,
  #           	  minGSSize = 3, 
  #                organism = "hsa",
  #                seed = TRUE)
  #kegg_res_df <- as.data.frame(kegg_res)
  #kegg_res_df$db <- "KEGG"
  
  set.seed(123)
  wp_res <- gseWP(geneList = gsea_list,
		  pvalueCutoff = 1,
                  organism = "Homo sapiens")
  wp_res_df <- as.data.frame(wp_res)
  wp_res_df$db <- "WP"
  
  set.seed(123)
  r_res <- gsePathway(geneList = gsea_list, 
		  pvalueCutoff = 1)
  r_res_df <- as.data.frame(r_res)
  r_res_df$db <- "Reactome"

#  if(df$celltype[1] != "NK"){
#  set.seed(123)
#  go_res <- gseGO(geneList=gsea_list, 
#             OrgDb = org.Hs.eg.db, 
#             ont ="BP", 
#	     eps = 1e-6,
#	     pvalueCutoff = 1,
#	     verbose = TRUE,
#             seed = TRUE)
#  go_res_df <- as.data.frame(go_res)
#  go_res_df$db <- "GOBP"
#  } else{
#  go_res <- NULL
#  go_res_df <- data.frame()
#  }
  
  #gsea_res <- list(kegg_res, wp_res, r_res, go_res)
  #names(gsea_res) <- c("KEGG", "WP", "Reactome", "GOBP")
  
  #gsea_res_df <- bind_rows(kegg_res_df, wp_res_df, r_res_df, go_res_df)
  
  gsea_res <- list(wp_res, r_res)
  names(gsea_res) <- c("WP", "Reactome")
  
  gsea_res_df <- bind_rows(wp_res_df, r_res_df)
  gsea_res_df$gene_name <- sapply(gsea_res_df$core_enrichment, function(x){
    genes <- unlist(strsplit(x, split = "/"))
    return(paste(annot_deg[annot_deg$ENTREZID %in% genes, "gene"], collapse = "/"))
  })

  return(list(gsea_res, gsea_res_df))
}

RunGSEAbyCelltype <- function(df, title){

fares <- lapply(as.list(unique(df$celltype)), function(ct){

	print(ct)
	ctdf <- df %>% dplyr::filter(celltype == ct)

	res <- runGSEA(ctdf)

	res_rds <- res[[1]]
	res_df <- res[[2]]

	if(nrow(res_df) == 0){
		res_df <- data.frame()
	}else{
		res_df <- res_df %>% mutate(celltype = ct, gsea_cluster = title)
	}
	return(list(res_rds, res_df))
})

df_all <- bind_rows(lapply(fares, `[[`, 2))
rds_all <- lapply(fares, `[[`, 1)
names(rds_all) <- unique(df$celltype)

return(list(rds_all, df_all))

}

## Main
both_cor <- import(snakemake@input[["both_cor"]])
female_cor <- import(snakemake@input[["female_cor"]]) 
male_cor <- import(snakemake@input[["male_cor"]])

print(dim(both_cor))
print(dim(female_cor))
print(dim(male_cor))

both_df <- import(snakemake@input[["both_df"]]) %>% dplyr::filter(!is.na(final_cluster), final_cluster != "") %>%
	left_join(both_cor, by = "feature") %>% mutate(type = "both")
female_df <- import(snakemake@input[["female_df"]]) %>% dplyr::filter(!is.na(final_cluster), final_cluster != "") %>%
	left_join(female_cor, by = "feature") %>% mutate(type = "female")
male_df <- import(snakemake@input[["male_df"]]) %>% dplyr::filter(!is.na(final_cluster), final_cluster != "") %>%
	left_join(male_cor, by = "feature") %>% mutate(type = "male")

print(summary(both_df))
print(summary(female_df))
print(summary(male_df))

annot_deg <- ensembldb::select(org.Hs.eg.db, keys = unique(c(both_df$gene, female_df$gene, male_df$gene)), keytype = "SYMBOL", columns = c("SYMBOL","ENTREZID")) %>% dplyr::select(gene = SYMBOL, ENTREZID)
print(head(annot_deg))

both_df <- both_df %>% dplyr::left_join(annot_deg, by = "gene")
female_df <- female_df %>% dplyr::left_join(annot_deg, by = "gene")
male_df <- male_df %>% dplyr::left_join(annot_deg, by = "gene")

df_list <- list(both_df, female_df, male_df)
names(df_list) <- c("both", "female", "male")
export(bind_rows(df_list), snakemake@output[["annot_cor"]])

## Run GSEA for each df
all_res <- lapply(as.list(names(df_list)), function(n){

df <- df_list[[n]]

res1 <- RunGSEAbyCelltype(df, "All cluster")

#clust_res <- lapply(as.list(unique(df$final_cluster)), function(cl){
#	
#	print(cl)
#	cldf <- df %>% dplyr::filter(final_cluster == cl)
#	res <- RunGSEAbyCelltype(cldf, cl) 
#
#	return(res)
#})
#
#df_all <- bind_rows(lapply(clust_res, `[[`, 2)) %>% mutate(type = n)
#rds_all <- lapply(clust_res, `[[`, 1)
#names(rds_all) <- unique(df$final_cluster)

df1 <- res1[[2]] %>% mutate(type = n)
rds1 <- res1[[1]]

return(list(rds1, df1))

})

rds <- lapply(all_res, `[[`, 1)
names(rds) <- c("both", "female", "male")

saveRDS(rds, snakemake@output[["res"]])
export(bind_rows(lapply(all_res, `[[`, 2)), snakemake@output[["res_df"]])

