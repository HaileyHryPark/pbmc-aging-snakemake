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
library(msigdbr)

cluster_level = c("Early\nincrease", "Early\ndecrease", "Continuous\ndecrease", "Irregular\nchange", "Late\nincrease", "Continuous\nincrease")
celltype_level = c("CD4 T", "CD8 T", "NK", "B", "Mono")

mitocarta <- readRDS(snakemake@input[["mitocarta"]])

### Functions
runMitoCarta <- function(features){
  
  genelist <- unique(features)
  print(genelist)
  
  set.seed(123)
  enr_mito <- enricher(genelist,
			TERM2GENE=mitocarta[[1]],
			TERM2NAME=mitocarta[[2]],
                         pvalueCutoff = 1, 
                         minGSSize = 3)
  if(is.null(enr_mito)){
        return(data.frame())
  }else{
    enr_mito_res <- enr_mito@result
  }
  print(head(enr_mito_res))
  
  if(nrow(enr_mito_res) == 0){
        return(data.frame())
  }
  
  print(head(enr_mito_res))

  return(enr_mito_res)
}

RunFAbyCluster <- function(df, title){

fares <- lapply(as.list(unique(df$final_cluster)), function(clust){

	print(clust)
	clustdf <- df %>% dplyr::filter(final_cluster == clust)

	fadf <- runMitoCarta(unique(clustdf$gene))

        if(nrow(fadf) == 0){
                return(data.frame())
        }else{
                fadf <- fadf %>% mutate(cluster = clust, term = paste(ID, Description, sep = "-"))
                return(fadf)
        }	
})

res <- bind_rows(fares)

return(res)

}

## Main
both_df <- import(snakemake@input[["both_df"]]) %>% dplyr::filter(!is.na(final_cluster))
female_df <- import(snakemake@input[["female_df"]]) %>% dplyr::filter(!is.na(final_cluster))
male_df <- import(snakemake@input[["male_df"]]) %>% dplyr::filter(!is.na(final_cluster))
female_df_exc <- female_df %>% dplyr::filter(!(feature %in% male_df$feature))
male_df_exc <- male_df %>% dplyr::filter(!(feature %in% female_df$feature))

print(table(female_df_exc$final_cluster))
print(table(male_df_exc$final_cluster))

df_list <- list(both_df, female_df, male_df, female_df_exc, male_df_exc)
names(df_list) <- c("both", "female", "male", "female_only", "male_only")

## Res 1
all_res <- lapply(as.list(names(df_list)), function(n){

df <- df_list[[n]]

res <- RunFAbyCluster(df, "All celltype")
res$fa_celltype <- "All celltype"

ct_res <- lapply(as.list(df %>% pull(celltype) %>% unique()), function(ct){

        print(ct)
        ctdf <- df %>% dplyr::filter(celltype == ct)
        res <- RunFAbyCluster(ctdf, ct)
        res$fa_celltype <- ct

        return(res)
})

return(rbind(res, bind_rows(ct_res)) %>% mutate(type = n))

})
export(bind_rows(all_res), snakemake@output[["res1"]])

### Exclude ribosomal genes
#annot_deg <- annot_deg %>% dplyr::filter(!grepl("^(RPS|RPL|MRPS|MRPL|MT-)", SYMBOL))
#
### Res 2
#all_res <- lapply(as.list(names(df_list)), function(n){
#
#df <- df_list[[n]]
#
#res <- RunFAbyCluster(df, "All celltype")
#res$fa_celltype <- "All celltype"
#
#ct_res <- lapply(as.list(df %>% pull(celltype) %>% unique()), function(ct){
#
#        print(ct)
#        ctdf <- df %>% dplyr::filter(celltype == ct)
#        res <- RunFAbyCluster(ctdf, ct)
#        res$fa_celltype <- ct
#
#        return(res)
#})
#
#return(rbind(res, bind_rows(ct_res)) %>% mutate(type = n))
#
#})
#
#export(bind_rows(all_res), snakemake@output[["res2"]])
