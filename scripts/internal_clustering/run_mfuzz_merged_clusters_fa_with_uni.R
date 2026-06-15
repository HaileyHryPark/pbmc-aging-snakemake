library(dplyr)
library(tidyverse)
library(ggpubr)
library(ggplot2)
library(colorspace)
library(ComplexHeatmap)
library(AnnotationHub)
library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)
library(ensembldb)
library(circlize)
library(rio)

cluster_level = c("Early\nincrease", "Early\ndecrease", "Continuous\ndecrease", "Early\nfluctuation", "Inverted\nUshape", "Continuous\nincrease", "Late\nincrease")
celltype_level = c("CD4 T", "CD8 T", "B", "NK", "Mono")

### Functions
runFA <- function(features, uni){
  
  genelist <- universe %>% dplyr::filter(gene %in% features, !is.na(ENTREZID)) %>% pull(ENTREZID) %>% unique() %>% as.character()
  print(genelist)
  
  print(head(uni))
  print(length(uni))
  
  set.seed(123)
  enr_kegg <- enrichKEGG(genelist,
                         organism = "hsa",
			 universe = uni,
                         pvalueCutoff = 1, 
                         minGSSize = 3)
  if(is.null(enr_kegg)){
    print("here1")
    enr_kegg_res <- NULL
  }else{
    enr_kegg_res <- enr_kegg@result
    enr_kegg_res$db <- "KEGG"
    enr_kegg_res <- enr_kegg_res[,-c(1,2)]
  }
  
  set.seed(123)
  enr_wp <- enrichWP(genelist,
                     organism = "Homo sapiens",
		     universe = uni,
                     pvalueCutoff = 1, 
                     minGSSize = 3)
  if(is.null(enr_wp)){
    print("here2")
    enr_wp_res <- NULL
  }else{
    enr_wp_res <- enr_wp@result
    enr_wp_res$db <- "WP"
  }
  
  set.seed(123)
  enr_r <- enrichPathway(genelist,
                         organism = "human",
			 universe = uni,
                         pvalueCutoff = 1, 
                         minGSSize = 3)
  if(is.null(enr_r)){
    print("here3")
    enr_r_res <- NULL
  }else{
    enr_r_res <- enr_r@result
    enr_r_res$db <- "Reactome"
  }

  res_list <- list(enr_kegg_res, enr_wp_res, enr_r_res)

  enr_all_res <- bind_rows(res_list[!unlist(lapply(res_list, is.null))])
  if(nrow(enr_all_res) == 0){
        return(data.frame())
  }
  
  enr_all_res$gene_name <- sapply(enr_all_res$geneID, function(x){
    genes <- unlist(strsplit(x, split = "/"))
    return(paste(universe[universe$ENTREZID %in% genes, "gene"], collapse = "/"))
  })
  return(enr_all_res)
}

runGOFA <- function(features, uni){
  
  genelist <- universe %>% dplyr::filter(gene %in% features, !is.na(ENTREZID)) %>% pull(ENTREZID) %>% unique() %>% as.character()
  print(genelist)
  
  print(head(uni))
  print(length(uni))
  
  set.seed(123)
  enr_go <- enrichGO(genelist,
                     OrgDb = org.Hs.eg.db,
		     universe = uni,
                     ont = "BP",
                     pvalueCutoff = 1, 
                     keyType = "ENTREZID",
                     minGSSize = 3)
  if(is.null(enr_go)){
    return(data.frame())
  }
  enr_go_res <- enr_go@result
  enr_go_res$db <- "GO"
  enr_go_res$gene_name <- sapply(enr_go_res$geneID, function(x){
    genes <- unlist(strsplit(x, split = "/"))
    return(paste(universe[universe$ENTREZID %in% genes, "gene"], collapse = "/"))
  })
  return(enr_go_res)
}

RunFAbyCluster <- function(df, title, uni){

fares <- lapply(as.list(unique(df$final_cluster)), function(clust){

	print(clust)
	clustdf <- df %>% dplyr::filter(final_cluster == clust)

	fadf <- rbind(runFA(unique(clustdf$gene), uni), runGOFA(unique(clustdf$gene), uni)) 

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
both_df <- import(snakemake@input[["both_df"]]) %>% dplyr::filter(!is.na(final_cluster), final_cluster != "")
female_df <- import(snakemake@input[["female_df"]]) %>% dplyr::filter(!is.na(final_cluster), final_cluster != "")
male_df <- import(snakemake@input[["male_df"]]) %>% dplyr::filter(!is.na(final_cluster), final_cluster != "")

universe <- import(snakemake@input[["uni"]])
print("Full universe dim")
print(dim(universe))

df_list <- list(both_df, female_df, male_df)
names(df_list) <- c("both", "female", "male")

## Res1
all_res <- lapply(as.list(names(df_list)), function(n){

df <- df_list[[n]]

ct_res <- lapply(as.list(df %>% pull(celltype) %>% unique()), function(ct){
	
	print(ct)
	ctdf <- df %>% dplyr::filter(celltype == ct)
	uni <- universe %>% dplyr::filter(celltype == ct, !is.na(ENTREZID)) %>% pull(ENTREZID) %>% unique() %>% as.character()
	print(head(uni))
	print(length(uni))
	res <- RunFAbyCluster(ctdf, ct, uni)
	res$fa_celltype <- ct

	return(res)
})

return(bind_rows(ct_res) %>% mutate(type = n))

})

export(bind_rows(all_res), snakemake@output[["res1"]])

