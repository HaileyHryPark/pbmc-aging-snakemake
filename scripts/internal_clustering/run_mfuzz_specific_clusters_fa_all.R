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

msig_hs <- msigdbr(species = "Homo sapiens")
msig_pos <- msig_hs %>% dplyr::filter(gs_cat == "C1") %>% dplyr::select(gs_name, entrez_gene)

### Functions
runFA <- function(features){
  
  genelist <- annot_deg %>% dplyr::filter(SYMBOL %in% features) %>% pull(ENTREZID) %>% unique() %>% as.integer()
  print(genelist)
  
  set.seed(123)
  enr_kegg <- enrichKEGG(genelist,
                         organism = "hsa",
                         pvalueCutoff = 1, 
                         keyType = "ncbi-geneid",
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
    return(paste(annot_deg[annot_deg$ENTREZID %in% genes, "SYMBOL"], collapse = "/"))
  })
  return(enr_all_res)
}

runGOFA <- function(features){
  
  genelist <- annot_deg %>% dplyr::filter(SYMBOL %in% features) %>% pull(ENTREZID) %>% unique() %>% as.integer()
  print(genelist)
  
  set.seed(123)
  enr_go <- enrichGO(genelist,
                     OrgDb = org.Hs.eg.db,
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
    return(paste(annot_deg[annot_deg$ENTREZID %in% genes, "SYMBOL"], collapse = "/"))
  })
  return(enr_go_res)
}

runC1FA <- function(features){
  
  genelist <- annot_deg %>% dplyr::filter(SYMBOL %in% features) %>% pull(ENTREZID) %>% unique() %>% as.integer()
  print(genelist)
  
  set.seed(123)
  enr_pos <- enricher(genelist,
                     TERM2GENE = msig_pos,
                     pvalueCutoff = 1, 
                     minGSSize = 3)
  if(is.null(enr_pos)){
    return(data.frame())
  }
  enr_pos_res <- enr_pos@result
  enr_pos_res$db <- "C1"
  enr_pos_res$gene_name <- sapply(enr_pos_res$geneID, function(x){
    genes <- unlist(strsplit(x, split = "/"))
    return(paste(annot_deg[annot_deg$ENTREZID %in% genes, "SYMBOL"], collapse = "/"))
  })
  return(enr_pos_res)
}

RunFAbySubcluster <- function(df, title){

fares <- lapply(as.list(unique(df$subcluster)), function(clust){

	print(clust)
	clustdf <- df %>% dplyr::filter(subcluster == clust)

	fadf <- bind_rows(runFA(unique(clustdf$gene)), runGOFA(unique(clustdf$gene)), runC1FA(unique(clustdf$gene))) 

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

annot_deg <- ensembldb::select(org.Hs.eg.db, keys = unique(c(both_df$gene, female_df$gene, male_df$gene)), keytype = "SYMBOL", columns = c("SYMBOL","ENTREZID"))
print(head(annot_deg))

subset1 <- intersect(female_df %>% dplyr::filter(final_cluster == "Continuous\nincrease") %>% pull(feature), 
          male_df %>% dplyr::filter(final_cluster == "Early\nincrease") %>% pull(feature))
subset2 <- intersect(female_df %>% dplyr::filter(final_cluster == "Early\nincrease") %>% pull(feature), 
          male_df %>% dplyr::filter(final_cluster == "Irregular\nchange") %>% pull(feature))

print(length(subset1))
print(length(subset2))

## DF with specific subclusters
df <- female_df %>% mutate(subcluster = ifelse(feature %in% subset1, "FCI_MEI",
						ifelse(feature %in% subset2, "FEI_MIC", NA))) %>%
			dplyr::filter(!is.na(subcluster))
print(dim(df))

## Res1
res <- RunFAbySubcluster(df, "All celltype")
res$fa_celltype <- "All celltype"

ct_res <- lapply(as.list(df %>% pull(celltype) %>% unique()), function(ct){
	
	print(ct)
	ctdf <- df %>% dplyr::filter(celltype == ct)
	res <- RunFAbySubcluster(ctdf, ct) 
	res$fa_celltype <- ct

	return(res)
})

export(rbind(res, bind_rows(ct_res)), snakemake@output[["res1"]])

## Exclude ribosomal genes
annot_deg <- annot_deg %>% dplyr::filter(!grepl("^(RPS|RPL|MRPS|MRPL|MT-)", SYMBOL))

## Res2
res <- RunFAbySubcluster(df, "All celltype")
res$fa_celltype <- "All celltype"

ct_res <- lapply(as.list(df %>% pull(celltype) %>% unique()), function(ct){
	
	print(ct)
	ctdf <- df %>% dplyr::filter(celltype == ct)
	res <- RunFAbySubcluster(ctdf, ct) 
	res$fa_celltype <- ct

	return(res)
})

export(rbind(res, bind_rows(ct_res)), snakemake@output[["res2"]])
