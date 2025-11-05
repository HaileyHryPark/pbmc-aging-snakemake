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

cluster_level = c("Early\nincrease", "Early\ndecrease", "Continuous\ndecrease", "Irregular\nchange", "Late\nincrease", "Continuous\nincrease", "Inverted\nU-shape")
celltype_level = c("CD4 T", "CD8 T", "NK", "B", "Mono")

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

RunFAbyCluster <- function(df, title){

fares <- lapply(as.list(unique(df$final_cluster)), function(clust){

	print(clust)
	clustdf <- df %>% dplyr::filter(final_cluster == clust)

	fadf <- rbind(runFA(unique(clustdf$gene)), runGOFA(unique(clustdf$gene))) 

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

df <- import(snakemake@input[["df"]]) %>% dplyr::filter(!is.na(final_cluster))

print(table(df$final_cluster))
print(df %>% dplyr::filter(is.na(final_cluster)))

annot_deg <- ensembldb::select(org.Hs.eg.db, keys = unique(df$gene), keytype = "SYMBOL", columns = c("SYMBOL","ENTREZID"))
print(head(annot_deg))


## Res1
res <- RunFAbyCluster(df, "All celltype")
res$fa_celltype <- "All celltype"

ct_res <- lapply(as.list(df %>% pull(celltype) %>% unique()), function(ct){
	
	print(ct)
	ctdf <- df %>% dplyr::filter(celltype == ct)
	res <- RunFAbyCluster(ctdf, ct) 
	res$fa_celltype <- ct

	return(res)
})

export(rbind(res, bind_rows(ct_res)), snakemake@output[["res1"]])

## Exclude ribosomal genes
annot_deg <- annot_deg %>% dplyr::filter(!grepl("^(RPS|RPL|MRPS|MRPL|MT-)", SYMBOL))

## Res2
res <- RunFAbyCluster(df, "All celltype")
res$fa_celltype <- "All celltype"

ct_res <- lapply(as.list(df %>% pull(celltype) %>% unique()), function(ct){
	
	print(ct)
	ctdf <- df %>% dplyr::filter(celltype == ct)
	res <- RunFAbyCluster(ctdf, ct) 
	res$fa_celltype <- ct

	return(res)
})

export(rbind(res, bind_rows(ct_res)), snakemake@output[["res2"]])

