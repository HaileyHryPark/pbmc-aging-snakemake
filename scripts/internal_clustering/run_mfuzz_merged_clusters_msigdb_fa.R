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

cluster_level = c("Early\nincrease", "Early\ndecrease", "Continuous\ndecrease", "Early\nfluctuation", "Late\nincrease", "Continuous\nincrease", "Inverted\nU-shape")
celltype_level = c("CD4 T", "CD8 T", "NK", "B", "Mono")

msig_hs <- msigdbr(species = "Homo sapiens")
print(unique(msig_hs$gs_cat))
print(unique(msig_hs$gs_subcat))
msig_imm <- msig_hs %>% dplyr::filter(gs_subcat == "IMMUNESIGDB") %>% dplyr::select(gs_name, entrez_gene)
msig_pos <- msig_hs %>% dplyr::filter(gs_cat == "C1") %>% dplyr::select(gs_name, entrez_gene)
msig_mir <- msig_hs %>% dplyr::filter(gs_subcat == "MIR:MIRDB") %>% dplyr::select(gs_name, entrez_gene)
msig_tft <- msig_hs %>% dplyr::filter(gs_subcat == "TFT:GTRD") %>% dplyr::select(gs_name, entrez_gene)

### Functions
runFA <- function(features){
  
  genelist <- annot_deg %>% dplyr::filter(SYMBOL %in% features) %>% pull(ENTREZID) %>% unique() %>% as.integer()
  print(genelist)
  
  set.seed(123)
  enr_imm <- enricher(genelist,
			TERM2GENE=msig_imm,
                         pvalueCutoff = 1, 
                         minGSSize = 3)
  if(is.null(enr_imm)){
    print("here1")
    enr_imm_res <- NULL
  }else{
    enr_imm_res <- enr_imm@result
    enr_imm_res$db <- "IMMUNESIGDB"
  }
  print(head(enr_imm_res))
  
  set.seed(123)
  enr_pos <- enricher(genelist,
			TERM2GENE=msig_pos,
                         pvalueCutoff = 1, 
                         minGSSize = 3)
  if(is.null(enr_pos)){
    print("here2")
    enr_pos_res <- NULL
  }else{
    enr_pos_res <- enr_pos@result
    enr_pos_res$db <- "C1"
  }
  print(head(enr_pos_res))
  
  set.seed(123)
  enr_mir <- enricher(genelist,
			TERM2GENE=msig_mir,
                         pvalueCutoff = 1, 
                         minGSSize = 3)
  if(is.null(enr_mir)){
    print("here3")
    enr_mir_res <- NULL
  }else{
    enr_mir_res <- enr_mir@result
    enr_mir_res$db <- "MIR"
  }
  print(head(enr_mir_res))

  set.seed(123)
  enr_tft <- enricher(genelist,
			TERM2GENE=msig_tft,
                         pvalueCutoff = 1, 
                         minGSSize = 3)
  if(is.null(enr_tft)){
    print("here3")
    enr_tft_res <- NULL
  }else{
    enr_tft_res <- enr_tft@result
    enr_tft_res$db <- "TFT"
  }
  print(head(enr_tft_res))

  res_list <- list(enr_imm_res, enr_pos_res, enr_mir_res, enr_tft_res)

  enr_all_res <- bind_rows(res_list[!unlist(lapply(res_list, is.null))])
  if(nrow(enr_all_res) == 0){
        return(data.frame())
  }
  
  enr_all_res$gene_name <- sapply(enr_all_res$geneID, function(x){
    genes <- unlist(strsplit(x, split = "/"))
    return(paste(annot_deg[annot_deg$ENTREZID %in% genes, "SYMBOL"], collapse = "/"))
  })

  print(head(enr_all_res))

  return(enr_all_res)
}

runGOFA <- function(features){
  
  genelist <- annot_deg %>% dplyr::filter(SYMBOL %in% features) %>% pull(ENTREZID) %>% unique() %>% as.integer()
  print(genelist)
  
  set.seed(123)
  enr_go <- enrichGO(genelist,
                     OrgDb = org.Hs.eg.db,
                     ont = "CC",
                     pvalueCutoff = 1, 
                     keyType = "ENTREZID",
                     minGSSize = 3)
  if(is.null(enr_go)){
    return(data.frame())
  }
  enr_go_res <- enr_go@result
  enr_go_res$db <- "GOCC"
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

## Main
df <- import(snakemake@input[["df"]]) %>% dplyr::filter(!is.na(final_cluster))

annot_deg <- ensembldb::select(org.Hs.eg.db, keys = unique(df$gene), keytype = "SYMBOL", columns = c("SYMBOL","ENTREZID"))
print(head(annot_deg))

## Res 1
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

## Res 1
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
