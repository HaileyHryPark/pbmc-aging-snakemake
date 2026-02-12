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

cluster_level = c("Early\nincrease", "Early\ndecrease", "Continuous\ndecrease", "Early\nfluctuation", "Inverted\nUshape", "Continuous\nincrease", "Late\nincrease")
celltype_level = c("CD4 T", "CD8 T", "B", "NK", "Mono")

msig_hs <- msigdbr(species = "Homo sapiens")
print(unique(msig_hs$gs_collection))
print(unique(msig_hs$gs_subcollection))
msig_pos <- msig_hs %>% dplyr::filter(gs_collection == "C1") %>% dplyr::select(gs_name, ncbi_gene)

### Functions
runFA <- function(features){
  
  genelist <- annot_deg %>% dplyr::filter(SYMBOL %in% features) %>% pull(ENTREZID) %>% unique() %>% as.integer()
  print(genelist)
  
  set.seed(123)
  enr_pos <- enricher(genelist,
			TERM2GENE=msig_pos,
                         pvalueCutoff = 1, 
			 qvalueCutoff = 0.5,
                         minGSSize = 3)
  if(is.null(enr_pos)){
    print("here2")
    enr_pos_res <- NULL
  }else{
    enr_pos_res <- enr_pos@result
    enr_pos_res$db <- "C1"
  }
  print(head(enr_pos_res))
  

  enr_all_res <- enr_pos_res

  if(is.null(enr_all_res)){
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

	fadf <- runFA(unique(clustdf$gene))

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

annot_deg <- ensembldb::select(org.Hs.eg.db, keys = unique(c(both_df$gene, female_df$gene, male_df$gene)), keytype = "SYMBOL", columns = c("SYMBOL","ENTREZID"))
print(head(annot_deg))

df_list <- list(both_df, female_df, male_df)
names(df_list) <- c("both", "female", "male")

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

