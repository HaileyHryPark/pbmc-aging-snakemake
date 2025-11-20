library(rio)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(ggfortify)
library(GSVA)
library(msigdbr)
library(purrr)

celltypes <- c("CD4 T", "CD8 T", "NK", "B", "Mono")

## Function
getGSScores <- function(mat, ct){

	gs <- lapply(as.list(unique(annotgo$ID)), function(id){
		genes <- msig_hs %>% filter(gs_exact_source == id, !is.na(gene_symbol), gene_symbol != "NA") %>% pull(gene_symbol) %>% unique()
		genes <- paste(ct, genes, sep = ".")
		genes <- intersect(genes, rownames(mat))
		print(genes[1:5])
		return(genes)
	})
	names(gs) <- unique(annotgo$ID)

	scores <- gsva(gsvaParam(mat, gs))
	scores <- t(scores) %>% as.data.frame() %>% rownames_to_column("sample_id") %>% pivot_longer(!sample_id, names_to = "gsid", values_to = "score") %>% mutate(celltype = ct)

	print(head(scores))

	return(scores)
}

mat <- import(snakemake@input[["data"]])
meta <- mat %>% select(sample_id = rowname, age, dataset, ethnicity, sex)
mat <- mat %>% column_to_rownames("rowname") %>% select(-c(age, dataset, ethnicity, sex)) %>% t()
print(mat[1:5,1:5])

annotgo <- import(snakemake@input[["annotgo"]]) %>% filter(Category %in% c("RNA biosynthetic process", "Proteostasis", "OXPHOS/Energy metabolism"))

msig_hs <- msigdbr(species = "Homo sapiens")

res <- lapply(as.list(celltypes), function(ct){

	mat_ct <- mat[startsWith(rownames(mat), ct), ]
	print(mat_ct[1:5,1:5])

	return(getGSScores(mat, ct))

})

all_res <- bind_rows(res) %>% left_join(annotgo %>% select(ID, Description, Category) %>% distinct(), by = join_by(gsid == ID)) %>% left_join(meta, by = "sample_id") %>% arrange(Category)

print(head(all_res))

export(all_res, snakemake@output[["scores"]])

pdf(snakemake@output[["plot"]], width = 15, height = 3)
lapply(as.list(unique(all_res$gsid)), function(id){

	data <- all_res %>% filter(gsid == id)
	p <- ggscatter(data, x = "age", y = "score", alpha = 0.1, size = 0.3) +
			facet_grid(.~celltype)+
			geom_smooth(method = "loess", se = F) +
			labs(title = data$Description[1]) +
			theme_linedraw(base_size = 15)+
			theme(panel.grid = element_blank())
	plot(p)
})
dev.off()



