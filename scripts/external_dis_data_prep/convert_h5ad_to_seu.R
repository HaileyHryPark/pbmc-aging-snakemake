library(Seurat)
library(anndata)

rds_path <- snakemake@output[["seu"]]

if (!file.exists(rds_path)) { 
	h5ad_path <- snakemake@input[["data"]]
	data <- read_h5ad(h5ad_path) 
	seu <- CreateSeuratObject(counts = t(data$X), meta.data = data$obs) 
	saveRDS(seu, rds_path)
}
