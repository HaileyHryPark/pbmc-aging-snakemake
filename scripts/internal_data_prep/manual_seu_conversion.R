## conda incompatibility with anndata and other related packages
## ran these lines in R for subset h5ad conversion to seurat obj
## using rbase conda environment
## pbs interactive job with at least 150GB mem to avoid crashing

.libPaths("resources/r_package")
library(Seurat)
library(anndata)

for (i in 1:21) { rds_path <- sprintf("data/internal_data_prep/aida_filtered_split%02d.rds", i); if (!file.exists(rds_path)) { h5ad_path <- sprintf("data/internal_data_prep/aida_filtered_split%02d.h5ad", i); data <- read_h5ad(h5ad_path); seu <- CreateSeuratObject(counts = t(data$X), meta.data = data$obs); saveRDS(seu, rds_path); } };
for (i in 1:6) { rds_path <- sprintf("data/internal_data_prep/marina_filtered_split%02d.rds", i); if (!file.exists(rds_path)) { h5ad_path <- sprintf("data/internal_data_prep/marina_filtered_split%02d.h5ad", i); data <- read_h5ad(h5ad_path); seu <- CreateSeuratObject(counts = t(data$X), meta.data = data$obs); saveRDS(seu, rds_path); } };
for (i in 1:10) { rds_path <- sprintf("data/internal_data_prep/onek1k_filtered_split%02d.rds", i); if (!file.exists(rds_path)) { h5ad_path <- sprintf("data/internal_data_prep/onek1k_filtered_split%02d.h5ad", i); data <- read_h5ad(h5ad_path); seu <- CreateSeuratObject(counts = t(data$X), meta.data = data$obs); saveRDS(seu, rds_path); } };
for (i in 1:4) { rds_path <- sprintf("data/internal_data_prep/perez_filtered_split%02d.rds", i); if (!file.exists(rds_path)) { h5ad_path <- sprintf("data/internal_data_prep/perez_filtered_split%02d.h5ad", i); data <- read_h5ad(h5ad_path); seu <- CreateSeuratObject(counts = t(data$X), meta.data = data$obs); saveRDS(seu, rds_path); } };

