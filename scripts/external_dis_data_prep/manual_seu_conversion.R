## conda incompatibility with anndata and other related packages
## ran these lines in R for subset h5ad conversion to seurat obj
## using rbase conda environment
## pbs interactive job with at least 150GB mem to avoid crashing

.libPaths("resources/r_package")
library(Seurat)
library(anndata)

for (i in 1:9) { rds_path <- sprintf("data/external_dis_data_prep/ren_filtered_split%02d.rds", i); if (!file.exists(rds_path)) { h5ad_path <- sprintf("data/external_dis_data_prep/ren_filtered_split%02d.h5ad", i); data <- read_h5ad(h5ad_path); seu <- CreateSeuratObject(counts = t(data$X), meta.data = data$obs); saveRDS(seu, rds_path); } };
for (i in 1:4) { rds_path <- sprintf("data/external_dis_data_prep/wellcome_filtered_split%02d.rds", i); if (!file.exists(rds_path)) { h5ad_path <- sprintf("data/external_dis_data_prep/wellcome_filtered_split%02d.h5ad", i); data <- read_h5ad(h5ad_path); seu <- CreateSeuratObject(counts = t(data$X), meta.data = data$obs); saveRDS(seu, rds_path); } };
for (i in 1:5) { rds_path <- sprintf("data/external_dis_data_prep/combat_filtered_split%02d.rds", i); if (!file.exists(rds_path)) { h5ad_path <- sprintf("data/external_dis_data_prep/combat_filtered_split%02d.h5ad", i); data <- read_h5ad(h5ad_path); seu <- CreateSeuratObject(counts = t(data$X), meta.data = data$obs); saveRDS(seu, rds_path); } };
for (i in 1:1) { rds_path <- sprintf("data/external_dis_data_prep/ch_filtered_split%02d.rds", i); if (!file.exists(rds_path)) { h5ad_path <- sprintf("data/external_dis_data_prep/ch_filtered_split%02d.h5ad", i); data <- read_h5ad(h5ad_path); seu <- CreateSeuratObject(counts = t(data$X), meta.data = data$obs); saveRDS(seu, rds_path); } };
for (i in 1:1) { rds_path <- sprintf("data/external_dis_data_prep/glaucoma_filtered_split%02d.rds", i); if (!file.exists(rds_path)) { h5ad_path <- sprintf("data/external_dis_data_prep/glaucoma_filtered_split%02d.h5ad", i); data <- read_h5ad(h5ad_path); seu <- CreateSeuratObject(counts = t(data$X), meta.data = data$obs); saveRDS(seu, rds_path); } };
for (i in 1:2) { rds_path <- sprintf("data/external_dis_data_prep/ra_filtered_split%02d.rds", i); if (!file.exists(rds_path)) { h5ad_path <- sprintf("data/external_dis_data_prep/ra_filtered_split%02d.h5ad", i); data <- read_h5ad(h5ad_path); seu <- CreateSeuratObject(counts = t(data$X), meta.data = data$obs); saveRDS(seu, rds_path); } };
for (i in 1:6) { rds_path <- sprintf("data/external_dis_data_prep/sle_filtered_split%02d.rds", i); if (!file.exists(rds_path)) { h5ad_path <- sprintf("data/external_dis_data_prep/sle_filtered_split%02d.h5ad", i); data <- read_h5ad(h5ad_path); seu <- CreateSeuratObject(counts = t(data$X), meta.data = data$obs); saveRDS(seu, rds_path); } };

