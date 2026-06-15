library(Seurat)
options(Seurat.object.assay.version = "v5")
library(rio)
library(tidyr)
library(dplyr)
library(tidyverse)
library(ggpubr)

options(future.globals.maxSize = 1e12)

### function to use

### Main
dataset <- snakemake@params[["dataset"]]

## Import data
seu <- readRDS(snakemake@input[["data"]])

celltypes <- c("CD4 Naive","CD4 TCM","CD4 TEM","CD8 Naive","CD8 TCM","CD8 TEM")

# subset out only interested cell types
seu <- subset(x = seu, subset = (predicted.celltype.l2 %in% celltypes))

# subset out only donors with enough cells
cell_counts <- seu@meta.data %>%
  count(donor_id, predicted.celltype.l2) %>%
  pivot_wider(names_from = "predicted.celltype.l2", values_from = n, values_fill = 0)
print(cell_counts)

eligible_donors <- cell_counts %>%
  filter(if_all(all_of(celltypes), ~ .x >= 3)) %>%
  pull(donor_id)

seu <- subset(seu, subset = donor_id %in% eligible_donors)


# normalize
seu <- NormalizeData(seu)

mat <- AverageExpression(seu, assays = "RNA", group.by = c("donor_id","predicted.celltype.l2"), features = rownames(seu), return.seurat = FALSE)$RNA
mat <- as.data.frame(t(mat)) %>% rownames_to_column("group") %>% separate(group, into = c("donor_id", "celltype"), sep = "_")
print((mat[1:10, 1:3]))

mat_long <- mat %>% 
	pivot_longer(cols = -c(celltype, donor_id), names_to = "gene", values_to = "expression") %>%
  mutate(celltype_gene = paste(celltype, gene, sep = ".")) %>%
  select(donor_id, celltype_gene, expression)
print(head(mat_long))

final_mat <- mat_long %>%
  pivot_wider(names_from = celltype_gene, values_from = expression) %>% as.data.frame() %>% column_to_rownames("donor_id")

print((final_mat[, 1:10]))

## Metadata
meta <- seu@meta.data %>% select(donor_id, age, sex, self_reported_ethnicity) 
meta <- meta[!duplicated(meta$donor_id),] %>% arrange(donor_id)

if(dataset == "onek1k"){
	meta$donor_id <- paste0("g", gsub("_", "-", meta$donor_id))
}else{
	meta$donor_id <- gsub("_", "-", meta$donor_id)
}
rownames(meta) <- meta$donor_id
print(head(meta))

cat("raw and metadata donor id set equal: ", setequal(rownames(meta), rownames(final_mat)), "\n")

final_mat$age <- as.integer(meta[rownames(final_mat), "age"])
final_mat$sex <- meta[rownames(final_mat), "sex"]
final_mat$ethnicity <- meta[rownames(final_mat), "self_reported_ethnicity"]
final_mat$dataset <- dataset
final_mat <- final_mat %>% rownames_to_column("rowname")
print(final_mat[,c("rowname","age","sex","dataset","ethnicity")])
cat("final_mat dim: ", dim(final_mat), "\n")

export(final_mat, snakemake@output[["pb"]])




