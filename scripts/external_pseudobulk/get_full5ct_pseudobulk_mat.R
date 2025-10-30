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

celltypes <- c("CD4 T","CD8 T","NK","B","Mono")

# subset out only interested cell types
seu <- subset(x = seu, subset = (predicted.celltype.l1 %in% celltypes))
seu$sample_id <- paste0("g", gsub("_", "-", seu$sample_id))

# normalize
seu <- NormalizeData(seu)

mat <- AverageExpression(seu, assays = "RNA", group.by = c("sample_id","predicted.celltype.l1"), features = rownames(seu), return.seurat = FALSE)$RNA
mat <- as.data.frame(t(mat)) %>% rownames_to_column("group") %>% separate(group, into = c("sample_id", "celltype"), sep = "_")
print((mat[1:10, 1:3]))

mat_long <- mat %>% 
	pivot_longer(cols = -c(celltype, sample_id), names_to = "gene", values_to = "expression") %>%
  mutate(celltype_gene = paste(celltype, gene, sep = ".")) %>%
  select(sample_id, celltype_gene, expression)
print(head(mat_long))

final_mat <- mat_long %>%
  pivot_wider(names_from = celltype_gene, values_from = expression) %>% as.data.frame() %>% column_to_rownames("sample_id")

print((final_mat[, 1:10]))

## Metadata
meta <- seu@meta.data %>% select(donor_id, sample_id, age, sex, disease, self_reported_ethnicity) 
meta <- meta[!duplicated(meta$sample_id),] %>% arrange(donor_id, sample_id)

rownames(meta) <- meta$sample_id
print(head(meta))

cat("raw and metadata donor id set equal: ", setequal(rownames(meta), rownames(final_mat)), "\n")
cat("raw and metadata donor id identical: ", identical(rownames(meta), rownames(final_mat)), "\n")

final_mat$donor_id <- meta[rownames(final_mat), "donor_id"]
final_mat$age <- as.character(meta[rownames(final_mat), "age"])
final_mat$sex <- meta[rownames(final_mat), "sex"]
final_mat$ethnicity <- meta[rownames(final_mat), "self_reported_ethnicity"]
final_mat$disease <- meta[rownames(final_mat), "disease"]
final_mat$dataset <- dataset
final_mat <- final_mat %>% rownames_to_column("sample_id")
print(final_mat[,c("sample_id","donor_id","age","sex","disease","dataset","ethnicity")])
cat("final_mat dim: ", dim(final_mat), "\n")

export(final_mat, snakemake@output[["pb"]])




