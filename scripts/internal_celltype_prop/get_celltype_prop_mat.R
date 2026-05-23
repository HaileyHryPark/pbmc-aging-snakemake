library(Seurat)
options(Seurat.object.assay.version = "v5")
library(rio)
library(tidyr)
library(dplyr)
library(tidyverse)
library(ggpubr)

options(future.globals.maxSize = 1e12)

### function to use

## Import data
seu <- readRDS(snakemake@input[["data"]])

## Celltype proportion
countdf1 <- as.data.frame.matrix(table(seu$donor_id, seu$predicted.celltype.l1))
total <- rowSums(countdf1)
propdf1 <- countdf1 / total
colnames(propdf1) <- gsub(pattern = "[ \\-\\.\\:]", "_", x = colnames(propdf1))
print(propdf1[1:5,1:5])

countdf2 <- as.data.frame.matrix(table(seu$donor_id, seu$predicted.celltype.l2))
total <- rowSums(countdf2)
propdf2 <- countdf2 / total
colnames(propdf2) <- gsub(pattern = "[ \\-\\.\\:]", "_", x = colnames(propdf2))
print(propdf2[1:5,1:5])

propdf <- cbind(propdf1, propdf2)
propdf[is.na(propdf)] <- 0
propdf <- propdf * 100
propdf <- propdf %>% as.data.frame()
propdf$donor_id <- rownames(propdf)
print(propdf[1:5,1:5])

meta <- seu@meta.data %>% distinct(donor_id, age, sex, ethnicity = self_reported_ethnicity) %>% mutate(dataset = snakemake@params[["dataset"]])
print(head(meta))

df <- merge(meta, propdf, by = "donor_id", all = TRUE)
print(head(df))

export(df, snakemake@output[["ctp"]])

table <- seu@meta.data %>% select(donor_id, age, sex, predicted.celltype.l1, predicted.celltype.l2) %>%
		mutate(dataset = snakemake@params[["dataset"]])
print(head(table))

export(table, snakemake@output[["meta"]])

table <- seu@meta.data %>% group_by(donor_id, sex, predicted.celltype.l1) %>% summarise(cell_num = n(), total_UMI = sum(nCount_RNA), total_gene = sum(nFeature_RNA), mean_UMI = mean(nCount_RNA), mean_gene = mean(nFeature_RNA)) %>% group_by(donor_id, sex) %>% mutate(proportion = cell_num / sum(cell_num)) %>% ungroup()
print(head(table))

export(table, snakemake@output[["meta2"]])

