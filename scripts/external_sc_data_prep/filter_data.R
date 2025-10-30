library(rio)
library(Seurat)
library(dplyr)

seu <- readRDS(snakemake@input[["seu"]])

filtered <- subset(x = seu,
                        subset = (nCount_RNA > 1000)
                        & (nCount_RNA < 20000)
                        & (nFeature_RNA > 300)
                        & (nFeature_RNA < 6000)
                        & (percent_mito < 10))

print(head(filtered@meta.data))
print(rownames(filtered))

saveRDS(filtered, snakemake@output[["filtered"]])
export(filtered@meta.data, snakemake@output[["meta"]])
