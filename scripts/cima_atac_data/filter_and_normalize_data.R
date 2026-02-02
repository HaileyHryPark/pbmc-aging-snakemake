.libPaths("resources/r_package")

library(Signac)
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(tidyverse)
library(rio)

pbmc <- readRDS(snakemake@input[["seu"]])
print(Assays(pbmc))
DefaultAssay(pbmc) <- "peaks"


## QC metrics
pbmc <- NucleosomeSignal(pbmc)
pbmc <- TSSEnrichment(pbmc)

pbmc$pct_reads_in_peaks <- pbmc$peak_region_fragments /
                           pbmc$passed_filters * 100

pbmc$blacklist_ratio <- FractionCountsInRegion(
  pbmc,
  assay = "peaks",
  regions = blacklist_hg38_unified
)

## Filter
pbmc <- subset(
  pbmc,
  subset =
    nCount_peaks > 9000 &
    nCount_peaks < 100000 &
    pct_reads_in_peaks > 40 &
    blacklist_ratio < 0.01 &
    nucleosome_signal < 4 &
    TSS.enrichment > 4
)

## Normalization
pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = "q0")
pbmc <- RunSVD(pbmc)

SaveH5Seurat(pbmc, snakemake@output[["res"]])

## Metadata
qc_cols <- c(
  "nCount_peaks",
  "nFeature_peaks",
  "TSS.enrichment",
  "nucleosome_signal",
  "pct_reads_in_peaks",
  "blacklist_ratio",
  "sample",
  "age",
  "sex"
)

qc_df <- pbmc@meta.data[, qc_cols, drop = FALSE] %>% rownames_to_column("cell_barcode")

export(qc_df, snakemake@output[["meta"]])
