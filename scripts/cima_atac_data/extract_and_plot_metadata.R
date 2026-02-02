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

## Plots
p1 <- VlnPlot(pbmc,
  features = c(
    "nCount_peaks",
    "TSS.enrichment",
    "blacklist_ratio",
    "nucleosome_signal",
    "pct_reads_in_peaks"
  ),
  pt.size = 0.1,
  ncol = 5
)
ggsave(snakemake@output[["plot1"]], p1, width = 9, height = 3)

p2 <- DensityScatter(pbmc, x = 'nCount_peaks', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
ggsave(snakemake@output[["plot2"]], p2, width = 7, height = 4)

pbmc$nucleosome_group <- ifelse(pbmc$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
p3 <- FragmentHistogram(object = pbmc, group.by = 'nucleosome_group')
ggsave(snakemake@output[["plot3"]], p3, width = 7, height = 4)


## Metadata
metadata_df <- pbmc@meta.data %>% rownames_to_column("cell_barcode")
export(metadata_df, snakemake@output[["meta"]])

