print(.libPaths())
.libPaths(c("resources/r_package"))
print(.libPaths())

options(future.globals.maxSize = 1e12)
options(timeout = 240)

library(Seurat)
library(tidyr)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(Azimuth)
library(SeuratData)

seu <- readRDS(snakemake@input[["data"]])
print(seu)

set.seed(123)
annotated <- RunAzimuth(seu, reference = "pbmcref")
print("azimuth annotated")

saveRDS(annotated, snakemake@output[["annotated"]])
print("rds saved")

p1 <-  DimPlot(annotated, group.by = "predicted.celltype.l1", label = FALSE, label.size = 3)
p2 <-  DimPlot(annotated, group.by = "predicted.celltype.l2", label = FALSE, label.size = 3)

ggsave(snakemake@output[["plot"]], ggarrange(p1, p2, nrow = 1, ncol = 2), width = 16, height = 5)

