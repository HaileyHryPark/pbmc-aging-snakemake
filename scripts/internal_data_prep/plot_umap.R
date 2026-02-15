library(Seurat)
library(rio)
library(tidyr)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(purrr)
library(svglite)

set.seed(123)
celltype_cols <- c("CD4 T" = "#D2533B", "CD8 T" = "#E6974D", "B" = "#79629E", "NK" = "#73AF68", "Mono" = "#5B83BF", "DC" = "grey88", "other T" = "grey88", "other" = "grey88")

group_var <- "predicted.celltype.l1"

all_embeddings <- list()
all_meta <- list()

onek1k <- snakemake@input[["onek1k"]]
aida <- snakemake@input[["aida"]]
perez <- snakemake@input[["perez"]]
marina <- snakemake@input[["marina"]]

objs <- c(onek1k, aida, perez, marina)

## Loop
for (i in seq_along(objs)) {
  f <- objs[i]
  cat("Loading:", f, "\n")

  obj <- readRDS(f)

  if (!"ref.umap" %in% Reductions(obj)) {
    stop(paste("Missing ref.umap in", f))
  }

  emb <- Embeddings(obj, reduction = "ref.umap")
  meta <- obj@meta.data

  if (!group_var %in% colnames(meta)) {
    stop(paste("Metadata column not found:", group_var, "in", f))
  }

  df <- as.data.frame(emb)
  df$cell_id <- rownames(df)
  df$celltype <- meta[[group_var]]

  all_embeddings[[i]] <- df

  rm(obj, emb, meta, df)
  gc()
}

## Merge
cat("Merging embeddings...\n")

umap_df <- bind_rows(all_embeddings)

rm(all_embeddings)
gc()

export(umap_df, snakemake@output[["coords"]])
umap_df <- umap_df %>% mutate(celltype = factor(celltype, levels = rev(c("NK", "CD4 T", "CD8 T", "Mono", "B", "DC", "other T", "other")))) %>% arrange(celltype)

## Plot
p <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = celltype)) +
  geom_point(size = 0.01, alpha = 1) +
  scale_color_manual(values = celltype_cols) +
  theme_void()+
  theme(legend.position = "none", axis.text = element_blank(), axis.ticks = element_blank())

ggsave(snakemake@output[["plot1"]], p, width = 2.5, height = 2.4)
ggsave(snakemake@output[["plot3"]], p, width = 2.5, height = 2.4, dpi = 1200)

## Downsample
MAX_CELLS <- 1e6

if (nrow(umap_df) > MAX_CELLS) {
  set.seed(123)
  umap_df <- umap_df[sample(seq_len(nrow(umap_df)), MAX_CELLS), ]
  cat("Downsampled to", MAX_CELLS, "cells\n")
}

p <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = celltype)) +
  geom_point(size = 0.01, alpha = 1) +
  scale_color_manual(values = celltype_cols) +
  theme_void()+
  theme(legend.position = "none", axis.text = element_blank(), axis.ticks = element_blank())

ggsave(snakemake@output[["plot2"]], p, width = 2.5, height = 2.4)
ggsave(snakemake@output[["plot4"]], p, width = 2.5, height = 2.4, dpi = 1200)


