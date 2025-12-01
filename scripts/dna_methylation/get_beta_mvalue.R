library(rio)
library(dplyr)
library(tidyverse)
library(minfi)
library(limma)
library(ENmix)
library(ChAMP)

col_agegroup <- c("#86B875","#4CB9CC","#CD99D8")
col_gender <- c("#E15566","#4981BF")

mSetNoobFlt <- readRDS(snakemake@input[["mset"]])
probes <- import(snakemake@input[["probes"]])

## Get beta values
beta <- getBeta(mSetNoobFlt)

## Merge replicate probes
replicate_probes <- probes[probes$probeID %in% rownames(beta), ] %>%
  add_count(prefix) %>%
  filter(n > 1)

unique_probes <- probes[probes$probeID %in% rownames(beta), ] %>%
  add_count(prefix) %>%
  filter(n == 1)

## Beta data for unique probes
beta_unique <- beta[unique_probes$probeID,] %>% as.data.frame()
rownames(beta_unique) <- sapply(strsplit(rownames(beta_unique),"_"), `[`, 1)

## Beta data for replicate probes - used mean value to merge
beta_replicate <- beta[replicate_probes$probeID,] %>% as.data.frame()
beta_replicate$probe_prefix <- sapply(strsplit(rownames(beta_replicate),"_"), `[`, 1)
beta_replicate_merged <- beta_replicate %>%
  group_by(probe_prefix) %>%
  mutate_all(funs(mean)) %>%
  distinct %>%
  column_to_rownames("probe_prefix")

identical(colnames(beta_unique), colnames(beta_replicate_merged))

beta_merged <- rbind(beta_unique, beta_replicate_merged)

## BMIQ Correction
beta_normalized <- champ.norm(beta = as.matrix(beta_merged),
                        method = "BMIQ", arraytype ="EPICv2")

beta_normalized <- beta_normalized %>%
  as.data.frame() %>%
  rownames_to_column("probe_prefix")

export(beta_normalized, snakemake@output[["beta"]])

beta_normalized <- beta_normalized %>% column_to_rownames("probe_prefix")

## Convert BMIQ corrected beta values to M value
Mvalue_normalized <- B2M(beta_normalized)

Mvalue_normalized <- Mvalue_normalized %>%
  as.data.frame() %>%
  rownames_to_column("probe_prefix")

export(Mvalue_normalized, snakemake@output[["mvalue"]])
