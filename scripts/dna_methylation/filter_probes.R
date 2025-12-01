library(rio)
library(dplyr)
library(tidyverse)
library(minfi)
library(limma)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)

col_agegroup <- c("#86B875","#4CB9CC","#CD99D8")
col_gender <- c("#E15566","#4981BF")

mSetNoob <- readRDS(snakemake@input[["mset"]])
print(mSetNoob)

## Probe dataframe
probes <- data.frame(probeID = rownames(mSetNoob))
probes <- probes %>% separate(probeID, c("prefix","suffix"), sep = "_", remove = FALSE) %>% 
  separate(suffix, c("space","strand","deamination","chemistry","index"), sep = "", remove = FALSE) %>%
  dplyr::select(-space)
probes$SNP <- FALSE
probes$Cross <- FALSE
probes$detection <- FALSE
probes$multimap <- FALSE
print(head(probes))

## SNP associated probes
gmSetNoob <- mapToGenome(mSetNoob)
gmSetNoob_noSNP <- dropLociWithSnps(gmSetNoob)
SNPprobes <- setdiff(rownames(gmSetNoob),rownames(gmSetNoob_noSNP))

probes[probes$probeID %in% SNPprobes, "SNP"] <- TRUE 

## Cross-reactive probes / non-specific probes
nsref <- import(snakemake@input[["nsref"]])
probes[probes$prefix %in% nsref$TargetID, "Cross"] <- TRUE 

## Multimap probes
mmref <- read_lines(snakemake@input[["mmref"]])
probes[probes$prefix %in% mmref, "multimap"] <- TRUE 

## detection p value > 0.01
rgSet <- readRDS(snakemake@input[["rgset"]])
detP <- detectionP(rgSet)
detection_probes <- rownames(detP[rowSums(detP < 0.01) == ncol(mSetNoob),])
probes[!probes$probeID %in% detection_probes, "detection"] <- TRUE

## Filtering
mSetNoobFlt <- mSetNoob[probes[!probes$SNP & 
                                 !probes$Cross &
                                 !probes$detection &
                                 !probes$multimap, "probeID"],]

saveRDS(mSetNoobFlt, snakemake@output[["filtered"]])
export(probes, snakemake@output[["probes"]])
