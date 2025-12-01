library(rio)
library(dplyr)
library(tidyverse)
library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

## Main
limma <- import(snakemake@input[["limma"]])

## The annotation object
annot <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19) %>% as.data.frame()

annotated <- limma %>% left_join(annot, by = join_by("probeID" == "Name"))
export(annotated, snakemake@output[["annotated"]])
