library(dplyr)
library(tidyverse)
library(ggpubr)
library(AnnotationHub)
library(org.Hs.eg.db)
library(ensembldb)
library(rio)



data <- import(snakemake@input[["data"]])
features <- colnames(data %>% dplyr::select(-c(rowname,age,sex,dataset,ethnicity)))

universe <- data.frame(feature = features)
universe <- universe %>% separate(feature, into = c("celltype", "gene"), sep = "\\.", remove = F)
annot <- ensembldb::select(org.Hs.eg.db, keys = unique(universe$gene), keytype = "SYMBOL", columns = c("SYMBOL","ENTREZID"))
print(head(annot))
universe <- merge(universe, annot, all.x = T, all.y = F, by.x = "gene", by.y = "SYMBOL")
print(head(universe))

export(universe, snakemake@output[["uni"]])

