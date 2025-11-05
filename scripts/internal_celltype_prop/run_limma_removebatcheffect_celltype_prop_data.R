library(rio)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(dplyr)
library(limma)

data <- import(snakemake@input[["data"]])

## Get all degs from deswan result (all both, gender specific)
limma_data <- data %>% select(-c(donor_id,age,sex,dataset,ethnicity))
print(dim(limma_data))

limma_res <- removeBatchEffect(t(limma_data), data$dataset)

limma_res <- as.data.frame(t(limma_res))
limma_res$donor_id <- data$donor_id
limma_res$age <- data$age
limma_res$sex <- data$sex
limma_res$dataset <- data$dataset
limma_res$ethnicity <- data$ethnicity
export(limma_res, snakemake@output[["res"]])

