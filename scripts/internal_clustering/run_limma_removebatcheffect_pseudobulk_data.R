library(rio)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(dplyr)
library(limma)

data <- import(snakemake@input[["data"]])
deg <- import(snakemake@input[["deg"]])
gender <- snakemake@params[["gender"]]

## Get all degs from deswan result (all both, gender specific)
features_to_include <- unique(deg$variable)

## Run limma
if(gender != "both"){
	data <- data %>% filter(sex == gender)
}
limma_data <- data %>% select(-c(rowname,age,sex,dataset,ethnicity))
print(dim(limma_data))

limma_res <- removeBatchEffect(t(limma_data), data$dataset)

limma_res <- as.data.frame(t(limma_res))
limma_res <- limma_res[, features_to_include]
limma_res$rowname <- data$rowname
limma_res$age <- data$age
limma_res$sex <- data$sex
limma_res$dataset <- data$dataset
limma_res$ethnicity <- data$ethnicity
export(limma_res, snakemake@output[["res"]])

