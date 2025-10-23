library(rio)
library(dplyr)
library(tidyverse)

data <- import(snakemake@input[["data"]])
deg <- import(snakemake@input[["deg"]])

gender <- unique(deg$gender)
print(gender)

features_to_include <- deg %>% pull(variable) %>% unique()

# Not including ethnicity here because probably not necessary in the downstream
subset <- data %>% select(rowname, age, sex, dataset, all_of(features_to_include))

## Decided to not subset for sex here so that I can do cross test later on, instead I will subset out for sex at the beginning of the cv scripts
#if(gender != "Both"){
#	subset <- subset %>% filter(sex == tolower(gender))
#}

export(subset, snakemake@output[["res"]])
