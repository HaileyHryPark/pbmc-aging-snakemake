library(rio)
library(tidyr)
library(dplyr)
library(tidyverse)
library(ggpubr)

## Functions
meta <- import(snakemake@input[["data"]]) %>% select(sample_id, donor_id, age, sex, disease, dataset, ethnicity)

lines <- lapply(as.list(unique(meta$dataset)), function(ds){
	m <- meta %>% filter(dataset == ds) %>% arrange(age)
	print(head(m))
	return(c(ds, "\n", paste(m$age, m$ethnicity, m$disease)))
})
writeLines(unlist(lines), snakemake@output[["age_res"]])
