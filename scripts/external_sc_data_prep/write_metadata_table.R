library(rio)
library(dplyr)
library(tidyverse)

## Function
make_summary_row <- function(df, dataset_name) {
  cols <- c("nCount_RNA", "nFeature_RNA", "percent_mito")
  row <- as.data.frame(
    t(
      sapply(cols, function(x) paste(capture.output(summary(df[[x]])), collapse = "; "))
    )
  ) 

  sum <- data.frame(dataset = dataset_name, 
		num_donor = length(unique(df$donor_id)), 
		num_control = df %>% filter(disease == "CT") %>% pull(donor_id) %>% unique() %>% length())
  sum <- sum %>% mutate(num_sc = num_donor - num_control, female_prop = df %>% filter(sex == "female") %>% pull(donor_id) %>% unique() %>% length() / num_donor * 100, num_sample = length(unique(df$sample_id)))
  sum <- sum %>% mutate(min_age = "50s", max_age = "110s")

  sum <- cbind(sum, row)
  colnames(sum) <- c("dataset", "num_donor", "num_control", "num_sc", "female_prop", "num_sample", "min_age", "max_age", cols)

  return(sum)
}

meta <- import(snakemake@input[["meta"]])

export(make_summary_row(meta, "Supercentenarian"), snakemake@output[["table"]])
