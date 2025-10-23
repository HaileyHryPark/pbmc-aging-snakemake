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
		num_healthy = df %>% filter(disease == "normal") %>% pull(donor_id) %>% unique() %>% length())
  sum <- sum %>% mutate(num_disease = num_donor - num_healthy, female_prop = df %>% filter(sex == "female") %>% pull(donor_id) %>% unique() %>% length() / num_donor * 100)
  if ("sample_id" %in% names(df)) {
	sum <- sum %>% mutate(num_sample = length(unique(df$sample_id)))
  } else {
	sum <- sum %>% mutate(num_sample = num_donor)
  }
  sum <- sum %>% mutate(min_age = min(df$age, na.rm = TRUE), max_age = max(df$age, na.rm = TRUE))

  sum <- cbind(sum, row)
  colnames(sum) <- c("dataset", "num_donor", "num_healthy", "num_disease", "female_prop", "num_sample", "min_age", "max_age", cols)
  return(sum)
}

onek1k <- import(snakemake@input[["onek1k"]])
aida <- import(snakemake@input[["aida"]])
marina <- import(snakemake@input[["marina"]])
perez <- import(snakemake@input[["perez"]])

# Example combining multiple datasets
summary_all <- rbind(
  make_summary_row(onek1k, "OneK1K"),
  make_summary_row(aida, "AIDA"),
  make_summary_row(marina, "M et al."),
  make_summary_row(perez, "Perez et al.")
)

export(summary_all, snakemake@output[["table"]])
