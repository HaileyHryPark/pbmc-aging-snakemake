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
  }else if ("sampleID" %in% names(df)) {
	sum <- sum %>% mutate(num_sample = length(unique(df$sampleID)))
  }else {
	sum <- sum %>% mutate(num_sample = num_donor)
  }
  sum <- sum %>% mutate(min_age = min(df$age, na.rm = TRUE), max_age = max(df$age, na.rm = TRUE))

  sum <- cbind(sum, row)
  colnames(sum) <- c("dataset", "num_donor", "num_healthy", "num_disease", "female_prop", "num_sample", "min_age", "max_age", cols)

  sum$others1 <- NA
  sum$others2 <- NA

  if(dataset_name == "Ren et al."){
    sum$others1 <- paste(
      names(table(unique(df[, c("donor_id", "tissue")])$tissue)),
      table(unique(df[, c("donor_id", "tissue")])$tissue),
      sep = ":",
      collapse = "; "
    )
    sum$others2 <- paste(
      names(table(unique(df[, c("donor_id", "Sample time")])$`Sample time`)),
      table(unique(df[, c("donor_id", "Sample time")])$`Sample time`),
      sep = ":",
      collapse = "; "
    )
  }else if(dataset_name == "Wellcome"){
    sum$others2 <- paste(
      names(table(unique(df[, c("donor_id", "Status")])$Status)),
      table(unique(df[, c("donor_id", "Status")])$Status),
      sep = ":",
      collapse = "; "
    )
    sum$others1 <- paste(
      names(table(unique(df[, c("donor_id", "Resample")])$Resample)),
      table(unique(df[, c("donor_id", "Resample")])$Resample),
      sep = ":",
      collapse = "; "
    )
  }else if(dataset_name == "COMBAT"){
    sum$others1 <- paste(
      names(table(unique(df[, c("donor_id", "Source")])$Source)),
      table(unique(df[, c("donor_id", "Source")])$Source),
      sep = ":",
      collapse = "; "
    )
  }

  return(sum)
}

data1 <- import(snakemake@input[["ofcn"]])
data2 <- import(snakemake@input[["ofcp"]])
data3 <- import(snakemake@input[["omcn"]])
data4 <- import(snakemake@input[["omcp"]])
data5 <- import(snakemake@input[["yfcn"]])
data6 <- import(snakemake@input[["yfcp"]])
data7 <- import(snakemake@input[["ymcn"]])
data8 <- import(snakemake@input[["ymcp"]])

# Example combining multiple datasets
summary_all <- bind_rows(
  make_summary_row(data1,"ofcn"), 
  make_summary_row(data2,"ofcp"),
  make_summary_row(data3,"omcn"),
  make_summary_row(data4,"omcp"),
  make_summary_row(data5,"yfcn"),
  make_summary_row(data6,"yfcp"),
  make_summary_row(data7,"ymcn"),
  make_summary_row(data8,"ymcp")
)

export(summary_all, snakemake@output[["table"]])
