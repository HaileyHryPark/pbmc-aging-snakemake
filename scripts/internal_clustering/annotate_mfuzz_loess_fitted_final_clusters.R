library(rio)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(colorspace)
library(purrr)

## Functions
celltype_colors <- c(
  "CD4 T" = "#D2533B",
  "CD8 T"    = "#E6974D",
  "NK"= "#73AF68",
  "B"   = "#79629E",
  "Mono" = "#5B83BF"
)

## Import data
clust_df <- import(snakemake@input[["final_cluster_df"]]) %>% 
	separate(feature, into = c("celltype", "gene"), sep = "\\.", remove = F) %>%
	filter(!is.na(merged_cluster))
print(head(clust_df))

gender <- snakemake@params[["gender"]]

if(gender == "both"){
	clust_df <- clust_df %>%
		mutate(final_cluster = ifelse(merged_cluster == "Merged_2", "Early\ndecrease", 
					ifelse(merged_cluster == "Merged_6", "Late\nincrease", 
					ifelse(merged_cluster == "Merged_4", "Continuous\nincrease", 
					ifelse(merged_cluster == "Merged_5", "Continuous\ndecrease", 
					ifelse(merged_cluster == "Merged_1", "Early\nfluctuation", 
					ifelse(merged_cluster == "Merged_3", "Early\nincrease", NA)))))))
}else if(gender == "female"){
	clust_df <- clust_df %>%
		mutate(final_cluster = ifelse(merged_cluster == "Merged_4", "Continuous\ndecrease", 
					ifelse(merged_cluster == "Merged_2", "Early\ndecrease", 
					ifelse(merged_cluster == "Merged_5", "Continuous\nincrease", 
					ifelse(merged_cluster == "Merged_1", "Early\nincrease", 
					ifelse(merged_cluster == "Merged_6", "Inverted\nUshape", 
					ifelse(merged_cluster == "Merged_3", "Late\nincrease", NA)))))))
}else if(gender == "male"){
	clust_df <- clust_df %>%
		mutate(final_cluster = ifelse(merged_cluster == "Merged_2", "Continuous\ndecrease", 
					ifelse(merged_cluster == "Merged_4", "Early\nincrease", 
					ifelse(merged_cluster == "Merged_1", "Early\ndecrease", 
					ifelse(merged_cluster == "Merged_3", "Early\nfluctuation", NA)))))
}

clust_df <- clust_df %>%
        mutate(final_cluster = factor(final_cluster, levels = c("Early\nincrease", "Early\ndecrease", "Continuous\ndecrease", "Early\nfluctuation", "Inverted\nUshape", "Continuous\nincrease", "Late\nincrease")))

print(head(clust_df))
print(table(clust_df$final_cluster))
print(table(clust_df$merged_cluster))

export(clust_df, snakemake@output[["annotated"]])

