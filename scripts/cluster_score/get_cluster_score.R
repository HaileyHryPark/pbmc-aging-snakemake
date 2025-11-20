library(rio)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(ggfortify)
library(singscore)
library(purrr)

## Function
getClusterScores <- function(mat, cl_df){

	clusters <- unique(cl_df$final_cluster)
	cluster_gs <- lapply(as.list(clusters), function(clust){
		gs <- cl_df %>% filter(final_cluster == clust) %>% pull(feature) %>% unique()
		score <- simpleScore(rankData, upSet = gs) %>% mutate(cluster = clust) %>% 
				select(cluster, score=TotalScore) %>%
				rownames_to_column("sample_id")
		return(score)
	})
	score_df <- bind_rows(cluster_gs)
	print(head(score_df))

	# Flip score direction to be aligned so that higher values = older age
#	for (cl in colnames(score_df)) {
#	  print(cl)
#	  cor_age <- cor(score_df[, cl], score_df$age, use = "pairwise.complete.obs")
#
#	  if (!is.na(cor_age) && cor_age < 0) {
#	    score_df[, cl] <- -1 * score_df[, cl]
#	  }
#	}
	return(score_df)
}

cohort <- snakemake@params[["cohort"]]
mat <- import(snakemake@input[["data"]])
clust_b <- import(snakemake@input[["clust_b"]]) %>% filter(final_cluster != "", !is.na(final_cluster))
clust_f <- import(snakemake@input[["clust_f"]]) %>% filter(final_cluster != "", !is.na(final_cluster))
clust_m <- import(snakemake@input[["clust_m"]]) %>% filter(final_cluster != "", !is.na(final_cluster))

print(mat[1:5,1:5])

if(cohort == "internal"){
	mat <- mat %>% column_to_rownames("rowname") %>% select(-c(age, dataset, ethnicity, sex)) %>% t()
}else{
	mat <- mat %>% column_to_rownames("sample_id") %>% select(-c(age, dataset, disease, donor_id, ethnicity, sex)) %>% t()
}

rankData <- rankGenes(mat)

dflist <- list(clust_b, clust_f, clust_m)
names(dflist) <- c("both", "female", "male")

res <- lapply(as.list(names(dflist)), function(n){

	return(getClusterScores(mat, dflist[[n]]) %>% mutate(gender = n))

})

export(bind_rows(res), snakemake@output[["scores"]])

