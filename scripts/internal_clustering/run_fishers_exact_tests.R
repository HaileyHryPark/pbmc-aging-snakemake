library(rio)
library(tidyverse)


clust_f <- import(snakemake@input[["clust_f"]])
gs <- readRDS(snakemake@input[["gs"]])
ias_gs <- gs$Serrano_Iron_Accumulation_Geneset

sink(snakemake@output[["res1"]])
# Gene level
print("Gene level")
a <- length(intersect(clust_f %>% filter(final_cluster == "Late\nincrease") %>% pull(gene) %>% unique(), ias_gs))
b <-  length(clust_f %>% filter(final_cluster == "Late\nincrease") %>% pull(gene) %>% unique()) - a
c <- length(ias_gs)- a
d <- length(clust_f %>% pull(gene) %>% unique()) - a - b - c

mat <- matrix(c(a, b, c, d), nrow = 2)
fisher.test(mat, alternative = "greater")

# Feature level
print("Feature level")
ias_features <- clust_f %>% filter(gene %in% ias_gs) %>% pull(feature)
print(ias_features)

a <- length(intersect(clust_f %>% filter(final_cluster == "Late\nincrease") %>% pull(feature), ias_features))
b <-  length(clust_f %>% filter(final_cluster == "Late\nincrease") %>% pull(feature)) - a
c <- length(ias_features) - a
d <- nrow(clust_f) - a - b - c

mat <- matrix(c(a, b, c, d), nrow = 2)
fisher.test(mat, alternative = "greater")

sink()
