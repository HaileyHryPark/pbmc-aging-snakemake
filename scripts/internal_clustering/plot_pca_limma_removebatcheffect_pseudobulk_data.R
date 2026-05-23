library(rio)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(dplyr)
library(limma)
library(svglite)

data <- import(snakemake@input[["data"]])
deg <- import(snakemake@input[["deg"]])
gender <- snakemake@params[["gender"]]

dataset_cols <- c("onek1k" = "#1b9e77", "aida" = "#d95f02", "perez" = "#7570b3", "marina" = "#e7298a")

## Get all degs from deswan result (all both, gender specific)
features_to_include <- unique(deg$variable)

## Run limma
if(gender != "both"){
	data <- data %>% filter(sex == gender)
}
limma_data <- data %>% select(-c(rowname,age,sex,dataset,ethnicity))
print(dim(limma_data))

## PCA BEFORE batch correction
pca_before <- prcomp(limma_data, scale. = TRUE)

pca_before_df <- as.data.frame(pca_before$x[,1:2])
pca_before_df$dataset <- data$dataset
pca_before_df$type <- "Before"
var_before <- summary(pca_before)$importance[2,1:2] * 100

limma_res <- removeBatchEffect(t(limma_data), data$dataset)
limma_res <- t(limma_res)

## PCA AFTER batch correction
pca_after <- prcomp(limma_res, scale. = TRUE)

pca_after_df <- as.data.frame(pca_after$x[,1:2])
pca_after_df$dataset <- data$dataset
pca_after_df$type <- "After"
var_after <- summary(pca_after)$importance[2,1:2] * 100

## Combine
pca_df <- bind_rows(pca_before_df, pca_after_df)

## Plot
p1 <- pca_df %>%
  filter(type == "Before") %>%
  ggplot(aes(PC1, PC2, color = dataset)) +
  geom_point(size = 1, alpha = 0.8) +
  labs(title = "Before batch correction",
	x = paste0("PC1 (", round(var_before[1],1), "%)"),
	y = paste0("PC2 (", round(var_before[2],1), "%)")) +
  scale_color_manual(values = dataset_cols) +
  theme_test() +
  theme(plot.title = element_text(hjust = 0.5))

p2 <- pca_df %>%
  filter(type == "After") %>%
  ggplot(aes(PC1, PC2, color = dataset)) +
  geom_point(size = 1, alpha = 0.8) +
  labs(title = "After batch correction",
	x = paste0("PC1 (", round(var_after[1],1), "%)"),
	y = paste0("PC2 (", round(var_after[2],1), "%)")) +
  scale_color_manual(values = dataset_cols) +
  theme_test() +
  theme(plot.title = element_text(hjust = 0.5))

p <- ggarrange(p1, p2, ncol = 2, nrow = 1)

ggsave(snakemake@output[["plots"]], p, width = 8, height = 3)

