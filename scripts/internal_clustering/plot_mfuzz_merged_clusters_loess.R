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
gender <- snakemake@params[["gender"]]
mat <- import(snakemake@input[["fit_res"]])

span <- import(snakemake@input[["span_res"]])
best_span <- span %>%
  group_by(feature) %>%
  slice_min(rmse, n = 1, with_ties = FALSE)
mat <- mat %>%
  left_join(best_span %>% select(feature, span), by = "feature")

clust_df <- import(snakemake@input[["var_cluster_df"]]) %>% 
	separate(variable, into = c("celltype", "gene"), sep = "\\.", remove = F) %>% 
	rename(feature=variable)
print(head(clust_df))

if(gender == "both"){
	clust_df <- clust_df %>%
		mutate(final_cluster = ifelse(merged_clusters == "MergedCluster_1", "Early\nincrease", ifelse(merged_clusters == "MergedCluster_2", "Early\ndecrease", ifelse(merged_clusters == "MergedCluster_3", "Continuous\ndecrease", ifelse(merged_clusters == "MergedCluster_4", "Irregular\nchange", "Late\nincrease")))))
}else if(gender == "female"){
	clust_df <- clust_df %>%
		mutate(final_cluster = ifelse(merged_clusters == "MergedCluster_1", "Continuous\ndecrease", ifelse(merged_clusters == "MergedCluster_2", "Early\ndecrease", ifelse(merged_clusters == "MergedCluster_3", "Early\nincrease", ifelse(merged_clusters == "MergedCluster_4", "Continuous\nincrease", ifelse(merged_clusters == "MergedCluster_5", "Inverted\nU-shape", "Late\nincrease"))))))
}else if(gender == "male"){
	clust_df <- clust_df %>%
		mutate(final_cluster = ifelse(merged_clusters == "MergedCluster_1", "Early\ndecrease", ifelse(merged_clusters == "MergedCluster_2", "Early\nincrease", ifelse(merged_clusters == "MergedCluster_3", "Continuous\ndecrease", "Irregular\nchange"))))
}

clust_df <- clust_df %>%
        mutate(final_cluster = factor(final_cluster, levels = c("Early\nincrease", "Early\ndecrease", "Continuous\ndecrease", "Irregular\nchange", "Late\nincrease", "Continuous\nincrease", "Inverted\nU-shape")))

print(head(clust_df))
print(table(clust_df$final_cluster))
print(table(clust_df$merged_clusters))

export(clust_df, snakemake@output[["annotated"]])

plotdata <- merge(mat, clust_df, all.x = T, by = "feature")
print(head(plotdata))
print(table(plotdata$final_cluster))
print(table(plotdata$merged_clusters))

n_features <- plotdata %>%
  distinct(final_cluster, feature) %>%
  count(final_cluster, name = "n_features")

cell_prop_summary <- plotdata %>%
  mutate(celltype = factor(celltype, levels = c("CD4 T","CD8 T","NK","B","Mono"))) %>%
  distinct(feature, final_cluster, celltype) %>%
  group_by(final_cluster, celltype) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(final_cluster) %>%
  mutate(prop = count / sum(count))

print(head(plotdata %>%
  mutate(celltype = factor(celltype, levels = c("CD4 T","CD8 T","NK","B","Mono"))) %>%
  distinct(feature, final_cluster, celltype)))
# Normalize proportions to match the x-scale (e.g., age range)
age_min <- min(plotdata$age)
age_max <- max(plotdata$age)
print(age_min)
print(age_max)

cell_prop_summary <- cell_prop_summary %>%
  mutate(
    xmin = age_min + (age_max - age_min) * cumsum(lag(prop, default = 0)),
    xmax = age_min + (age_max - age_min) * cumsum(prop),
    ymin = 1.7,
    ymax = 1.9
  )
print(head(cell_prop_summary))

pdf(snakemake@output[["plot"]], width = 15, height = 4)
ggplot(plotdata, aes(x=age, y=fitted, group=feature)) +
	geom_line(color = ifelse(gender == "both", "grey", ifelse(gender == "female", "#E15566", "#4981BF")), size = 0.3, alpha = 0.5) +
	geom_smooth(aes(group = 1), method = "gam", formula = y ~ s(x, k = 10),
              color = "black", size = 2, linetype = "solid") +
	geom_text(
	    data = n_features,
	    aes(x = Inf, y = -Inf, label = paste0("n=", n_features)),
	    hjust = 1.1, vjust = -0.5, size = 5, inherit.aes = FALSE
	) +
	geom_rect(
	    data = cell_prop_summary,
	    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = celltype),
	    inherit.aes = FALSE
	) +
	scale_fill_manual(values = celltype_colors) + 
	facet_wrap(~final_cluster, nrow = 1, ncol = 7, drop = FALSE)+
	theme_linedraw(base_size = 16)+
	theme(panel.grid = element_blank(), strip.text = element_text(size=17), legend.position = "bottom")+
	ylim(-1.7,1.9)+
	xlim(19,93)+
	labs(x = "Age", y = "Scaled expression")

dev.off()
