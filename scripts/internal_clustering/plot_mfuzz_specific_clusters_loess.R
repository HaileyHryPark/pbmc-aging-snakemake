library(rio)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(colorspace)
library(purrr)

celltype_colors <- c(
  "CD4 T" = "#D2533B",
  "CD8 T"    = "#E6974D",
  "NK"= "#73AF68",
  "B"   = "#79629E",
  "Mono" = "#5B83BF"
)

## Functions
make_smoothed_heatmap_gam_clustered <- function(smooth, features, gender_to_clust, label_features = NULL, title = "") {

  sub <- smooth %>% filter(feature %in% features)

  ## Cluster first by female features
  mat_wide <- sub %>% filter(gender == gender_to_clust) %>% select(-gender) %>% pivot_wider(names_from = age, values_from = fitted) %>%
    column_to_rownames("feature") %>%
    as.matrix()

  # Remove features with NAs (to avoid clustering issues)
  mat_wide <- mat_wide[complete.cases(mat_wide), ]

  # -------------------------------
  # Step 2: Hierarchical clustering
  # -------------------------------
  dist_mat <- dist(mat_wide)
  hc <- hclust(dist_mat)
  feature_order <- rownames(mat_wide)[hc$order]

  # -------------------------------
  # Step 4: Selective y-axis labeling
  # -------------------------------
  # Mark only the label_features (leave others blank)
  #y_labels <- ifelse(levels(sub$feature) %in% label_features,
  #                   levels(sub$feature),
  #                   "")

  # -------------------------------
  # Step 5: Split by sex and plot
  # -------------------------------
  sub_female <- sub %>% filter(gender == "female", feature %in% feature_order) %>% mutate(feature = factor(feature, levels = feature_order))
  sub_male <- sub %>% filter(gender == "male", feature %in% feature_order) %>% mutate(feature = factor(feature, levels = feature_order))

  fill_scale <- scale_fill_continuous_diverging(palette = "Blue-Red", limits = c(-1.5, 1.5), name = "z-score")
  heatmap_theme <- theme_minimal(base_size = 15) +
    theme(panel.grid = element_blank(), 
          axis.ticks.y = element_blank(),
          legend.position = "bottom")

  all_features_ordered <- levels(sub_female$feature)
  y_breaks <- all_features_ordered
  y_labels <- ifelse(all_features_ordered %in% label_features,
                     all_features_ordered, "") # only label selected
  feature_positions <- which(levels(sub_female$feature) %in% label_features)

  p_female <- ggplot(sub_female, aes(x = age, y = feature, fill = fitted)) +
    geom_raster(interpolate = TRUE) +
    geom_hline(yintercept = feature_positions, color = "black", linetype = "dashed", size = 0.3) +
    labs(x = "Age", y = "Features") +
    xlim(19,93)+
    scale_y_discrete(breaks = y_breaks, labels = y_labels) +
    fill_scale + heatmap_theme

  p_male <- ggplot(sub_male, aes(x = age, y = feature, fill = fitted)) +
    geom_raster(interpolate = TRUE) +
    geom_hline(yintercept = feature_positions, color = "black", linetype = "dashed", size = 0.3) +
    labs(x = "Age", y = "Features") +
    xlim(19,93)+
    scale_y_discrete(breaks = y_breaks, labels = y_labels) +
    fill_scale + heatmap_theme

  p_both <- ggarrange(p_female, p_male, ncol = 2, nrow = 1)

  return(p_both)
}


## Import data
clust_df_f <- import(snakemake@input[["clust_df_f"]]) %>% filter(!is.na(final_cluster))
clust_df_m <- import(snakemake@input[["clust_df_m"]]) %>% filter(!is.na(final_cluster))

## Features in specific subsets
subset1 <- intersect(clust_df_f %>% filter(final_cluster == "Continuous\nincrease") %>% pull(feature),
	clust_df_m %>% filter(final_cluster == "Early\nincrease") %>% pull(feature))
subset2 <- intersect(clust_df_f %>% filter(final_cluster == "Early\nincrease") %>% pull(feature), 
	clust_df_m %>% filter(final_cluster == "Irregular\nchange") %>% pull(feature))

## Import fit and spans
mat_f <- import(snakemake@input[["fit_res_f"]]) %>% filter(feature %in% c(subset1, subset2))
span_f <- import(snakemake@input[["span_res_f"]]) %>% filter(feature %in% c(subset1, subset2))
span_f <- span_f %>% group_by(feature) %>% slice_min(rmse, n = 1, with_ties = FALSE)

mat_m <- import(snakemake@input[["fit_res_m"]]) %>% filter(feature %in% c(subset1, subset2))
span_m <- import(snakemake@input[["span_res_m"]]) %>% filter(feature %in% c(subset1, subset2))
span_m <- span_m %>% group_by(feature) %>% slice_min(rmse, n = 1, with_ties = FALSE)

mat_all <- bind_rows(mat_f %>% mutate(gender = "female"), mat_m %>% mutate(gender = "male"))

mat_f <- mat_f %>%
  left_join(span_f %>% select(feature, span), by = "feature")
mat_m <- mat_m %>%
  left_join(span_m %>% select(feature, span), by = "feature")

plotdata_f <- merge(mat_f, clust_df_f, all.x = T, by = "feature")
print(head(plotdata_f))
print(table(plotdata_f$final_cluster))
print(table(plotdata_f$merged_clusters))

plotdata_m <- merge(mat_m, clust_df_m, all.x = T, by = "feature")
print(head(plotdata_m))
print(table(plotdata_m$final_cluster))
print(table(plotdata_m$merged_clusters))

n_features <- plotdata_f %>%
  distinct(final_cluster, feature) %>%
  count(final_cluster, name = "n_features")

cell_prop_summary <- plotdata_f %>%
  mutate(celltype = factor(celltype, levels = c("CD4 T","CD8 T","NK","B","Mono"))) %>%
  distinct(feature, final_cluster, celltype) %>%
  group_by(final_cluster, celltype) %>%
  summarise(count = n(), .groups = 'drop')

print(head(plotdata_f %>%
  mutate(celltype = factor(celltype, levels = c("CD4 T","CD8 T","NK","B","Mono"))) %>%
  distinct(feature, final_cluster, celltype)))

pdf(snakemake@output[["plot1"]], width = 6, height = 6)
lapply(list(subset1, subset2), function(sub){

pf <- ggplot(plotdata_f %>% filter(feature %in% sub), aes(x=age, y=fitted, group=feature)) +
	geom_line(color = "#E15566", size = 0.3, alpha = 0.5) +
	geom_smooth(aes(group = 1), method = "gam", formula = y ~ s(x, k = 10),
              color = "black", size = 2, linetype = "solid") +
	theme_linedraw(base_size = 15)+
	theme(panel.grid = element_blank(), strip.text = element_text(size=17), legend.position = "bottom")+
	ylim(-1.7,1.7)+
	xlim(19,93)+
	labs(x = "Age", y = "Scaled expression")

pm <- ggplot(plotdata_m %>% filter(feature %in% sub), aes(x=age, y=fitted, group=feature)) +
	geom_line(color = "#4981BF", size = 0.3, alpha = 0.5) +
	geom_smooth(aes(group = 1), method = "gam", formula = y ~ s(x, k = 10),
              color = "black", size = 2, linetype = "solid") +
	theme_linedraw(base_size = 15)+
	theme(panel.grid = element_blank(), strip.text = element_text(size=17), legend.position = "bottom")+
	ylim(-1.7,1.7)+
	xlim(19,93)+
	labs(x = "Age", y = "Scaled expression")

p <- ggarrange(pf, pm, ncol = 2, nrow = 1)

ph <- make_smoothed_heatmap_gam_clustered(mat_all, sub, "male")

plot(ggarrange(p, ph, ncol = 1, nrow = 2, heights = c(2,2)))

})

dev.off()

pdf(snakemake@output[["plot2"]], width = 6, height = 3.5)

pie(cell_prop_summary %>% filter(final_cluster == "Continuous\nincrease") %>% pull(count), col = celltype_colors, main = "Subset1", labels = rep("", length(celltype_colors)))
pie(cell_prop_summary %>% filter(final_cluster == "Early\nincrease") %>% pull(count), col = celltype_colors, main = "Subset2", labels = rep("", length(celltype_colors)))

dev.off()
