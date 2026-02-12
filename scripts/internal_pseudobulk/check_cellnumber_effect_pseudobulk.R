library(rio)
library(tidyverse)
library(ggpubr)

non_feature_cols <- c("age", "sex", "dataset", "ethnicity")

avg <- import(snakemake@input[["avg"]]) %>% select(-non_feature_cols) %>% column_to_rownames()
agg <- import(snakemake@input[["agg"]]) %>% select(-non_feature_cols) %>% column_to_rownames()
cellcount <- import(snakemake@input[["cellcount"]])


plots <- lapply(list(avg, agg), function(mat){

	pca <- prcomp(mat, scale. = T)
	pc1 <- pca$x[,1]
	
	pc_df <- data.frame(donor_id = gsub("g", "", gsub("-", "_", rownames(pca$x))), PC1 = pc1, total_expr = rowSums(mat)) %>% 
		left_join(cellcount, by = "donor_id")
	print(head(pc_df))

	p1 <- ggplot(pc_df, aes(n, PC1, color = dataset)) + 
		geom_point(size = 1, alpha = .8) +
		geom_smooth(method = "lm", se = F) +
		scale_x_log10() +
		theme_test(base_size = 13) + 
		labs(x = "Number of cells per sample (log10)", y = "PC1")
	p2 <- ggplot(pc_df, aes(n, total_expr, color = dataset)) + 
		geom_point(size = 1, alpha = .8) +
		geom_smooth(method = "lm", se = F) +
		scale_x_log10() +
		scale_y_log10() +
		theme_test(base_size = 13) + 
		labs(x = "Number of cells per sample (log10)", y = "Total expression per sample (log10)")
	return(list(p1, p2))
})

ggsave(snakemake@output[["plot1"]], ggarrange(plotlist = lapply(plots, `[[`, 1), ncol = 2, nrow = 1), width = 10, height = 3)
ggsave(snakemake@output[["plot2"]], ggarrange(plotlist = lapply(plots, `[[`, 2), ncol = 2, nrow = 1), width = 10, height = 3)


