library(rio)
library(dplyr)
library(tidyverse)
library(purrr)
library(ggpubr)


## Functions
celltypes <- c("CD4_T", "CD8_t", "NK", "B", "Mono")
celltype_col <- c(
  "CD4 T" = "#D2533B",
  "CD8 T"    = "#E6974D",
  "NK"= "#73AF68",
  "B"   = "#79629E",
  "Mono" = "#5B83BF"
)

both <- import(snakemake@input[["degb"]]) %>% distinct(celltype, variable) %>% count(celltype, name = "deg_both_n") 
female <- import(snakemake@input[["degf"]]) %>% distinct(celltype, variable) %>% count(celltype, name = "deg_female_n")
male <- import(snakemake@input[["degm"]]) %>% distinct(celltype, variable) %>% count(celltype, name = "deg_male_n")

deg_summary <- full_join(both, female, by = "celltype") %>%
  full_join(male, by = "celltype") %>%
  replace(is.na(.), 0)
print(head(deg_summary))

table <- import(snakemake@input[["table"]])
print(head(table))

meta_summary_by_sex <- table %>% group_by(sex, predicted.celltype.l1) %>% 
	summarise(total_cells = sum(cell_num, na.rm = T),
		total_UMIs = sum(total_UMI, na.rm=T),
		total_genes = sum(total_gene, na.rm=T),
		mean_tcells = mean(cell_num, na.rm=T),
		mean_tUMIs = mean(total_UMI, na.rm=T),
		mean_tgenes = mean(total_gene, na.rm=T),
		mean_UMIs = mean(mean_UMI, na.rm=T),
		mean_genes = mean(mean_gene, na.rm=T),
		mean_cells = mean(proportion, na.rm=T),
		.groups = "drop"
		)

meta_summary_all <- table %>% group_by(predicted.celltype.l1) %>%
	summarise(total_cells = sum(cell_num, na.rm = T),
		total_UMIs = sum(total_UMI, na.rm=T),
		total_genes = sum(total_gene, na.rm=T),
		mean_tcells = mean(cell_num, na.rm=T),
		mean_tUMIs = mean(total_UMI, na.rm=T),
		mean_tgenes = mean(total_gene, na.rm=T),
		mean_UMIs = mean(mean_UMI, na.rm=T),
		mean_genes = mean(mean_gene, na.rm=T),
		mean_cells = mean(proportion, na.rm=T),
		.groups = "drop"
		) %>% mutate(sex = "both")

meta_summary <- bind_rows(meta_summary_by_sex, meta_summary_all) %>% select(predicted.celltype.l1, sex, everything())
print(head(meta_summary))

meta_summary_wide <- meta_summary %>% pivot_wider(id_cols = predicted.celltype.l1, names_from = sex, values_from = c(total_cells, total_UMIs, total_genes, mean_tcells, mean_tUMIs, mean_tgenes, mean_UMIs, mean_genes, mean_cells), names_glue = "{.value}_{sex}")
print(head(meta_summary_wide))

all_summary <- deg_summary %>% left_join(meta_summary_wide, by = c("celltype" = "predicted.celltype.l1"))
print(head(all_summary))
print(dim(all_summary))

pdf(snakemake@output[["plot"]], width = 30, height = 3)

lapply(list("both","female","male"), function(g){ 
	df <- all_summary 
	plots <- lapply(list("total_cells","total_UMIs","total_genes", "mean_tcells", "mean_tUMIs", "mean_tgenes", "mean_UMIs", "mean_genes", "mean_cells"), function(var){
		p1 <- ggscatter(df, x = paste(var,g, sep = "_"), y = paste("deg",g,"n",sep="_"), color = "celltype", label = "celltype", repel = TRUE, size = 2, add = "reg.line", cor.coef = TRUE, conf.int = FALSE) +
			scale_color_manual(values = celltype_col) + 
			theme_test(base_size = 15) + theme(legend.position = "none")
		return(p1)
	})
	ggarrange(plotlist = plots, ncol = 9, nrow = 1)
})
dev.off()

