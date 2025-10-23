library(rio)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(svglite)

pred_b <- import(snakemake@input[["pred_b"]]) %>% mutate(color_name = "grey")
pred_f <- import(snakemake@input[["pred_f"]]) %>% mutate(color_name = "#E15566")
pred_m <- import(snakemake@input[["pred_m"]]) %>% mutate(color_name = "#4981BF")

plots <- lapply(list(pred_b, pred_f, pred_m), function(data){
	type = unique(data$color_name)
	ggscatter(data, x = "actual_age", y = "predicted_age", color = type, cor.coef = TRUE, alpha = 0.7, cor.coeff.args = list(method = "pearson", label.sep = "\n"), cor.coef.size = 2) +
		geom_smooth(method = "lm", color = "black", se = F) + theme_test() + theme(legend.position = "none") + xlim(c(0,100)) + ylim(c(0,100))
})

pdf(snakemake@output[["predplot"]], width = 6.5, height = 2)
ggarrange(plotlist=plots, ncol=3, nrow=1)
dev.off()
