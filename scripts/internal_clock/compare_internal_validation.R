library(rio)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(svglite)
library(Metrics)

pred_b <- import(snakemake@input[["pred_b"]]) %>% mutate(color_name = "grey", model = "Model B")
pred_f <- import(snakemake@input[["pred_f"]]) %>% mutate(color_name = "#E15566", model = "Model F")
pred_m <- import(snakemake@input[["pred_m"]]) %>% mutate(color_name = "#4981BF", model = "Model M")

plots <- lapply(list(pred_b, pred_f, pred_m), function(data){
	type = unique(data$color_name)
	ggscatter(data, x = "actual_age", y = "predicted_age", color = type, cor.coef = TRUE, alpha = 0.7, cor.coeff.args = list(method = "pearson", label.sep = "\n"), cor.coef.size = 2) +
		geom_smooth(method = "lm", color = "black", se = F) + theme_test() + theme(legend.position = "none") + xlim(c(0,100)) + ylim(c(0,100))
})

pdf(snakemake@output[["predplot"]], width = 6.5, height = 2)
ggarrange(plotlist=plots, ncol=3, nrow=1)
dev.off()

female_res <- bind_rows(pred_b %>% filter(sex == "female"), pred_f)
male_res <- bind_rows(pred_b %>% filter(sex == "male"), pred_m)

female_res_met <- female_res %>% group_by(model, fold) %>% 
	summarise(RMSE = rmse(actual_age, predicted_age), MAE = mae(actual_age, predicted_age), .groups = "drop")
male_res_met <- male_res %>% group_by(model, fold) %>% 
	summarise(RMSE = rmse(actual_age, predicted_age), MAE = mae(actual_age, predicted_age), .groups = "drop")
print(head(female_res_met))
print(head(male_res_met))

pdf(snakemake@output[["resplot"]], width = 3, height = 2)

p1 <- ggplot(female_res_met, aes(x = model, y = RMSE)) + 
	geom_boxplot(aes(color = model), width = 0.6) +
	scale_color_manual(values = c("Model B" = "black", "Model F" = "#E15566")) +
	geom_point(position = position_jitter(width = 0.2), size = 1, alpha = 0.8, color = "#E15566") +
	stat_compare_means(method = "wilcox", paired = FALSE, label = "p", comparisons = list(c("Model B", "Model F")), , tip.length = 0, bracket.size = 0.7, vjust = 0.1) + 
	ylim(7.8, 13)+
	xlab("")+
	theme_classic(base_size = 15)+
	theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))

p2 <- ggplot(male_res_met, aes(x = model, y = RMSE)) + 
	geom_boxplot(aes(color = model), width = 0.6) +
	scale_color_manual(values = c("Model B" = "black", "Model M" = "#4981BF")) +
	geom_point(position = position_jitter(width = 0.2), size = 1, alpha = 0.8, color = "#4981BF") +
	stat_compare_means(method = "wilcox", paired = FALSE, label = "p", comparisons = list(c("Model B", "Model M")), , tip.length = 0, bracket.size = 0.7, vjust = 0.1) + 
	ylim(5.8, 13)+
	xlab("")+
	theme_classic(base_size = 15)+
	theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))

ggarrange(p1, p2, ncol = 2, nrow = 1)

p1 <- ggplot(female_res_met, aes(x = model, y = MAE)) + 
	geom_boxplot(aes(color = model), width = 0.6) +
	scale_color_manual(values = c("Model B" = "black", "Model F" = "#E15566")) +
	geom_point(position = position_jitter(width = 0.2), size = 1, alpha = 0.8, color = "#E15566") +
	stat_compare_means(method = "wilcox", paired = FALSE, label = "p", comparisons = list(c("Model B", "Model F")), , tip.length = 0, bracket.size = 0.7, vjust = -0.4) + 
	ylim(5.9, 10.5)+
	xlab("")+
	theme_classic(base_size = 15)+
	theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))

p2 <- ggplot(male_res_met, aes(x = model, y = MAE)) + 
	geom_boxplot(aes(color = model), width = 0.6) +
	scale_color_manual(values = c("Model B" = "black", "Model M" = "#4981BF")) +
	geom_point(position = position_jitter(width = 0.2), size = 1, alpha = 0.8, color = "#4981BF") +
	stat_compare_means(method = "wilcox", paired = FALSE, label = "p", comparisons = list(c("Model B", "Model M")), , tip.length = 0, bracket.size = 0.7, vjust = -0.4) + 
	ylim(5.9, 10.5)+
	xlab("")+
	theme_classic(base_size = 15)+
	theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))

ggarrange(p1, p2, ncol = 2, nrow = 1)

dev.off()



