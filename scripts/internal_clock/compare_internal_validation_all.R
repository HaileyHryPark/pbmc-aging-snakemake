library(rio)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(svglite)
library(Metrics)

pred_b1 <- import(snakemake@input[["pred_b1"]]) %>% mutate(color_name = "grey", model = "Model B", alg = "MLP")
pred_f1 <- import(snakemake@input[["pred_f1"]]) %>% mutate(color_name = "#E15566", model = "Model F", alg = "MLP")
predct_f1 <- import(snakemake@input[["predct_f1"]]) %>% mutate(color_name = "#E15566", model = "Model F", alg = "MLP")
pred_m1 <- import(snakemake@input[["pred_m1"]]) %>% mutate(color_name = "#4981BF", model = "Model M", alg = "MLP")
predct_m1 <- import(snakemake@input[["predct_m1"]]) %>% mutate(color_name = "#E15566", model = "Model M", alg = "MLP")
pred_b2 <- import(snakemake@input[["pred_b2"]]) %>% mutate(color_name = "grey", model = "Model B", alg = "Elastic Net")
pred_f2 <- import(snakemake@input[["pred_f2"]]) %>% mutate(color_name = "#E15566", model = "Model F", alg = "Elastic Net")
predct_f2 <- import(snakemake@input[["predct_f2"]]) %>% mutate(color_name = "#E15566", model = "Model F", alg = "Elastic Net")
pred_m2 <- import(snakemake@input[["pred_m2"]]) %>% mutate(color_name = "#4981BF", model = "Model M", alg = "Elastic Net")
predct_m2 <- import(snakemake@input[["predct_m2"]]) %>% mutate(color_name = "#4981BF", model = "Model M", alg = "Elastic Net")
pred_b3 <- import(snakemake@input[["pred_b3"]]) %>% mutate(color_name = "grey", model = "Model B", alg = "XGBoost")
pred_f3 <- import(snakemake@input[["pred_f3"]]) %>% mutate(color_name = "#E15566", model = "Model F", alg = "XGBoost")
predct_f3 <- import(snakemake@input[["predct_f3"]]) %>% mutate(color_name = "#E15566", model = "Model F", alg = "XGBoost")
pred_m3 <- import(snakemake@input[["pred_m3"]]) %>% mutate(color_name = "#4981BF", model = "Model M", alg = "XGBoost")
predct_m3 <- import(snakemake@input[["predct_m3"]]) %>% mutate(color_name = "#4981BF", model = "Model M", alg = "XGBoost")

preds <- bind_rows(pred_m1, pred_m2, pred_m3, pred_b1, pred_b2, pred_b3, pred_f1, pred_f2, pred_f3, predct_f1, predct_f2, predct_f3, predct_m1, predct_m2, predct_m3)
preds <- preds %>% mutate(agediff = predicted_age - actual_age)

female_preds <- preds %>% filter(sex == "female")
male_preds <- preds %>% filter(sex == "male")

both_met <- preds %>% filter(model == "Model B") %>% group_by(alg, model, fold) %>% 
	summarise(RMSE = rmse(actual_age, predicted_age), MAE = mae(actual_age, predicted_age), r = cor(actual_age, predicted_age), .groups = "drop") %>% 
	mutate(alg = factor(alg, levels = c("Elastic Net", "XGBoost", "MLP")), type = "all donors")
female_met <- female_preds %>% group_by(alg, model, fold) %>% 
	summarise(RMSE = rmse(actual_age, predicted_age), MAE = mae(actual_age, predicted_age), r = cor(actual_age, predicted_age), .groups = "drop") %>% 
	mutate(alg = factor(alg, levels = c("Elastic Net", "XGBoost", "MLP")), type = "female donors")
male_met <- male_preds %>% group_by(alg, model, fold) %>% 
	summarise(RMSE = rmse(actual_age, predicted_age), MAE = mae(actual_age, predicted_age), r = cor(actual_age, predicted_age), .groups = "drop") %>%
	mutate(alg = factor(alg, levels = c("Elastic Net", "XGBoost", "MLP")), type = "male donors")
print(head(female_met))
print(head(male_met))

met_all <- bind_rows(both_met, female_met, male_met)
head(met_all)
met_all <- met_all %>% select(alg, model, type, fold, RMSE, MAE, r)
export(met_all, snakemake@output[["fold"]])

met_all_sum <- met_all %>% group_by(alg, model, type) %>% 
	summarise(median_RMSE = median(RMSE), median_MAE = median(MAE), median_r = median(r), .group = "drop")
head(met_all_sum)
export(met_all_sum, snakemake@output[["sum"]])

pdf(snakemake@output[["plots"]], width = 9, height = 3.5)

p1 <- ggplot(female_preds, aes(x = model, y = abs(agediff))) + 
	geom_point(position = position_jitter(width = 0.2), size = 1, alpha = 0.3, color = "#E15566") +
	geom_boxplot(aes(color = model), width = 0.6) +
	scale_color_manual(values = c("Model B" = "black", "Model F" = "#E15566", "Model M" = "#4981BF")) +
	stat_compare_means(method = "wilcox", paired = FALSE, label = "p", comparisons = list(c("Model B", "Model F")), tip.length = 0, bracket.size = 0.7, vjust = -0.4, step.increase = 0.4) + 
	facet_wrap(~alg) +
	xlab("")+
	scale_y_continuous(expand = expansion(mult = c(0,.3)))+
	theme_linedraw(base_size = 15)+
	theme(panel.grid = element_blank(), legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))

p2 <- ggplot(male_preds, aes(x = model, y = abs(agediff))) + 
	geom_point(position = position_jitter(width = 0.2), size = 1, alpha = 0.3, color = "#4981BF") +
	geom_boxplot(aes(color = model), width = 0.6) +
	scale_color_manual(values = c("Model B" = "black", "Model F" = "#E15566", "Model M" = "#4981BF")) +
	stat_compare_means(method = "wilcox", paired = FALSE, label = "p", comparisons = list(c("Model B", "Model M")), tip.length = 0, bracket.size = 0.7, vjust = -0.4, step.increase = 0.4) + 
	facet_wrap(~alg) +
	xlab("")+
	scale_y_continuous(expand = expansion(mult = c(0,.3)))+
	theme_linedraw(base_size = 15)+
	theme(panel.grid = element_blank(), legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))

ggarrange(p1, p2, ncol = 2, nrow = 1)

p1 <- ggplot(female_met, aes(x = model, y = RMSE)) + 
	geom_boxplot(aes(color = model), width = 0.6) +
	scale_color_manual(values = c("Model B" = "grey", "Model F" = "#E15566", "Model M" = "#4981BF")) +
	geom_point(position = position_jitter(width = 0.2), size = 0.5, alpha = 1, color = "black") +
	stat_compare_means(method = "wilcox", paired = FALSE, label = "p", comparisons = list(c("Model B", "Model F"), c("Model M", "Model F")), tip.length = 0, bracket.size = 0.7, vjust = -0.4, step.increase = 0.2) + 
	facet_wrap(~alg) +
	xlab("")+
	scale_y_continuous(expand = expansion(mult = c(0.05,.3)))+
	theme_linedraw(base_size = 15)+
	theme(panel.grid = element_blank(), legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))

p2 <- ggplot(male_met %>% mutate(model = factor(model, levels = c("Model B", "Model M", "Model F"))), aes(x = model, y = RMSE)) + 
	geom_boxplot(aes(color = model), width = 0.6) +
	scale_color_manual(values = c("Model B" = "grey", "Model F" = "#E15566", "Model M" = "#4981BF")) +
	geom_point(position = position_jitter(width = 0.2), size = 0.5, alpha = 1, color = "black") +
	stat_compare_means(method = "wilcox", paired = FALSE, label = "p", comparisons = list(c("Model B", "Model M"), c("Model M", "Model F")), tip.length = 0, bracket.size = 0.7, vjust = -0.4, step.increase = 0.2) + 
	facet_wrap(~alg) +
	xlab("")+
	scale_y_continuous(expand = expansion(mult = c(0.05,.3)))+
	theme_linedraw(base_size = 15)+
	theme(panel.grid = element_blank(), legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))

ggarrange(p1, p2, ncol = 2, nrow = 1)

p1 <- ggplot(female_met, aes(x = model, y = MAE)) + 
	geom_boxplot(aes(color = model), width = 0.6) +
	scale_color_manual(values = c("Model B" = "grey", "Model F" = "#E15566", "Model M" = "#4981BF")) +
	geom_point(position = position_jitter(width = 0.2), size = 0.5, alpha = 1, color = "black") +
	stat_compare_means(method = "wilcox", paired = FALSE, label = "p", comparisons = list(c("Model B", "Model F"), c("Model M", "Model F")), tip.length = 0, bracket.size = 0.7, vjust = -0.4, step.increase = 0.2) + 
	facet_wrap(~alg) +
	xlab("")+
	scale_y_continuous(expand = expansion(mult = c(0.05,.3)))+
	theme_linedraw(base_size = 15)+
	theme(panel.grid = element_blank(), legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))

p2 <- ggplot(male_met %>% mutate(model = factor(model, levels = c("Model B", "Model M", "Model F"))), aes(x = model, y = MAE)) + 
	geom_boxplot(aes(color = model), width = 0.6) +
	scale_color_manual(values = c("Model B" = "grey", "Model F" = "#E15566", "Model M" = "#4981BF")) +
	geom_point(position = position_jitter(width = 0.2), size = 0.5, alpha = 1, color = "black") +
	stat_compare_means(method = "wilcox", paired = FALSE, label = "p", comparisons = list(c("Model B", "Model M"), c("Model M", "Model F")), tip.length = 0, bracket.size = 0.7, vjust = -0.4, step.increase = 0.2) + 
	facet_wrap(~alg) +
	xlab("")+
	scale_y_continuous(expand = expansion(mult = c(0.05,.3)))+
	theme_linedraw(base_size = 15)+
	theme(panel.grid = element_blank(), legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))

ggarrange(p1, p2, ncol = 2, nrow = 1)

p1 <- ggplot(female_met, aes(x = model, y = r)) + 
	geom_boxplot(aes(color = model), width = 0.6) +
	scale_color_manual(values = c("Model B" = "grey", "Model F" = "#E15566", "Model M" = "#4981BF")) +
	geom_point(position = position_jitter(width = 0.2), size = 0.5, alpha = 1, color = "black") +
	stat_compare_means(method = "wilcox", paired = FALSE, label = "p", comparisons = list(c("Model B", "Model F"), c("Model M", "Model F")), tip.length = 0, bracket.size = 0.7, vjust = -0.4, step.increase = 0.2) + 
	facet_wrap(~alg) +
	xlab("")+
	scale_y_continuous(expand = expansion(mult = c(0.05,.3)))+
	theme_linedraw(base_size = 15)+
	theme(panel.grid = element_blank(), legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))

p2 <- ggplot(male_met %>% mutate(model = factor(model, levels = c("Model B", "Model M", "Model F"))), aes(x = model, y = r)) + 
	geom_boxplot(aes(color = model), width = 0.6) +
	scale_color_manual(values = c("Model B" = "grey", "Model F" = "#E15566", "Model M" = "#4981BF")) +
	geom_point(position = position_jitter(width = 0.2), size = 0.5, alpha = 1, color = "black") +
	stat_compare_means(method = "wilcox", paired = FALSE, label = "p", comparisons = list(c("Model B", "Model M"), c("Model M", "Model F")), tip.length = 0, bracket.size = 0.7, vjust = -0.4, step.increase = 0.2) + 
	facet_wrap(~alg) +
	xlab("")+
	scale_y_continuous(expand = expansion(mult = c(0.05,.3)))+
	theme_linedraw(base_size = 15)+
	theme(panel.grid = element_blank(), legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))

ggarrange(p1, p2, ncol = 2, nrow = 1)

lapply(list("Model B", "Model F", "Model M"), function(m){
p1 <- ggplot(female_met %>% filter(model == m), aes(x = alg, y = RMSE)) + 
	geom_boxplot(aes(color = alg), width = 0.6) +
	scale_color_manual(values = c("Elastic Net" = "#4DBBD5", "XGBoost" = "#E64B35", "MLP" = "#00A087")) +
	geom_point(position = position_jitter(width = 0.2), size = 0.5, alpha = 1, color = "black") +
	stat_compare_means(method = "wilcox", paired = FALSE, label = "p", comparisons = list(c("Elastic Net", "XGBoost"), c("XGBoost", "MLP"), c("Elastic Net", "MLP")), tip.length = 0, bracket.size = 0.7, vjust = -0.4, step.increase = 0.3) + 
	facet_wrap(~model) +
	xlab("")+
	scale_y_continuous(expand = expansion(mult = c(0.05,.25)))+
	theme_linedraw(base_size = 15)+
	theme(panel.grid = element_blank(), legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))

p2 <- ggplot(female_met %>% filter(model == m), aes(x = alg, y = r)) + 
	geom_boxplot(aes(color = alg), width = 0.6) +
	scale_color_manual(values = c("Elastic Net" = "#4DBBD5", "XGBoost" = "#E64B35", "MLP" = "#00A087")) +
	geom_point(position = position_jitter(width = 0.2), size = 0.5, alpha = 1, color = "black") +
	stat_compare_means(method = "wilcox", paired = FALSE, label = "p", comparisons = list(c("Elastic Net", "XGBoost"), c("XGBoost", "MLP"), c("Elastic Net", "MLP")), tip.length = 0, bracket.size = 0.7, vjust = -0.4, step.increase = 0.3) + 
	facet_wrap(~model) +
	xlab("")+
	scale_y_continuous(expand = expansion(mult = c(0.05,.25)))+
	theme_linedraw(base_size = 15)+
	theme(panel.grid = element_blank(), legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))

ggarrange(p1, p2, ncol = 4, nrow = 1)

})

p1 <- ggplot(female_met %>% filter(alg == "MLP"), aes(x = model, y = RMSE)) + 
	geom_boxplot(aes(color = model), width = 0.6) +
	scale_color_manual(values = c("Model B" = "grey", "Model F" = "#E15566", "Model M" = "#4981BF")) +
	geom_point(position = position_jitter(width = 0.2), size = 0.5, alpha = 1, color = "black") +
	stat_compare_means(method = "wilcox", paired = FALSE, label = "p", comparisons = list(c("Model B", "Model F"), c("Model M", "Model F")), tip.length = 0, bracket.size = 0.7, vjust = -0.4, step.increase = 0.2) + 
	facet_wrap(~type) +
	xlab("")+
	scale_y_continuous(expand = expansion(mult = c(0.05,.3)))+
	theme_linedraw(base_size = 15)+
	theme(panel.grid = element_blank(), legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))

p2 <- ggplot(male_met %>% filter(alg == "MLP") %>% mutate(model = factor(model, levels = c("Model B", "Model M", "Model F"))), aes(x = model, y = RMSE)) + 
	geom_boxplot(aes(color = model), width = 0.6) +
	scale_color_manual(values = c("Model B" = "grey", "Model F" = "#E15566", "Model M" = "#4981BF")) +
	geom_point(position = position_jitter(width = 0.2), size = 0.5, alpha = 1, color = "black") +
	stat_compare_means(method = "wilcox", paired = FALSE, label = "p", comparisons = list(c("Model B", "Model M"), c("Model M", "Model F")), tip.length = 0, bracket.size = 0.7, vjust = -0.4, step.increase = 0.2) + 
	facet_wrap(~type) +
	xlab("")+
	scale_y_continuous(expand = expansion(mult = c(0.05,.3)))+
	theme_linedraw(base_size = 15)+
	theme(panel.grid = element_blank(), legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))

p3 <- ggplot(female_met %>% filter(alg == "MLP"), aes(x = model, y = MAE)) + 
	geom_boxplot(aes(color = model), width = 0.6) +
	scale_color_manual(values = c("Model B" = "grey", "Model F" = "#E15566", "Model M" = "#4981BF")) +
	geom_point(position = position_jitter(width = 0.2), size = 0.5, alpha = 1, color = "black") +
	stat_compare_means(method = "wilcox", paired = FALSE, label = "p", comparisons = list(c("Model B", "Model F"), c("Model M", "Model F")), tip.length = 0, bracket.size = 0.7, vjust = -0.4, step.increase = 0.2) + 
	facet_wrap(~type) +
	xlab("")+
	scale_y_continuous(expand = expansion(mult = c(0.05,.3)))+
	theme_linedraw(base_size = 15)+
	theme(panel.grid = element_blank(), legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))

p4 <- ggplot(male_met %>% filter(alg == "MLP") %>% mutate(model = factor(model, levels = c("Model B", "Model M", "Model F"))), aes(x = model, y = MAE)) + 
	geom_boxplot(aes(color = model), width = 0.6) +
	scale_color_manual(values = c("Model B" = "grey", "Model F" = "#E15566", "Model M" = "#4981BF")) +
	geom_point(position = position_jitter(width = 0.2), size = 0.5, alpha = 1, color = "black") +
	stat_compare_means(method = "wilcox", paired = FALSE, label = "p", comparisons = list(c("Model B", "Model M"), c("Model M", "Model F")), tip.length = 0, bracket.size = 0.7, vjust = -0.4, step.increase = 0.2) + 
	facet_wrap(~type) +
	xlab("")+
	scale_y_continuous(expand = expansion(mult = c(0.05,.3)))+
	theme_linedraw(base_size = 15)+
	theme(panel.grid = element_blank(), legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))

p5 <- ggplot(female_met %>% filter(alg == "MLP"), aes(x = model, y = r)) + 
	geom_boxplot(aes(color = model), width = 0.6) +
	scale_color_manual(values = c("Model B" = "grey", "Model F" = "#E15566", "Model M" = "#4981BF")) +
	geom_point(position = position_jitter(width = 0.2), size = 0.5, alpha = 1, color = "black") +
	stat_compare_means(method = "wilcox", paired = FALSE, label = "p", comparisons = list(c("Model B", "Model F"), c("Model M", "Model F")), tip.length = 0, bracket.size = 0.7, vjust = -0.4, step.increase = 0.2) + 
	facet_wrap(~type) +
	xlab("")+
	scale_y_continuous(expand = expansion(mult = c(0.05,.3)))+
	theme_linedraw(base_size = 15)+
	theme(panel.grid = element_blank(), legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))

p6 <- ggplot(male_met %>% filter(alg == "MLP") %>% mutate(model = factor(model, levels = c("Model B", "Model M", "Model F"))), aes(x = model, y = r)) + 
	geom_boxplot(aes(color = model), width = 0.6) +
	scale_color_manual(values = c("Model B" = "grey", "Model F" = "#E15566", "Model M" = "#4981BF")) +
	geom_point(position = position_jitter(width = 0.2), size = 0.5, alpha = 1, color = "black") +
	stat_compare_means(method = "wilcox", paired = FALSE, label = "p", comparisons = list(c("Model B", "Model M"), c("Model M", "Model F")), tip.length = 0, bracket.size = 0.7, vjust = -0.4, step.increase = 0.2) + 
	facet_wrap(~type) +
	xlab("")+
	scale_y_continuous(expand = expansion(mult = c(0.05,.3)))+
	theme_linedraw(base_size = 15)+
	theme(panel.grid = element_blank(), legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))

ggarrange(p1, p2, p3, p4, p5, p6, ncol = 4, nrow = 1)
dev.off()
