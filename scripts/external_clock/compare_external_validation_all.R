library(rio)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(svglite)
library(Metrics)
library(svglite)

## Function
PlotMetricsCohort <- function(df){

both_h_df <- df %>% filter(disease == "normal", model == "Model B") %>% 
	group_by(fold) %>% 
	summarise(RMSE = rmse(actual_age, predicted_age), MAE = mae(actual_age, predicted_age), r = cor(actual_age,predicted_age), .groups = "drop") %>%
	mutate(type = "Both healthy", sex = "Both")
both_d_df <- df %>% filter(disease != "normal", model == "Model B") %>% 
	group_by(fold) %>% 
	summarise(RMSE = rmse(actual_age, predicted_age), MAE = mae(actual_age, predicted_age), r = cor(actual_age,predicted_age), .groups = "drop") %>%
	mutate(type = "Both disease", sex = "Both")
female_h_df <- df %>% filter(sex == "female", disease == "normal", model == "Model F") %>% 
	group_by(fold) %>% 
	summarise(RMSE = rmse(actual_age, predicted_age), MAE = mae(actual_age, predicted_age), r = cor(actual_age,predicted_age), .groups = "drop") %>%
	mutate(type = "Female healthy", sex = "Female")
female_d_df <- df %>% filter(sex == "female", disease != "normal", model == "Model F") %>% 
	group_by(fold) %>% 
	summarise(RMSE = rmse(actual_age, predicted_age), MAE = mae(actual_age, predicted_age), r = cor(actual_age,predicted_age), .groups = "drop") %>%
	mutate(type = "Female disease", sex = "Female")
male_h_df <- df %>% filter(sex == "male", disease == "normal", model == "Model M") %>% 
	group_by(fold) %>% 
	summarise(RMSE = rmse(actual_age, predicted_age), MAE = mae(actual_age, predicted_age), r = cor(actual_age,predicted_age), .groups = "drop") %>%
	mutate(type = "Male healthy", sex = "Male")
male_d_df <- df %>% filter(sex == "male", disease != "normal", model == "Model M") %>% 
	group_by(fold) %>% 
	summarise(RMSE = rmse(actual_age, predicted_age), MAE = mae(actual_age, predicted_age), r = cor(actual_age,predicted_age), .groups = "drop") %>%
	mutate(type = "Male disease", sex = "Male")

all <- bind_rows(both_h_df, both_d_df, female_h_df, female_d_df, male_h_df, male_d_df) %>% 
	mutate(sex = factor(sex, levels = c("Both", "Female", "Male")), 
		type = factor(type, levels = c("Both healthy", "Both disease", "Female healthy", "Female disease", "Male healthy", "Male disease")))

p1 <- ggplot(all, aes(x = type, y = RMSE)) + 
	geom_boxplot(aes(color = sex), width = 0.6) +
	scale_color_manual(values = c("Both" = "grey", "Female" = "#E15566", "Male" = "#4981BF")) +
	geom_point(position = position_jitter(width = 0.2), size = 0.5, alpha = 1, color = "black") +
	facet_wrap(~sex, scale = "free") +
	xlab("")+
	scale_y_continuous(expand = expansion(mult = c(0.05,0.3)))+
	theme_linedraw(base_size = 15)+
	theme(panel.grid = element_blank(), legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))
#p1 <- p1 + stat_compare_means(method = "wilcox", paired = FALSE, label = "p", tip.length = 0, bracket.size = 0.7, vjust = -0.4)
p2 <- ggplot(all, aes(x = type, y = MAE)) + 
	geom_boxplot(aes(color = sex), width = 0.6) +
	scale_color_manual(values = c("Both" = "grey", "Female" = "#E15566", "Male" = "#4981BF")) +
	geom_point(position = position_jitter(width = 0.2), size = 0.5, alpha = 1, color = "black") +
	facet_wrap(~sex, scale = "free") +
	xlab("")+
	scale_y_continuous(expand = expansion(mult = c(0.05,0.3)))+
	theme_linedraw(base_size = 15)+
	theme(panel.grid = element_blank(), legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))
#p2 <- p2 + stat_compare_means(method = "wilcox", paired = FALSE, label = "p", tip.length = 0, bracket.size = 0.7, vjust = -0.4)
p3 <- ggplot(all, aes(x = type, y = r)) + 
	geom_boxplot(aes(color = sex), width = 0.6) +
	scale_color_manual(values = c("Both" = "grey", "Female" = "#E15566", "Male" = "#4981BF")) +
	geom_point(position = position_jitter(width = 0.2), size = 0.5, alpha = 1, color = "black") +
	facet_wrap(~sex, scale = "free") +
	xlab("")+
	scale_y_continuous(expand = expansion(mult = c(0.05,0.3)))+
	theme_linedraw(base_size = 15)+
	theme(panel.grid = element_blank(), legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))
#p3 <- p3 + stat_compare_means(method = "wilcox", paired = FALSE, label = "p", tip.length = 0, bracket.size = 0.7, vjust = -0.4)

pall <- ggarrange(p1, p2, p3, ncol=3, nrow=1)

return(pall)

}

PlotMetricsGender <- function(df){

female_df <- df %>% filter(sex == "female")
male_df <- df %>% filter(sex == "male")

female_met <- female_df %>% group_by(alg, model, fold) %>% 
	summarise(RMSE = rmse(actual_age, predicted_age), MAE = mae(actual_age, predicted_age), r = cor(actual_age,predicted_age), .groups = "drop") %>%
	mutate(alg = factor(alg, levels = c("Elastic Net", "XGBoost", "MLP")), type = "female donors")
male_met <- male_df %>% group_by(alg, model, fold) %>% 
	summarise(RMSE = rmse(actual_age, predicted_age), MAE = mae(actual_age, predicted_age), r = cor(actual_age,predicted_age), .groups = "drop") %>%
	mutate(alg = factor(alg, levels = c("Elastic Net", "XGBoost", "MLP")), type = "male donors")
print(head(female_met))
print(head(male_met))

## metrics comparison of models (5 folds) in either female donors / male donors
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
	scale_color_manual(values = c("Model B" = "grey", "Model M" = "#4981BF", "Model F" = "#E15566")) +
	geom_point(position = position_jitter(width = 0.2), size = 0.5, alpha = 1, color = "black") +
	stat_compare_means(method = "wilcox", paired = FALSE, label = "p", comparisons = list(c("Model B", "Model M"), c("Model M", "Model F")), tip.length = 0, bracket.size = 0.7, vjust = -0.4, step.increase = 0.2) + 
	facet_wrap(~alg) +
	xlab("")+
	scale_y_continuous(expand = expansion(mult = c(0.05,.3)))+
	theme_linedraw(base_size = 15)+
	theme(panel.grid = element_blank(), legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))

p3 <- ggplot(female_met, aes(x = model, y = MAE)) + 
	geom_boxplot(aes(color = model), width = 0.6) +
	scale_color_manual(values = c("Model B" = "grey", "Model F" = "#E15566", "Model M" = "#4981BF")) +
	geom_point(position = position_jitter(width = 0.2), size = 0.5, alpha = 1, color = "black") +
	stat_compare_means(method = "wilcox", paired = FALSE, label = "p", comparisons = list(c("Model B", "Model F"), c("Model M", "Model F")), tip.length = 0, bracket.size = 0.7, vjust = -0.4, step.increase = 0.2) + 
	facet_wrap(~alg) +
	xlab("")+
	scale_y_continuous(expand = expansion(mult = c(0.05,.3)))+
	theme_linedraw(base_size = 15)+
	theme(panel.grid = element_blank(), legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))

p4 <- ggplot(male_met %>% mutate(model = factor(model, levels = c("Model B", "Model M", "Model F"))), aes(x = model, y = MAE)) + 
	geom_boxplot(aes(color = model), width = 0.6) +
	scale_color_manual(values = c("Model B" = "grey", "Model M" = "#4981BF", "Model F" = "#E15566")) +
	geom_point(position = position_jitter(width = 0.2), size = 0.5, alpha = 1, color = "black") +
	stat_compare_means(method = "wilcox", paired = FALSE, label = "p", comparisons = list(c("Model B", "Model M"), c("Model M", "Model F")), tip.length = 0, bracket.size = 0.7, vjust = -0.4, step.increase = 0.2) + 
	facet_wrap(~alg) +
	xlab("")+
	scale_y_continuous(expand = expansion(mult = c(0.05,.3)))+
	theme_linedraw(base_size = 15)+
	theme(panel.grid = element_blank(), legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))

p5 <- ggplot(female_met, aes(x = model, y = r)) + 
	geom_boxplot(aes(color = model), width = 0.6) +
	scale_color_manual(values = c("Model B" = "grey", "Model F" = "#E15566", "Model M" = "#4981BF")) +
	geom_point(position = position_jitter(width = 0.2), size = 0.5, alpha = 1, color = "black") +
	stat_compare_means(method = "wilcox", paired = FALSE, label = "p", comparisons = list(c("Model B", "Model F"), c("Model M", "Model F")), tip.length = 0, bracket.size = 0.7, vjust = -0.4, step.increase = 0.2) + 
	facet_wrap(~alg) +
	xlab("")+
	scale_y_continuous(expand = expansion(mult = c(0.05,.3)))+
	theme_linedraw(base_size = 15)+
	theme(panel.grid = element_blank(), legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))

p6 <- ggplot(male_met %>% mutate(model = factor(model, levels = c("Model B", "Model M", "Model F"))), aes(x = model, y = r)) + 
	geom_boxplot(aes(color = model), width = 0.6) +
	scale_color_manual(values = c("Model B" = "grey", "Model M" = "#4981BF", "Model F" = "#E15566")) +
	geom_point(position = position_jitter(width = 0.2), size = 0.5, alpha = 1, color = "black") +
	stat_compare_means(method = "wilcox", paired = FALSE, label = "p", comparisons = list(c("Model B", "Model M"), c("Model M", "Model F")), tip.length = 0, bracket.size = 0.7, vjust = -0.4, step.increase = 0.2) + 
	facet_wrap(~alg) +
	xlab("")+
	scale_y_continuous(expand = expansion(mult = c(0.05,.3)))+
	theme_linedraw(base_size = 15)+
	theme(panel.grid = element_blank(), legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))

pall1 <- ggarrange(p1, p2, p3, p4, p5, p6, ncol=2, nrow=3)

## Algorithm comparison
models <- lapply(list("Model B", "Model F", "Model M"), function(m){
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

p2 <- ggplot(male_met %>% filter(model == m), aes(x = alg, y = r)) +
        geom_boxplot(aes(color = alg), width = 0.6) +
        scale_color_manual(values = c("Elastic Net" = "#4DBBD5", "XGBoost" = "#E64B35", "MLP" = "#00A087")) +
        geom_point(position = position_jitter(width = 0.2), size = 0.5, alpha = 1, color = "black") +
        stat_compare_means(method = "wilcox", paired = FALSE, label = "p", comparisons = list(c("Elastic Net", "XGBoost"), c("XGBoost", "MLP"), c("Elastic Net", "MLP")), tip.length = 0, bracket.size = 0.7, vjust = -0.4, step.increase = 0.3) +
        facet_wrap(~model) +
        xlab("")+
        scale_y_continuous(expand = expansion(mult = c(0.05,.25)))+
        theme_linedraw(base_size = 15)+
        theme(panel.grid = element_blank(), legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))

return(ggarrange(p1, p2, ncol = 4, nrow = 1))

})

pall2 <- ggarrange(plotlist = models, ncol = 1, nrow = 3)

## Gender-combined or Gender-specific model comparison for MLP models in both female and male donors
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
	scale_color_manual(values = c("Model B" = "grey", "Model M" = "#4981BF", "Model F" = "#E15566")) +
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
	scale_color_manual(values = c("Model B" = "grey", "Model M" = "#4981BF", "Model F" = "#E15566")) +
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
	scale_color_manual(values = c("Model B" = "grey", "Model M" = "#4981BF", "Model F" = "#E15566")) +
	geom_point(position = position_jitter(width = 0.2), size = 0.5, alpha = 1, color = "black") +
	stat_compare_means(method = "wilcox", paired = FALSE, label = "p", comparisons = list(c("Model B", "Model M"), c("Model M", "Model F")), tip.length = 0, bracket.size = 0.7, vjust = -0.4, step.increase = 0.2) + 
	facet_wrap(~type) +
	xlab("")+
	scale_y_continuous(expand = expansion(mult = c(0.05,.3)))+
	theme_linedraw(base_size = 15)+
	theme(panel.grid = element_blank(), legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))

pall3 <- ggarrange(p1, p2, p3, p4, p5, p6, ncol=4, nrow=3)


plot(pall1)
plot(pall2)
plot(pall3)

}

PlotScatterGender <- function(df){

lapply(list("Elastic Net", "XGBoost", "MLP"), function(a){

df_alg <- df %>% filter(alg == a)

plots <- lapply(as.list(1:5), function(f){
	max_predicted_age <- max(df_alg %>% pull(predicted_age))

	p1 <- ggscatter(df_alg %>% filter(fold == f, model == "Model B"), x = "actual_age", y = "predicted_age", color = "grey", cor.coef = TRUE, alpha = 0.7, cor.coeff.args = list(method = "pearson"), cor.coef.size = 4) +
		geom_smooth(method = "lm", color = "black", se = F) + theme_test(base_size = 15) + theme(legend.position = "none") + xlim(c(0,max(100, max_predicted_age))) + ylim(c(0,max(100, max_predicted_age)))
	p2 <- ggscatter(df_alg %>% filter(fold == f, model == "Model F", sex == "female"), x = "actual_age", y = "predicted_age", color = "#E15566", cor.coef = TRUE, alpha = 0.7, cor.coeff.args = list(method = "pearson"), cor.coef.size = 4) +
		geom_smooth(method = "lm", color = "black", se = F) + theme_test(base_size = 15) + theme(legend.position = "none") + xlim(c(0,max(100, max_predicted_age))) + ylim(c(0,max(100, max_predicted_age)))
	p3 <- ggscatter(df_alg %>% filter(fold == f, model == "Model M", sex == "male"), x = "actual_age", y = "predicted_age", color = "#4981BF", cor.coef = TRUE, alpha = 0.7, cor.coeff.args = list(method = "pearson"), cor.coef.size = 4) +
		geom_smooth(method = "lm", color = "black", se = F) + theme_test(base_size = 15) + theme(legend.position = "none") + xlim(c(0,max(100, max_predicted_age))) + ylim(c(0,max(100, max_predicted_age)))
	return(ggarrange(p1, p2, p3, ncol = 1, nrow = 3))

})

plot(ggarrange(plotlist = plots, ncol = 5, nrow = 1))
})

}

## Main
pred_b1 <- import(snakemake@input[["pred_b1"]]) %>% mutate(color_name = "grey", model = "Model B", alg = "MLP")
pred_f1 <- import(snakemake@input[["pred_f1"]]) %>% mutate(color_name = "#E15566", model = "Model F", alg = "MLP")
pred_m1 <- import(snakemake@input[["pred_m1"]]) %>% mutate(color_name = "#4981BF", model = "Model M", alg = "MLP")
pred_b2 <- import(snakemake@input[["pred_b2"]]) %>% mutate(color_name = "grey", model = "Model B", alg = "Elastic Net")
pred_f2 <- import(snakemake@input[["pred_f2"]]) %>% mutate(color_name = "#E15566", model = "Model F", alg = "Elastic Net")
pred_m2 <- import(snakemake@input[["pred_m2"]]) %>% mutate(color_name = "#4981BF", model = "Model M", alg = "Elastic Net")
pred_b3 <- import(snakemake@input[["pred_b3"]]) %>% mutate(color_name = "grey", model = "Model B", alg = "XGBoost")
pred_f3 <- import(snakemake@input[["pred_f3"]]) %>% mutate(color_name = "#E15566", model = "Model F", alg = "XGBoost")
pred_m3 <- import(snakemake@input[["pred_m3"]]) %>% mutate(color_name = "#4981BF", model = "Model M", alg = "XGBoost")

preds <- bind_rows(pred_m1, pred_m2, pred_m3, pred_b1, pred_b2, pred_b3, pred_f1, pred_f2, pred_f3) %>% mutate(alg = factor(alg, c("Elastic Net", "XGBoost", "MLP")))
print(preds %>% distinct(sample_id, dataset, disease) %>% count(dataset, disease) %>% group_by(disease) %>% mutate(prop = prop.table(n)))

# External data with actual age
preds_age <- preds %>% mutate(actual_age = as.integer(actual_age)) %>% filter(!is.na(actual_age))
preds_age_table <- preds_age %>% distinct(sample_id, donor_id, dataset, disease) 
print(preds_age_table %>% count(dataset, disease) %>% group_by(disease) %>% mutate(prop = prop.table(n)))
print(preds_age_table[!duplicated(preds_age_table$donor_id),])

## This was for individual-level validation
#preds_age <- preds_age %>% filter(sample_id %in% preds_age_table[!duplicated(preds_age_table$donor_id), "sample_id"])
#print(preds_age %>% distinct(sample_id, donor_id, dataset, disease) %>% count(dataset, disease) %>% group_by(disease) %>% mutate(prop = prop.table(n)))

preds_age <- preds_age %>% mutate(agediff = predicted_age - actual_age)

## tables with all age available data (both healthy and disease)
preds <- preds_age
female_preds <- preds %>% filter(sex == "female")
male_preds <- preds %>% filter(sex == "male")

both_met <- preds %>% filter(model == "Model B") %>% group_by(alg, model, fold) %>% 
	summarise(RMSE = rmse(actual_age, predicted_age), MAE = mae(actual_age, predicted_age), r = cor(actual_age, predicted_age), .groups = "drop") %>% 
	mutate(alg = factor(alg, levels = c("Elastic Net", "XGBoost", "MLP")), type = "all donors")
female_met <- female_preds %>% group_by(alg, model, fold) %>% 
	summarise(RMSE = rmse(actual_age, predicted_age), MAE = mae(actual_age, predicted_age), r = cor(actual_age,predicted_age), .groups = "drop") %>%
	mutate(alg = factor(alg, levels = c("Elastic Net", "XGBoost", "MLP")), type = "female donors")
male_met <- male_preds %>% group_by(alg, model, fold) %>% 
	summarise(RMSE = rmse(actual_age, predicted_age), MAE = mae(actual_age, predicted_age), r = cor(actual_age,predicted_age), .groups = "drop") %>%
	mutate(alg = factor(alg, levels = c("Elastic Net", "XGBoost", "MLP")), type = "male donors")
print(head(female_met))
print(head(male_met))

met_all1 <- bind_rows(both_met, female_met, male_met) %>% mutate(type2 = "all")

## tables with only healthy
preds <- preds_age %>% filter(disease == "normal")
female_preds <- preds %>% filter(sex == "female")
male_preds <- preds %>% filter(sex == "male")

both_met <- preds %>% filter(model == "Model B") %>% group_by(alg, model, fold) %>% 
	summarise(RMSE = rmse(actual_age, predicted_age), MAE = mae(actual_age, predicted_age), r = cor(actual_age, predicted_age), .groups = "drop") %>% 
	mutate(alg = factor(alg, levels = c("Elastic Net", "XGBoost", "MLP")), type = "all donors")
female_met <- female_preds %>% group_by(alg, model, fold) %>% 
	summarise(RMSE = rmse(actual_age, predicted_age), MAE = mae(actual_age, predicted_age), r = cor(actual_age,predicted_age), .groups = "drop") %>%
	mutate(alg = factor(alg, levels = c("Elastic Net", "XGBoost", "MLP")), type = "female donors")
male_met <- male_preds %>% group_by(alg, model, fold) %>% 
	summarise(RMSE = rmse(actual_age, predicted_age), MAE = mae(actual_age, predicted_age), r = cor(actual_age,predicted_age), .groups = "drop") %>%
	mutate(alg = factor(alg, levels = c("Elastic Net", "XGBoost", "MLP")), type = "male donors")
print(head(female_met))
print(head(male_met))

met_all2 <- bind_rows(both_met, female_met, male_met) %>% mutate(type2 = "healthy")

## tables with only disease
preds <- preds_age %>% filter(disease != "normal")
female_preds <- preds %>% filter(sex == "female")
male_preds <- preds %>% filter(sex == "male")

both_met <- preds %>% filter(model == "Model B") %>% group_by(alg, model, fold) %>% 
	summarise(RMSE = rmse(actual_age, predicted_age), MAE = mae(actual_age, predicted_age), r = cor(actual_age, predicted_age), .groups = "drop") %>% 
	mutate(alg = factor(alg, levels = c("Elastic Net", "XGBoost", "MLP")), type = "all donors")
female_met <- female_preds %>% group_by(alg, model, fold) %>% 
	summarise(RMSE = rmse(actual_age, predicted_age), MAE = mae(actual_age, predicted_age), r = cor(actual_age,predicted_age), .groups = "drop") %>%
	mutate(alg = factor(alg, levels = c("Elastic Net", "XGBoost", "MLP")), type = "female donors")
male_met <- male_preds %>% group_by(alg, model, fold) %>% 
	summarise(RMSE = rmse(actual_age, predicted_age), MAE = mae(actual_age, predicted_age), r = cor(actual_age,predicted_age), .groups = "drop") %>%
	mutate(alg = factor(alg, levels = c("Elastic Net", "XGBoost", "MLP")), type = "male donors")
print(head(female_met))
print(head(male_met))

met_all3 <- bind_rows(both_met, female_met, male_met) %>% mutate(type2 = "disease")
met_all <- bind_rows(met_all1, met_all2, met_all3) %>% select(alg, model, type, type2, fold, RMSE, MAE, r)
export(met_all, snakemake@output[["fold"]])

met_all_sum <- met_all %>% group_by(alg, model, type) %>%
        summarise(median_RMSE = median(RMSE), median_MAE = median(MAE), median_r = median(r), .group = "drop")
head(met_all_sum)
export(met_all_sum, snakemake@output[["sum"]])

pdf(snakemake@output[["plot1"]], width = 9, height = 10.5)

PlotMetricsGender(preds_age)
PlotMetricsGender(preds_age %>% filter(disease == "normal"))
PlotMetricsGender(preds_age %>% filter(disease != "normal"))

dev.off()

pdf(snakemake@output[["plot2"]], width = 15, height = 8.5)

PlotScatterGender(preds_age)
PlotScatterGender(preds_age %>% filter(disease == "normal"))
PlotScatterGender(preds_age %>% filter(disease != "normal"))

dev.off()

## By disease
pdf(snakemake@output[["plot3"]], width = 9, height = 15)
lapply(list("Elastic Net", "XGBoost", "MLP"), function(a){

ch_p <- lapply(list("Model B", "Model F", "Model M"), function(m){

ch_df <- preds_age %>% filter(alg == a, dataset == "ch", model == m) %>% mutate(disease = factor(disease, levels = c("normal", "clonal hematopoiesis")))
if(m == "Model F"){
	ch_df <- ch_df %>% filter(sex == "female")
}else if(m == "Model M"){
	ch_df <- ch_df %>% filter(sex == "male")
}
p <- ggplot(ch_df, aes(x = disease, y = agediff)) +
	geom_boxplot(color = "black")+
	geom_jitter(aes(color = disease), position=position_jitter(0.2), alpha = 0.5)+
	stat_compare_means(method = "wilcox", paired = FALSE, label = "p", comparisons = list(c("normal", "clonal hematopoiesis")), tip.length = 0, bracket.size = 0.7, vjust = -0.4, step.increase = .1) + 
	facet_wrap(~fold, ncol = 5) +
	xlab("")+
	scale_y_continuous(expand = expansion(mult = c(0,.35)))+
	ggtitle(paste(a, m))+
	theme_linedraw(base_size = 15)+
	theme(panel.grid = element_blank(), legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))

return(p)
})
plot(ggarrange(plotlist = ch_p, ncol = 1, nrow = 3))
	
glaucoma_p <- lapply(list("Model B", "Model F", "Model M"), function(m){

glaucoma_df <- preds_age %>% filter(alg == a, dataset == "glaucoma", model == m) %>% mutate(disease = factor(disease, levels = c("normal", "open-angle glaucoma")))
if(m == "Model F"){
	glaucoma_df <- glaucoma_df %>% filter(sex == "female")
}else if(m == "Model M"){
	glaucoma_df <- glaucoma_df %>% filter(sex == "male")
}
p <- ggplot(glaucoma_df, aes(x = disease, y = agediff)) +
	geom_boxplot(color = "black")+
	geom_jitter(aes(color = disease), position=position_jitter(0.2), alpha = 0.5)+
	stat_compare_means(method = "wilcox", paired = FALSE, label = "p", comparisons = list(c("normal", "open-angle glaucoma")), tip.length = 0, bracket.size = 0.7, vjust = -0.4, step.increase = .1) + 
	facet_wrap(~fold, ncol = 5) +
	xlab("")+
	ggtitle(paste(a, m))+
	scale_y_continuous(expand = expansion(mult = c(0,.35)))+
	theme_linedraw(base_size = 15)+
	theme(panel.grid = element_blank(), legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))

return(p)
})
plot(ggarrange(plotlist = glaucoma_p, ncol = 1, nrow = 3))
	
ren_p <- lapply(list("Model B", "Model F", "Model M"), function(m){

ren_df <- preds_age %>% filter(alg == a, dataset == "ren", model == m) %>% mutate(disease = factor(disease, levels = c("normal", "mild_progression", "mild_convalescence", "severe_progression", "severe_convalescence")))
if(m == "Model F"){
	ren_df <- ren_df %>% filter(sex == "female")
}else if(m == "Model M"){
	ren_df <- ren_df %>% filter(sex == "male")
}
p <- ggplot(ren_df, aes(x = disease, y = agediff)) +
	geom_boxplot(color = "black")+
	geom_jitter(aes(color = disease), position=position_jitter(0.2), alpha = 0.5)+
	stat_compare_means(method = "wilcox", paired = FALSE, label = "p", comparisons = list(c("normal", "mild_progression"), c("mild_progression", "mild_convalescence"), c("severe_progression", "severe_convalescence"), c("normal", "severe_progression")), tip.length = 0, bracket.size = 0.7, vjust = -0.4, step.increase = .15) + 
	facet_wrap(~fold, ncol = 5) +
	xlab("")+
	ggtitle(paste(a, m))+
	scale_y_continuous(expand = expansion(mult = c(0,.8)))+
	theme_linedraw(base_size = 15)+
	theme(panel.grid = element_blank(), legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))

return(p)
})
plot(ggarrange(plotlist = ren_p, ncol = 1, nrow = 3))
	
wellcome_p <- lapply(list("Model B", "Model F", "Model M"), function(m){

wellcome_df <- preds_age %>% filter(alg == a, dataset == "wellcome", model == m) %>% mutate(disease = factor(disease, levels = c("normal", "mild", "moderate", "severe", "critical")))
if(m == "Model F"){
	wellcome_df <- wellcome_df %>% filter(sex == "female")
}else if(m == "Model M"){
	wellcome_df <- wellcome_df %>% filter(sex == "male")
}
p <- ggplot(wellcome_df, aes(x = disease, y = agediff)) +
	geom_boxplot(color = "black")+
	geom_jitter(aes(color = disease), position=position_jitter(0.2), alpha = 0.5)+
	stat_compare_means(method = "wilcox", paired = FALSE, label = "p", comparisons = list(c("normal", "mild"), c("mild", "moderate"), c("moderate", "severe"), c("severe", "critical"), c("normal", "moderate"), c("normal", "severe"), c("normal", "critical")), tip.length = 0, bracket.size = 0.7, vjust = -0.4, step.increase = .5) + 
	facet_wrap(~fold, ncol = 5) +
	xlab("")+
	ggtitle(paste(a, m))+
	scale_y_continuous(expand = expansion(mult = c(0,1)))+
	theme_linedraw(base_size = 15)+
	theme(panel.grid = element_blank(), legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))

return(p)
})
plot(ggarrange(plotlist = wellcome_p, ncol = 1, nrow = 3))
	
})
dev.off()


mlp_comp <- PlotMetricsCohort(preds_age %>% filter(alg == "MLP"))
ggsave(snakemake@output[["plot4"]], mlp_comp, width = 12, height =3) 
