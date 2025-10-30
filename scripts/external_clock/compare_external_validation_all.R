library(rio)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(svglite)
library(Metrics)

## Function
PlotMetricsGender <- function(df){

female_df <- df %>% filter(sex == "female")
male_df <- df %>% filter(sex == "male")

female_met <- female_df %>% group_by(alg, model, fold) %>% 
	summarise(RMSE = rmse(actual_age, predicted_age), MAE = mae(actual_age, predicted_age), .groups = "drop")
male_met <- male_df %>% group_by(alg, model, fold) %>% 
	summarise(RMSE = rmse(actual_age, predicted_age), MAE = mae(actual_age, predicted_age), .groups = "drop")
print(head(female_met))
print(head(male_met))

p1 <- ggplot(female_met, aes(x = model, y = RMSE)) + 
	geom_boxplot(aes(color = model), width = 0.6) +
	scale_color_manual(values = c("Model B" = "black", "Model F" = "#E15566", "Model M" = "#4981BF")) +
	geom_point(position = position_jitter(width = 0.2), size = 1, alpha = 0.5, color = "#E15566") +
	stat_compare_means(method = "wilcox", paired = FALSE, label = "p", comparisons = list(c("Model B", "Model F"), c("Model M", "Model F")), tip.length = 0, bracket.size = 0.7, vjust = -0.4, step_increase = 2) + 
	facet_wrap(~alg) +
	xlab("")+
	scale_y_continuous(expand = expansion(mult = c(0,.35)))+
	theme_linedraw(base_size = 15)+
	theme(panel.grid = element_blank(), legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))

p2 <- ggplot(male_met, aes(x = model, y = RMSE)) + 
	geom_boxplot(aes(color = model), width = 0.6) +
	scale_color_manual(values = c("Model B" = "black", "Model M" = "#4981BF", "Model F" = "#E15566")) +
	geom_point(position = position_jitter(width = 0.2), size = 1, alpha = 0.5, color = "#4981BF") +
	stat_compare_means(method = "wilcox", paired = FALSE, label = "p", comparisons = list(c("Model B", "Model M"), c("Model M", "Model F")), tip.length = 0, bracket.size = 0.7, vjust = -0.4, step_increase = 2) + 
	facet_wrap(~alg) +
	xlab("")+
	scale_y_continuous(expand = expansion(mult = c(0,.35)))+
	theme_linedraw(base_size = 15)+
	theme(panel.grid = element_blank(), legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))

p3 <- ggplot(female_met, aes(x = model, y = MAE)) + 
	geom_boxplot(aes(color = model), width = 0.6) +
	scale_color_manual(values = c("Model B" = "black", "Model F" = "#E15566", "Model M" = "#4981BF")) +
	geom_point(position = position_jitter(width = 0.2), size = 1, alpha = 0.5, color = "#E15566") +
	stat_compare_means(method = "wilcox", paired = FALSE, label = "p", comparisons = list(c("Model B", "Model F"), c("Model M", "Model F")), tip.length = 0, bracket.size = 0.7, vjust = -0.4, step_increase = 2) + 
	facet_wrap(~alg) +
	xlab("")+
	scale_y_continuous(expand = expansion(mult = c(0,.35)))+
	theme_linedraw(base_size = 15)+
	theme(panel.grid = element_blank(), legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))

p4 <- ggplot(male_met, aes(x = model, y = MAE)) + 
	geom_boxplot(aes(color = model), width = 0.6) +
	scale_color_manual(values = c("Model B" = "black", "Model M" = "#4981BF", "Model F" = "#E15566")) +
	geom_point(position = position_jitter(width = 0.2), size = 1, alpha = 0.5, color = "#4981BF") +
	stat_compare_means(method = "wilcox", paired = FALSE, label = "p", comparisons = list(c("Model B", "Model M"), c("Model M", "Model F")), tip.length = 0, bracket.size = 0.7, vjust = -0.4, step_increase = 2) + 
	facet_wrap(~alg) +
	xlab("")+
	scale_y_continuous(expand = expansion(mult = c(0,.35)))+
	theme_linedraw(base_size = 15)+
	theme(panel.grid = element_blank(), legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))

ggarrange(p1, p2, p3, p4, ncol=2, nrow=2)

}

PlotScatterGender <- function(df){

lapply(list("Elastic Net", "XGBoost", "MLP"), function(a){

df_alg <- df %>% filter(alg == a)

plots <- lapply(as.list(1:5), function(f){

	p1 <- ggscatter(df_alg %>% filter(fold == f, model == "Model B"), x = "actual_age", y = "predicted_age", color = "grey", cor.coef = TRUE, alpha = 0.7, cor.coeff.args = list(method = "pearson"), cor.coef.size = 2.5) +
		geom_smooth(method = "lm", color = "black", se = F) + theme_test() + theme(legend.position = "none") + xlim(c(0,100)) + ylim(c(0,100))
	p2 <- ggscatter(df_alg %>% filter(fold == f, model == "Model F", sex == "female"), x = "actual_age", y = "predicted_age", color = "#E15566", cor.coef = TRUE, alpha = 0.7, cor.coeff.args = list(method = "pearson"), cor.coef.size = 2.5) +
		geom_smooth(method = "lm", color = "black", se = F) + theme_test() + theme(legend.position = "none") + xlim(c(0,100)) + ylim(c(0,100))
	p3 <- ggscatter(df_alg %>% filter(fold == f, model == "Model M", sex == "male"), x = "actual_age", y = "predicted_age", color = "#4981BF", cor.coef = TRUE, alpha = 0.7, cor.coeff.args = list(method = "pearson"), cor.coef.size = 2.5) +
		geom_smooth(method = "lm", color = "black", se = F) + theme_test() + theme(legend.position = "none") + xlim(c(0,100)) + ylim(c(0,100))
	return(ggarrange(p1, p2, p3, ncol = 1, nrow = 3))

})

plot(ggarrange(plotlist = plots, ncol = 5, nrow = 1))
})

}

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

preds_age <- preds_age %>% filter(sample_id %in% preds_age_table[!duplicated(preds_age_table$donor_id), "sample_id"])
print(preds_age %>% distinct(sample_id, donor_id, dataset, disease) %>% count(dataset, disease) %>% group_by(disease) %>% mutate(prop = prop.table(n)))

preds_age <- preds_age %>% mutate(agediff = predicted_age - actual_age)

pdf(snakemake@output[["plot1"]], width = 9, height = 6)

PlotMetricsGender(preds_age)
PlotMetricsGender(preds_age %>% filter(disease == "normal"))
PlotMetricsGender(preds_age %>% filter(disease != "normal"))

dev.off()

pdf(snakemake@output[["plot2"]], width = 11, height = 6)

PlotScatterGender(preds_age)
PlotScatterGender(preds_age %>% filter(disease == "normal"))
PlotScatterGender(preds_age %>% filter(disease != "normal"))

dev.off()

pdf(snakemake@output[["plot3"]], width = 9, height = 9)
lapply(list("Elastic Net", "XGBoost", "MLP"), function(a){

ch_p <- lapply(list("Model B", "Model F", "Model M"), function(m){

ch_df <- preds_age %>% filter(alg == a, dataset == "ch", model == m) %>% mutate(disease = factor(disease, levels = c("normal", "clonal hematopoiesis")))
p <- ggplot(df, aes(x = disease, y = agediff)) +
	geom_boxplot(color = "black")+
	geom_jitter(aes(color = disease), position=position_jitter(0.2), alpha = 0.5)+
	stat_compare_means(method = "wilcox", paired = FALSE, label = "p", comparisons = list(c("normal", "clonal hematopoiesis")), tip.length = 0, bracket.size = 0.7, vjust = -0.4, step_increase = 2) + 
	facet_wrap(~fold, ncol = 5) +
	xlab("")+
	scale_y_continuous(expand = expansion(mult = c(0,.35)))+
	theme_linedraw(base_size = 15)+
	theme(panel.grid = element_blank(), legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))

return(p)
})
plot(ggarrange(plotlist = ch_p, ncol = 1, nrow = 3))
	
glaucoma_p <- lapply(list("Model B", "Model F", "Model M"), function(m){

glaucoma_df <- preds_age %>% filter(alg == a, dataset == "glaucoma", model == m) %>% mutate(disease = factor(disease, levels = c("normal", "open-angle glaucoma")))
p <- ggplot(df, aes(x = disease, y = agediff)) +
	geom_boxplot(color = "black")+
	geom_jitter(aes(color = disease), position=position_jitter(0.2), alpha = 0.5)+
	stat_compare_means(method = "wilcox", paired = FALSE, label = "p", comparisons = list(c("normal", "open-angle glaucoma")), tip.length = 0, bracket.size = 0.7, vjust = -0.4, step_increase = 2) + 
	facet_wrap(~fold, ncol = 5) +
	xlab("")+
	scale_y_continuous(expand = expansion(mult = c(0,.35)))+
	theme_linedraw(base_size = 15)+
	theme(panel.grid = element_blank(), legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))

return(p)
})
plot(ggarrange(plotlist = glaucoma_p, ncol = 1, nrow = 3))
	
ren_p <- lapply(list("Model B", "Model F", "Model M"), function(m){

ren_df <- preds_age %>% filter(alg == a, dataset == "ren", model == m) %>% mutate(disease = factor(disease, levels = c("normal", "mild_progression", "mild_convalescence", "severe_progression", "severe_convalescence")))
p <- ggplot(df, aes(x = disease, y = agediff)) +
	geom_boxplot(color = "black")+
	geom_jitter(aes(color = disease), position=position_jitter(0.2), alpha = 0.5)+
	stat_compare_means(method = "wilcox", paired = FALSE, label = "p", comparisons = list(c("normal", "mild_progression"), c("mild_progression", "mild_convalescence"), c("severe_progression", "severe_convalescence"), c("normal", "severe_progression")), tip.length = 0, bracket.size = 0.7, vjust = -0.4, step_increase = 2) + 
	facet_wrap(~fold, ncol = 5) +
	xlab("")+
	scale_y_continuous(expand = expansion(mult = c(0,1)))+
	theme_linedraw(base_size = 15)+
	theme(panel.grid = element_blank(), legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))

return(p)
})
plot(ggarrange(plotlist = ren_p, ncol = 1, nrow = 3))
	
wellcome_p <- lapply(list("Model B", "Model F", "Model M"), function(m){

wellcome_df <- preds_age %>% filter(alg == a, dataset == "wellcome", model == m) %>% mutate(disease = factor(disease, levels = c("normal", "mild", "moderate", "severe", "critical")))
p <- ggplot(df, aes(x = disease, y = agediff)) +
	geom_boxplot(color = "black")+
	geom_jitter(aes(color = disease), position=position_jitter(0.2), alpha = 0.5)+
	stat_compare_means(method = "wilcox", paired = FALSE, label = "p", comparisons = list(c("normal", "mild"), c("mild", "moderate"), c("moderate", "severe"), c("severe", "critical"), c("normal", "moderate"), c("normal", "severe"), c("normal", "critical")), tip.length = 0, bracket.size = 0.7, vjust = -0.4, step_increase = 2) + 
	facet_wrap(~fold, ncol = 5) +
	xlab("")+
	scale_y_continuous(expand = expansion(mult = c(0,2)))+
	theme_linedraw(base_size = 15)+
	theme(panel.grid = element_blank(), legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))

return(p)
})
plot(ggarrange(plotlist = wellcome_p, ncol = 1, nrow = 3))
	
})
dev.off()
