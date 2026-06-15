library(rio)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(ggfortify)


## Function
plotPredScatter <- function(df, var, title){
	max <- max(c(df$actual_age, df[,var]))
	b <- ggscatter(df %>% filter(model == "Model B"), x = "actual_age", y = var, color = "grey", size = 1, alpha = 0.7, 
		cor.coef = TRUE, cor.coeff.args = list(method = "pearson"), cor.coef.size = 4) +
                geom_smooth(method = "lm", color = "black", se = F) + 
		theme_test(base_size = 15) + theme(legend.position = "none") + 
		ggtitle(title) +
		xlim(c(0,max)) + ylim(c(0,max))
	f <- ggscatter(df %>% filter(model == "Model F"), x = "actual_age", y = var, color = "#E15566", size = 1, alpha = 0.7, 
		cor.coef = TRUE, cor.coeff.args = list(method = "pearson"), cor.coef.size = 4) +
                geom_smooth(method = "lm", color = "black", se = F) + 
		theme_test(base_size = 15) + theme(legend.position = "none") + 
		ggtitle(title) +
		xlim(c(0,max)) + ylim(c(0,max))
	m <- ggscatter(df %>% filter(model == "Model M"), x = "actual_age", y = var, color = "#4981BF", size = 1, alpha = 0.7, 
		cor.coef = TRUE, cor.coeff.args = list(method = "pearson"), cor.coef.size = 4) +
                geom_smooth(method = "lm", color = "black", se = F) + 
		theme_test(base_size = 15) + theme(legend.position = "none") + 
		ggtitle(title) +
		xlim(c(0,max)) + ylim(c(0,max))
	ggarrange(b,f,m, ncol = 3, nrow = 1)

}

both <- import(snakemake@input[["both"]]) %>% mutate(model = "Model B")
female <- import(snakemake@input[["female"]]) %>% mutate(model = "Model F")
male <- import(snakemake@input[["male"]]) %>% mutate(model = "Model M")

data <- bind_rows(both, female, male)
print(unique(data$dataset))

data <- data %>% mutate(disease = ifelse(dataset %in% c("immage", "soundlife"), "normal", disease))

pdf(snakemake@output[["plot1"]], width = 9, height = 3)

plotPredScatter(data %>% filter(cohort == "internal"), "predicted_age", "Internal")
plotPredScatter(data %>% filter(cohort == "internal"), "c_predicted_age", "Internal")
plotPredScatter(data %>% filter(cohort == "external"), "predicted_age", "External")
plotPredScatter(data %>% filter(cohort == "external"), "c_predicted_age", "External")
plotPredScatter(data %>% filter(cohort == "external", disease == "normal"), "predicted_age", "External healthy")
plotPredScatter(data %>% filter(cohort == "external", disease == "normal"), "c_predicted_age", "External healthy")
plotPredScatter(data %>% filter(cohort == "external", disease != "normal"), "predicted_age", "External disease")
plotPredScatter(data %>% filter(cohort == "external", disease != "normal"), "c_predicted_age", "External disease")
plotPredScatter(data %>% filter(dataset == "immage"), "predicted_age", "Immage")
plotPredScatter(data %>% filter(dataset == "immage"), "c_predicted_age", "Immage")
plotPredScatter(data %>% filter(dataset == "soundlife"), "predicted_age", "Soundlife")
plotPredScatter(data %>% filter(dataset == "soundlife"), "c_predicted_age", "Soundlife")
plotPredScatter(data %>% filter(cohort == "external", disease == "normal", dataset != "immage", dataset != "soundlife"), "predicted_age", "Non-allen External healthy")
plotPredScatter(data %>% filter(cohort == "external", disease == "normal", dataset != "immage", dataset != "soundlife"), "c_predicted_age", "Non-allen External healthy")

dev.off()
