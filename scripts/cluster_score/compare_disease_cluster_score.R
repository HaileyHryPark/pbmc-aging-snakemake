.libPaths(new = "/scratch/users/nus/e0859928/Snakemake/onek1k-analysis-snakemake/resources/r_package")

library(rio)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(ggfortify)
library(fmsb)

age_disease_col <- c("grey", "#B24745")

## Functions
PlotBoxplot <- function(data, ds, dis){

	p1 <- ggplot(data %>% filter(dataset == ds) %>% mutate(disease = factor(disease, levels = c("normal", dis))), aes(x = disease, y = score)) + 
		geom_boxplot(aes(color = disease), width = 0.6) +
		geom_point(position = position_jitter(width = 0.2), size = 0.5, alpha = 1, color = "black")+
		scale_color_manual(values = age_disease_col)+
		stat_compare_means(method = "wilcox", paired = FALSE, label = "p", comparisons = list(c("normal", dis)), tip.length = 0, bracket.size = 0.7, vjust = -0.4) +
		facet_grid(cluster) +
		scale_y_continuous(expand = expansion(mult = c(0.05,.15))) +
		theme_linedraw(base_size = 15) +
		theme(panel.grid = element_blank(), legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))

	plot(p1)
}

PlotRadarplot <- function(data){

	med_df1 <- data %>% filter(dataset == ds) %>% mutate(disease = factor(disease, levels = c("normal", dis))) %>% group_by(disease, cluster) %>% summarise(med = median(score, na.rm = T)) %>% 
		pivot_wider(names_from = cluster, values_from = med) %>% column_to_rownames("disease")
	print(med_df1)
	
	n <- ncol(med_df1)
	med_p1 <- radarchart(rbind(rep(0.5,n), rep(-0.1,n), med_df1), axistype = 0, pcol = age_disease_col, plwd=5 , vlcex=0.8, cglcol="grey", cglty=1, axislabcol="grey")
	plot(med_p1)
}


## Main
g <- snakemake@params[["gender"]]

print(c)
print(g)

if(g == "female"){
	cluster_levels <- c("Early\nincrease", "Early\ndecrease", "Continuous\ndecrease", "Continuous\nincrease", "Late\nincrease")
}else if(g == "male"){
	cluster_levels <- c("Early\nincrease", "Early\ndecrease", "Continuous\ndecrease", "Irregular\nchange")
}else{
	cluster_levels <- c("Early\nincrease", "Early\ndecrease", "Continuous\ndecrease", "Irregular\nchange", "Late\nincrease")
}

data <- import(snakemake@input[["data"]]) %>% filter(cohort == c)

score <- import(snakemake@input[["score"]]) %>% filter(sample_id %in% data$sample_id)

data <- data %>% right_join(score, by = "sample_id") %>% 
	mutate(cluster = factor(cluster, levels = cluster_levels))

print(data[!complete.cases(data),])

pdf(snakemake@output[["boxplot"]], width = 3*length(unique(data$cluster)), height = 9)

PlotBoxplot(data, "")

dev.off()

pdf(snakemake@output[["radarplot"]], width = 12, height = 6)

PlotRadarplot(data)
if(c == "external"){
	PlotRadarplot(data %>% filter(disease == "normal"))
	PlotRadarplot(data %>% filter(disease != "normal"))
}

dev.off()
