.libPaths(new = "/scratch/users/nus/e0859928/Snakemake/onek1k-analysis-snakemake/resources/r_package")

library(rio)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(ggfortify)
library(fmsb)

age_accel1_col <- c("Decelerated" = "#79AF97", "Accelerated" = "#B24745")
age_accel2_col <- c("Decelerated" = "#79AF97", "Intermediate" = "grey", "Accelerated" = "#B24745")

## Functions
PlotBoxplot <- function(data){

	p1 <- ggplot(data, aes(x = age_accel1, y = score)) + 
		geom_boxplot(aes(color = age_accel1), width = 0.6) +
		geom_point(position = position_jitter(width = 0.2), size = 0.5, alpha = 1, color = "black")+
		scale_color_manual(values = age_accel1_col)+
		stat_compare_means(method = "wilcox", paired = FALSE, label = "p", comparisons = list(c("Decelerated", "Accelerated")), tip.length = 0, bracket.size = 0.7, vjust = -0.4) +
		facet_grid(age_group~cluster) +
		scale_y_continuous(expand = expansion(mult = c(0.05,.15))) +
		theme_linedraw(base_size = 15) +
		theme(panel.grid = element_blank(), legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))

	p2 <- ggplot(data, aes(x = age_accel2, y = score)) + 
		geom_boxplot(aes(color = age_accel2), width = 0.6) +
		geom_point(position = position_jitter(width = 0.2), size = 0.5, alpha = 1, color = "black")+
		scale_color_manual(values = age_accel2_col)+
		stat_compare_means(method = "wilcox", paired = FALSE, label = "p", comparisons = list(c("Decelerated", "Intermediate"), c("Intermediate", "Accelerated"), c("Decelerated", "Accelerated")), tip.length = 0, bracket.size = 0.7, vjust = -0.4, step.increase = 0.2) +
		facet_grid(age_group~cluster) +
		scale_y_continuous(expand = expansion(mult = c(0.05,.2))) +
		theme_linedraw(base_size = 15) +
		theme(panel.grid = element_blank(), legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))

	plot(p1)
	plot(p2)
}

PlotRadarplot <- function(data){

	med_df1 <- data %>% group_by(age_accel1, cluster) %>% summarise(med = median(score, na.rm = T)) %>% 
		pivot_wider(names_from = cluster, values_from = med) %>% column_to_rownames("age_accel1")
	print(med_df1)
	
	med_df2 <- data %>% group_by(age_accel2, cluster) %>% summarise(med = median(score, na.rm = T)) %>% 
		pivot_wider(names_from = cluster, values_from = med) %>% column_to_rownames("age_accel2")
	print(med_df2)
	
	n <- ncol(med_df1)
	med_p1 <- radarchart(rbind(rep(0.5,n), rep(-0.1,n), med_df1), axistype = 0, pcol = age_accel1_col, plwd=5 , vlcex=0.8, cglcol="grey", cglty=1, axislabcol="grey")
	med_p2 <- radarchart(rbind(rep(0.5,n) , rep(-0.1,n) , med_df2), axistype = 0, pcol = age_accel2_col, plwd=5 , vlcex=0.8, cglcol="grey", cglty=1, axislabcol="grey")
	ggarrange(med_p1, med_p2, ncol = 2, nrow = 1)

}


## Main
c <- snakemake@params[["cohort"]]
g <- snakemake@params[["gender"]]

print(c)
print(g)

if(g == "female"){
	cluster_levels <- c("Early\nincrease", "Early\ndecrease", "Continuous\ndecrease", "Continuous\nincrease", "Late\nincrease")
}else if(g == "male"){
	cluster_levels <- c("Early\nincrease", "Early\ndecrease", "Continuous\ndecrease", "Early\nfluctuation")
}else{
	cluster_levels <- c("Early\nincrease", "Early\ndecrease", "Continuous\ndecrease", "Early\nfluctuation", "Late\nincrease")
}

data <- import(snakemake@input[["data"]]) %>% filter(cohort == c)

score <- import(snakemake@input[["score"]]) %>% filter(sample_id %in% data$sample_id)

data <- data %>% right_join(score, by = "sample_id") %>% 
	mutate(cluster = factor(cluster, levels = cluster_levels),
		age_group = factor(age_group, levels = c("<40", "40-60", ">60")), 
		age_accel1 = factor(age_accel1, levels = c("Decelerated", "Accelerated")), 
		age_accel2 = factor(age_accel2, levels = c("Decelerated", "Intermediate", "Accelerated")))

print(data[!complete.cases(data),])

pdf(snakemake@output[["boxplot"]], width = 3*length(unique(data$cluster)), height = 9)

PlotBoxplot(data)
if(c == "external"){
	PlotBoxplot(data %>% filter(disease == "normal"))
	PlotBoxplot(data %>% filter(disease != "normal"))
}

dev.off()

pdf(snakemake@output[["radarplot"]], width = 12, height = 6)

PlotRadarplot(data)
if(c == "external"){
	PlotRadarplot(data %>% filter(disease == "normal"))
	PlotRadarplot(data %>% filter(disease != "normal"))
}

dev.off()
