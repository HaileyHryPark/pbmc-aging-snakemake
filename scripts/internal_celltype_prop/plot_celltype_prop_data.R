library(rio)
library(dplyr)
library(tidyverse)
library(purrr)
library(ggpubr)


## Functions
raw <- import(snakemake@input[["raw"]])
data <- import(snakemake@input[["data"]])
celltypes <- colnames(data %>% select(-c(donor_id, age, sex, dataset, ethnicity)))

pdf(snakemake@output[["plot1"]], width = 9, height = 5)
lapply(celltypes, function(f){

	plots <- lapply(list("both","female","male"), function(g){ 
		df <- data
		if(g != "both"){
			df <- df %>% filter(sex == g)
		}

		p1 <- ggscatter(df, x = "age", y = f, alpha = 0.7, size = 0.5, cor.coef = T, cor.method = "spearman") +
			geom_smooth(method = "loess", method.args = list(degree = 1), se = F) +
			ggtitle(paste(f, g)) +
			theme_test() + theme(legend.position = "none")
		p2 <- ggscatter(df, x = "age", y = f, color = "dataset", alpha = 0.7, size = 0.5, cor.coef = T, cor.method = "spearman") +
			ggtitle(paste(f, g)) +
			theme_test() + theme(legend.position = "none")
		p <- ggarrange(p1, p2, ncol = 1, nrow = 2, widths = c(1,1))
		return(p)
	})
	ggarrange(plotlist = plots, ncol = 3, nrow = 1)
})
dev.off()

pdf(snakemake@output[["plot2"]], width = 9, height = 5)
lapply(celltypes, function(f){

	plots <- lapply(list("both","female","male"), function(g){ 
		df <- raw 
		if(g != "both"){
			df <- df %>% filter(sex == g)
		}

		p1 <- ggscatter(df, x = "age", y = f, alpha = 0.7, size = 0.5, cor.coef = T, cor.method = "spearman") +
			geom_smooth(method = "loess", method.args = list(degree = 1), se = F) +
			ggtitle(paste(f, g)) +
			theme_test() + theme(legend.position = "none")
		p2 <- ggscatter(df, x = "age", y = f, color = "dataset", alpha = 0.7, size = 0.5, cor.coef = T, cor.method = "spearman") +
			ggtitle(paste(f, g)) +
			theme_test() + theme(legend.position = "none")
		p <- ggarrange(p1, p2, ncol = 1, nrow = 2, widths = c(1,1))
		return(p)
	})
	ggarrange(plotlist = plots, ncol = 3, nrow = 1)
})
dev.off()
