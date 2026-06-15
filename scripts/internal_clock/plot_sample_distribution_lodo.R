library(rio)
library(tidyr)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(svglite)

## Functions

sexcolor <- c("female" = "#E96D79", "male" = "#5995CB")
dataset_cols <- c("onek1k" = "#1b9e77", "aida" = "#d95f02", "perez" = "#7570b3", "marina" = "#e7298a")

meta <- import(snakemake@input[["data"]]) %>% select(rowname, age, sex, dataset, ethnicity) %>% 
	mutate(lodo_group = ifelse(dataset %in% c("perez", "marina"), "perez_marina", dataset))

print(table(meta$lodo_group))


plotLODOdistribution <- function(df){

p <- lapply(as.list(unique(df$lodo_group)), function(g){

	m <- df %>% filter(lodo_group != g)

	p1 <- ggplot(m, aes(x = age)) +
        	geom_histogram(aes(fill = sex), color = "white", binwidth = 5, boundary = 50) +
        	scale_fill_manual(values = sexcolor, name = "Sex") +
        	xlab("Age")+
        	ylab("Number of samples") +
		ggtitle(g) +
		ylim(0,150) +
		xlim(18,100) +
        	theme_classic()+
		theme(panel.grid = element_blank(), strip.text = element_text(size=17))
	p2 <- ggplot(m, aes(x = age)) +
        	geom_histogram(aes(fill = dataset), color = "white", binwidth = 5, boundary = 50) +
        	scale_fill_manual(values = dataset_cols, name = "Dataset") +
        	xlab("Age")+
        	ylab("Number of samples") +
		ggtitle(g) +
		ylim(0,150) +
		xlim(18,100) +
        	theme_classic()+
		theme(panel.grid = element_blank(), strip.text = element_text(size=17))
	return(ggarrange(p1, p2, ncol = 1, nrow = 2))
})

return(ggarrange(plotlist=p, ncol = length(unique(df$lodo_group)), nrow = 1))

}

ggsave(snakemake@output[["plot1"]], plotLODOdistribution(meta %>% filter(sex == "female")), width = 10, height = 4.5)
ggsave(snakemake@output[["plot2"]], plotLODOdistribution(meta %>% filter(sex == "male")), width = 10, height = 4.5)
ggsave(snakemake@output[["plot3"]], plotLODOdistribution(meta), width = 10, height = 4.5)
