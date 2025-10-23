library(rio)
library(tidyr)
library(dplyr)
library(tidyverse)
library(ggpubr)

## Functions

sexcolor <- c("female" = "#E96D79", "male" = "#5995CB")

meta <- import(snakemake@input[["data"]]) %>% select(rowname, age, sex, dataset, ethnicity) %>% 
	mutate(ethnicity = ifelse(ethnicity %in% c("European", "Caucasian"), "European/Caucasian", ethnicity))

p1 <- ggplot(meta, aes(x = age)) +
        geom_histogram(aes(fill = sex), color = "white", binwidth = 5, boundary = 50) +
        scale_fill_manual(values = sexcolor, name = "gender") +
        xlab("Age")+
        ylab("Number of samples")+
        theme_classic()

p2 <- ggplot(meta, aes(x = age)) +
        geom_histogram(aes(fill = ethnicity), color = "white", binwidth = 5, boundary = 50) +
        xlab("Age")+
        ylab("Number of samples")+
        theme_classic()

ggsave(snakemake@output[["plot"]], ggarrange(p1,p2,ncol=1,nrow=2), width = 4, height = 4)

pdf(snakemake@output[["plot_ds"]], width = 3.5, height = 4)

lapply(as.list(unique(meta$dataset)), function(ds){
	m <- meta %>% filter(dataset == ds)

	p1 <- ggplot(m, aes(x = age)) +
        	geom_histogram(aes(fill = sex), color = "white", binwidth = 5, boundary = 50) +
        	scale_fill_manual(values = sexcolor, name = "gender") +
        	xlab("Age")+
        	ylab("Number of samples") +
		ylim(0,150) +
		xlim(18,100) +
        	theme_classic()
	p2 <- ggplot(m, aes(x = age)) +
        	geom_histogram(aes(fill = ethnicity), color = "white", binwidth = 5, boundary = 50) +
        	xlab("Age")+
        	ylab("Number of samples") +
		ylim(0,150) +
		xlim(18,100) +
        	theme_classic()
	print(ggarrange(p1, p2, ncol = 1, nrow = 2))
})

dev.off()
