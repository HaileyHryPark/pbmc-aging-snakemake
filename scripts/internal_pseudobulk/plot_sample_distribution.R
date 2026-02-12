library(rio)
library(tidyr)
library(dplyr)
library(tidyverse)
library(ggpubr)

## Functions

sexcolor <- c("female" = "#E96D79", "male" = "#5995CB")

meta <- import(snakemake@input[["data"]]) %>% select(rowname, age, sex, dataset, ethnicity) %>% 
	mutate(ethnicity_org = ethnicity, ethnicity = ifelse(ethnicity_org %in% c("European", "Caucasian"), "White", ifelse(ethnicity_org == "Hispanic or Latin American", "Hispanic", "Asian")))

print(table(meta$ethnicity_org))

p1 <- ggplot(meta, aes(x = age)) +
        geom_histogram(aes(fill = sex), color = "white", binwidth = 5, boundary = 50) +
        scale_fill_manual(values = sexcolor, name = "Sex") +
        xlab("Age")+
        ylab("Number of samples")+
        theme_classic(base_size = 12)+
	theme(panel.grid = element_blank(), strip.text = element_text(size=17))

p2 <- ggplot(meta, aes(x = age)) +
        geom_histogram(aes(fill = ethnicity), color = "white", binwidth = 5, boundary = 50) +
        scale_fill_manual(values = c("White" = "#ACA233", "Asian" = "#53B7DA", "Hispanic" = "#DE6FE5"), name = "Ethnic group") +
        xlab("Age")+
        ylab("Number of samples")+
        theme_classic(base_size = 12)+
	theme(panel.grid = element_blank(), strip.text = element_text(size=17))

ggsave(snakemake@output[["plot1"]], ggarrange(p1,p2,ncol=1,nrow=2), width = 3.5, height = 4)

p3 <- ggplot(meta, aes(x = age)) +
        geom_histogram(aes(fill = ethnicity_org), color = "white", binwidth = 5, boundary = 50) +
        xlab("Age")+
        ylab("Number of samples")+
        theme_classic(base_size = 25)+
	theme(panel.grid = element_blank(), strip.text = element_text(size=17))

ggsave(snakemake@output[["plot2"]], p3, width = 10, height = 5)

pdf(snakemake@output[["plot_ds"]], width = 3.5, height = 4)

lapply(as.list(unique(meta$dataset)), function(ds){
	m <- meta %>% filter(dataset == ds)

	p1 <- ggplot(m, aes(x = age)) +
        	geom_histogram(aes(fill = sex), color = "white", binwidth = 5, boundary = 50) +
        	scale_fill_manual(values = sexcolor, name = "Sex") +
        	xlab("Age")+
        	ylab("Number of samples") +
		ylim(0,150) +
		xlim(18,100) +
        	theme_classic()+
		theme(panel.grid = element_blank(), strip.text = element_text(size=17))
	p2 <- ggplot(m, aes(x = age)) +
        	geom_histogram(aes(fill = ethnicity_org), color = "white", binwidth = 5, boundary = 50) +
        	xlab("Age")+
        	ylab("Number of samples") +
		ylim(0,150) +
		xlim(18,100) +
        	theme_classic()+
		theme(panel.grid = element_blank(), strip.text = element_text(size=17))
	print(ggarrange(p1, p2, ncol = 1, nrow = 2))
})

dev.off()
