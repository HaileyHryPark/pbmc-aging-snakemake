library(rio)
library(tidyr)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(svglite)

## Functions

sexcolor <- c("female" = "#E96D79", "male" = "#5995CB")
datasetcolor <- c("immage" = "#71a659", "soundlife" = "#8975ca", "marina" = "#c5783e", "perez" = "#cb5683")

meta <- import(snakemake@input[["data"]]) %>% select(rowname, age, sex, dataset) %>% mutate(age = as.integer(age))
print(head(meta))

p1 <- ggplot(meta, aes(x = age)) +
        geom_histogram(aes(fill = sex), color = "white", binwidth = 5, boundary = 50) +
        scale_fill_manual(values = sexcolor, name = "Sex") +
        xlab("Age")+
        ylab("Number of samples")+
	ylim(0,150) +
	xlim(18,100) +
        theme_classic()+
	theme(panel.grid = element_blank(), strip.text = element_text(size=17))

p2 <- ggplot(meta, aes(x = age)) +
        geom_histogram(aes(fill = dataset), color = "white", binwidth = 5, boundary = 50) +
        scale_fill_manual(values = datasetcolor, name = "Dataset") +
        xlab("Age")+
        ylab("Number of samples")+
	ylim(0,150) +
	xlim(18,100) +
        theme_classic()+
	theme(panel.grid = element_blank(), strip.text = element_text(size=17))

p <- ggarrange(p1, p2, ncol = 1, nrow = 2)

ggsave(snakemake@output[["plot"]], p, width = 3.5, height = 5)

