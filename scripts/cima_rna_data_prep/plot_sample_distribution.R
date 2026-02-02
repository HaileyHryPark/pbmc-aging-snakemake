library(rio)
library(tidyr)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(svglite)

## Functions

sexcolor <- c("Female" = "#E96D79", "Male" = "#5995CB")

meta <- import(snakemake@input[["data"]]) %>% select(rowname, age, sex, dataset)

p1 <- ggplot(meta, aes(x = age)) +
        geom_histogram(aes(fill = sex), color = "white", binwidth = 5, boundary = 50) +
        scale_fill_manual(values = sexcolor, name = "Sex") +
        xlab("Age")+
        ylab("Number of samples")+
        theme_classic(base_size = 12)+
	theme(panel.grid = element_blank(), strip.text = element_text(size=17))

ggsave(snakemake@output[["plot"]], p1, width = 3.5, height = 2.2)

