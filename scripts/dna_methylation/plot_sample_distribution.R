library(rio)
library(tidyverse)
library(ggpubr)

metadata <- import(snakemake@input[["metadata"]])

gender_cols <- c("female" = "#E96D79", "male" = "#5995CB")

p1 <- ggplot(metadata, aes(x = Age)) +
  geom_histogram(aes(fill = Gender), color = "white", binwidth = 5, boundary = 50) +
  scale_fill_manual(values = gender_cols, name = "Gender") +
  xlab("Age")+
  ylab("Number of samples")+
  theme_classic(base_size = 14)+
  theme(panel.grid = element_blank(), strip.text = element_text(size=17))

ggsave(snakemake@output[["plot"]], p1, width = 4, height = 2.2)

