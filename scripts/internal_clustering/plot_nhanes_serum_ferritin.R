library(tidyverse)
library(ggpubr)
library(svglite)
library(rio)

## NHANES 2017-2018 Data
meta_data <- haven::read_xpt(snakemake@input[["meta"]])
ferr_data <- haven::read_xpt(snakemake@input[["sf"]])

dflist <- list(meta_data[,c(1,4,5,8)], ferr_data[,c(1,3)])

nhanes_2017 <- dflist %>% purrr::reduce(full_join, by = "SEQN")
nhanes_2017$RIAGENDR <- factor(nhanes_2017$RIAGENDR, labels = c("Male","Female"))
nhanes_2017$RIDRETH3 <- factor(nhanes_2017$RIDRETH3, labels = c("MEXICAN","HISPANIC","WHITE","BLACK","ASIAN","OTHERS"))
nhanes_2017$AGEGROUP <- cut(nhanes_2017$RIDAGEYR, breaks = c(0,seq(11,80, by = 10), Inf), labels = c(paste(c(0,seq(11,80, by = 10)), seq(10, 80, by = 10),sep = "-")), right = FALSE)

export(nhanes_2017, snakemake@output[["data"]])

print(table(nhanes_2017$RIAGENDR, nhanes_2017$AGEGROUP))

## serum ferritin
p <- ggplot(nhanes_2017 %>% filter(!is.na(LBDFERSI)), aes(x = AGEGROUP, y = LBDFERSI, fill = RIAGENDR)) +
  geom_boxplot(outliers = FALSE) +
  facet_grid(~RIAGENDR) + 
  scale_fill_manual(values = c("Female" = "#E15566", "Male" = "#4981BF")) +
  theme_linedraw(base_size = 15) +
  theme(panel.grid.major = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  ylim(0,1000) +
  labs(y = "Serum Ferritin (ug/L)")

ggsave(snakemake@output[["plot"]], p, width = 5, height = 3.5)
