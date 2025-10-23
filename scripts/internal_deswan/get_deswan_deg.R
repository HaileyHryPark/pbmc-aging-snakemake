library(rio)
library(dplyr)
library(tidyverse)
library(ggpubr)

q <- import(snakemake@input[["q"]])

q <- q %>% separate(variable, c("celltype", "gene"), remove = FALSE, sep ="\\.")

q_long <- q %>%
  pivot_longer(cols = -c(gender, variable, celltype, gene), names_to = "age_threshold", values_to = "qvalue") %>%
  mutate(age_threshold = as.numeric(gsub("X", "", age_threshold)), gender = factor(gender, levels = c("Both","Female","Male")))

print(head(q_long))

q_long <- q_long %>% as.data.frame() %>% filter(qvalue < 0.05)
export(q_long %>% filter(gender == "Both"), snakemake@output[["degb"]])
export(q_long %>% filter(gender == "Female"), snakemake@output[["degf"]])
export(q_long %>% filter(gender == "Male"), snakemake@output[["degm"]])




