library(ggpubr)
library(rio)
library(tidyverse)
library(svglite)

df <- import(snakemake@input[["cellcount"]])

p1 <- ggplot(df, aes(x = dataset, y = n, fill = dataset)) +
  geom_violin(color = "black", width = 0.6) +
  geom_boxplot(fill = "white", width = 0.2) +
  theme_linedraw(base_size = 15) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), panel.grid.major = element_blank(), legend.position = "none") +
  labs(x = "", y = "Total cell number by donor")

p2 <- ggplot(df, aes(x = dataset, y = log10(n), fill = dataset)) +
  geom_violin(color = "black", width = 0.6) +
  geom_boxplot(fill = "white", width = 0.2) +
  theme_linedraw(base_size = 15) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), panel.grid.major = element_blank(), legend.position = "none") +
  labs(x = "", y = "Total cell number by donor (log10)")

ggsave(snakemake@output[["plot"]], ggarrange(p1, p2, nrow = 1, ncol = 2), width = 6, height = 4)
