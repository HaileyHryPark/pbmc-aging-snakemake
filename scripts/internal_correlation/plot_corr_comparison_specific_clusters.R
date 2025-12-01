library(rio)
library(tidyverse)
library(dplyr)
library(ggpubr)

annot_cor <- import(snakemake@input[["annot_cor"]])

cor_fci <- annot_cor %>% filter(type == "female", final_cluster == "Continuous\nincrease")
cor_mei <- annot_cor %>% filter(type == "male", final_cluster == "Early\nincrease")

mei_fci <- intersect(cor_fci$feature, cor_mei$feature)

cor_mei_fci <- bind_rows(cor_fci %>% filter(feature %in% mei_fci),
                         cor_mei %>% filter(feature %in% mei_fci))

pdf(snakemake@output[["plots"]], width = 4, height = 3)

p1 <- ggplot(cor_mei_fci, aes(x = type, y = rho)) + 
  geom_violin(aes(fill = type)) +
  geom_boxplot(fill = "white", width = 0.3) + 
#  geom_point(size = 0.5) +
#  geom_line(aes(group = feature), linewidth = 0.3) +
  stat_compare_means(paired = T, comparisons = list(c("female", "male")), label = "p", tip.length = 0, bracket.size = 0.7, vjust = -0.4) + 
  scale_fill_manual(values = c("female" = "#E15566", "male" = "#4981BF"))+
  xlab("") +
  scale_y_continuous(expand = expansion(mult = c(0,.15)))+
  theme_classic(base_size = 15) +
  theme(legend.position = "none")

p2 <- ggplot(cor_mei_fci, aes(x = type, y = r)) + 
  geom_violin(aes(fill = type)) +
  geom_boxplot(fill = "white", width = 0.3) + 
  geom_point(size = 0.5) +
  geom_line(aes(group = feature), linewidth = 0.3) +
  stat_compare_means(paired = T, comparisons = list(c("female", "male")), label = "p", tip.length = 0, bracket.size = 0.7, vjust = -0.4) + 
  scale_fill_manual(values = c("female" = "#E15566", "male" = "#4981BF"))+
  xlab("") +
  scale_y_continuous(expand = expansion(mult = c(0,.15)))+
  theme_classic(base_size = 15) +
  theme(legend.position = "none")

ggarrange(p1, p2, ncol = 2, nrow = 1)

p3 <- cor_mei_fci %>% filter(type == "female") %>% select(feature, rho, r) %>% 
  pivot_longer(!feature, names_to = "method", values_to = "correlation") %>% 
  mutate(method = ifelse(method == "rho", "spearman", "pearson")) %>% 
  ggplot(aes(x = method, y = correlation)) + 
  geom_violin(fill = "#E15566") +
  geom_boxplot(fill = "white", width = 0.3) + 
  geom_point(size = 0.5) +
  geom_line(aes(group = feature), linewidth = 0.3) +
  stat_compare_means(paired = T, comparisons = list(c("spearman", "pearson")), label = "p", tip.length = 0, bracket.size = 0.7, vjust = -0.4) + 
  xlab("") +
  scale_y_continuous(expand = expansion(mult = c(0,.15)))+
  theme_classic(base_size = 15) +
  theme(legend.position = "none")

p4 <- cor_mei_fci %>% filter(type == "male") %>% select(feature, rho, r) %>% 
  pivot_longer(!feature, names_to = "method", values_to = "correlation") %>% 
  mutate(method = ifelse(method == "rho", "spearman", "pearson")) %>% 
  ggplot(aes(x = method, y = correlation)) + 
  geom_violin(fill = "#4981BF") +
  geom_boxplot(fill = "white", width = 0.3) + 
  geom_point(size = 0.5) +
  geom_line(aes(group = feature), linewidth = 0.3) +
  stat_compare_means(paired = T, comparisons = list(c("spearman", "pearson")), label = "p", tip.length = 0, bracket.size = 0.7, vjust = -0.4) + 
  xlab("") +
  scale_y_continuous(expand = expansion(mult = c(0,.15)))+
  theme_classic(base_size = 15) +
  theme(legend.position = "none")

ggarrange(p3, p4, ncol = 2, nrow = 1)

p5 <- cor_mei_fci %>% filter(type == "female") %>% select(feature, rho.b, rho.a) %>% 
  pivot_longer(!feature, names_to = "range", values_to = "correlation") %>% 
  mutate(range = ifelse(range == "rho.b", "Before age 60", "After age 60")) %>% 
  mutate(range = factor(range, levels = c("Before age 60", "After age 60"))) %>% 
  ggplot(aes(x = range, y = correlation)) + 
  geom_violin(fill = "#E15566") +
  geom_boxplot(fill = "white", width = 0.3) + 
  geom_point(size = 0.5) +
  geom_line(aes(group = feature), linewidth = 0.3) +
  stat_compare_means(paired = T, comparisons = list(c("Before age 60", "After age 60")), label = "p", tip.length = 0, bracket.size = 0.7, vjust = -0.4) + 
  xlab("") +
  scale_y_continuous(expand = expansion(mult = c(0,.15)))+
  theme_classic(base_size = 15) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))

p6 <- cor_mei_fci %>% filter(type == "male") %>% select(feature, rho.b, rho.a) %>% 
  pivot_longer(!feature, names_to = "range", values_to = "correlation") %>% 
  mutate(range = ifelse(range == "rho.b", "Before age 60", "After age 60")) %>% 
  mutate(range = factor(range, levels = c("Before age 60", "After age 60"))) %>% 
  ggplot(aes(x = range, y = correlation)) + 
  geom_violin(fill = "#4981BF") +
  geom_boxplot(fill = "white", width = 0.3) + 
  geom_point(size = 0.5) +
  geom_line(aes(group = feature), linewidth = 0.3) +
  stat_compare_means(paired = T, comparisons = list(c("Before age 60", "After age 60")), label = "p", tip.length = 0, bracket.size = 0.7, vjust = -0.4) + 
  xlab("") +
  scale_y_continuous(expand = expansion(mult = c(0,.15)))+
  theme_classic(base_size = 15) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))

ggarrange(p5, p6, ncol = 2, nrow = 1)

p7 <- ggplot(cor_mei_fci, aes(x = type, y = rho.b)) + 
  geom_violin(aes(fill = type)) +
  geom_boxplot(fill = "white", width = 0.3) + 
#  geom_point(size = 0.5) +
#  geom_line(aes(group = feature), linewidth = 0.3) +
  stat_compare_means(paired = T, comparisons = list(c("female", "male")), label = "p", tip.length = 0, bracket.size = 0.7, vjust = -0.4) + 
  scale_fill_manual(values = c("female" = "#E15566", "male" = "#4981BF"))+
  xlab("") +
  scale_y_continuous(expand = expansion(mult = c(0,.15)))+
  theme_classic(base_size = 15) +
  theme(legend.position = "none")

p8 <- ggplot(cor_mei_fci, aes(x = type, y = rho.a)) + 
  geom_violin(aes(fill = type)) +
  geom_boxplot(fill = "white", width = 0.3) + 
#  geom_point(size = 0.5) +
#  geom_line(aes(group = feature), linewidth = 0.3) +
  stat_compare_means(paired = T, comparisons = list(c("female", "male")), label = "p", tip.length = 0, bracket.size = 0.7, vjust = -0.4) + 
  scale_fill_manual(values = c("female" = "#E15566", "male" = "#4981BF"))+
  xlab("") +
  scale_y_continuous(expand = expansion(mult = c(0,.15)))+
  theme_classic(base_size = 15) +
  theme(legend.position = "none")

ggarrange(p7, p8, ncol = 2, nrow = 1)

dev.off()
