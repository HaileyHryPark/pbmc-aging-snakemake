library(rio)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(svglite)

q <- import(snakemake@input[["q"]])

q <- q %>% separate(variable, c("celltype", "gene"), remove = FALSE, sep ="\\.")

q_long <- q %>%
  pivot_longer(cols = -c(gender, variable, celltype, gene), names_to = "age_threshold", values_to = "qvalue") %>%
  mutate(age_threshold = as.numeric(gsub("X", "", age_threshold)), gender = factor(gender, levels = c("Both","Female","Male")))

print(head(q_long))

all_combinations <- q_long %>% 
  distinct(gender, celltype, age_threshold)

deg_summary <- q_long %>%
  filter(qvalue < 0.05) %>%
  group_by(gender, celltype, age_threshold) %>%
  summarise(n_DEGs = n(), .groups = "drop") %>%
  right_join(all_combinations, by = c("gender", "celltype", "age_threshold")) %>%
  mutate(n_DEGs = replace_na(n_DEGs, 0))

deg_total <- deg_summary %>%
  group_by(gender, age_threshold) %>%
  summarise(n_DEGs = sum(n_DEGs), .groups = "drop") %>%
  mutate(celltype = "All")

deg_plot_data <- bind_rows(deg_summary, deg_total) %>% mutate(celltype = factor(celltype, levels = c("All","CD4 T","CD8 T","B","NK","Mono")))

p <- deg_plot_data %>% filter(celltype == "All") %>% ggplot(aes(x = age_threshold, y = n_DEGs, color = gender)) +
	geom_line(linewidth = 1) +
	geom_point(size = 1.5) +
	scale_color_manual(values = c("Both" = "grey", "Female" = "#E15566", "Male" = "#4981BF")) +
  	labs(
    		x = "Age (years)",
    		y = "Number of DEFs",
    		color = "Sex"
  	) +
	xlim(20, 100) +
	ylim(0, 4050) +
  	theme_classic(base_size = 15)
		
ggsave(snakemake@output[["plot1"]], p, width = 5.3, height = 3.5)


sexcolor <- c("female" = "#E96D79", "male" = "#5995CB")

data <- import(snakemake@input[["data"]]) %>% select(rowname, age, sex, dataset) %>% mutate(age = as.integer(age))

p1 <- ggplot(data, aes(x = age)) +
        geom_histogram(aes(fill = sex), color = "white", binwidth = 5, boundary = 50) +
        scale_fill_manual(values = sexcolor, name = "Sex") +
        xlab("Age")+
        ylab("Number of samples")+
        ylim(0,200) +
        xlim(18,100) +
        theme_classic()+
        theme(panel.grid = element_blank(), strip.text = element_text(size=17))

ggsave(snakemake@output[["plot2"]], p1, width = 3.5, height = 2.2)

