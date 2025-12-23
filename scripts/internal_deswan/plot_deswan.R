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

pdf(snakemake@output[["plots"]], width = 4.8, height = 3.5)

lapply(as.list(unique(deg_plot_data$celltype)), function(ct){
	p1 <- deg_plot_data %>% filter(celltype == ct) %>% ggplot( aes(x = age_threshold, y = n_DEGs, color = gender)) +
  		geom_line(linewidth = 1) +
		geom_point(size = 1.5) +
  		scale_color_manual(values = c("Both" = "grey", "Female" = "#E15566", "Male" = "#4981BF")) +
  		labs(
    			title = ct, 
    			x = "Age (years)",
    			y = "Number of DEGs",
    			color = "Gender"
  		) +
  		theme_classic(base_size = 15)+
  		theme(panel.grid = element_blank(), plot.title = element_text(hjust=0.5), strip.text = element_text(size=17))
		
	p2 <- deg_plot_data %>% filter(celltype == ct, gender != "Both") %>% ggplot( aes(x = age_threshold, y = n_DEGs, color = gender)) +
  		geom_line(linewidth = 1) +
		geom_point(size = 1.5) +
  		scale_color_manual(values = c("Female" = "#E15566", "Male" = "#4981BF")) +
  		labs(
    			title = ct, 
    			x = "Age Threshold",
    			y = "Number of DEGs",
    			color = "Gender"
  		) +
  		theme_classic(base_size = 15)+
  		theme(panel.grid = element_blank(), plot.title = element_text(hjust=0.5), strip.text = element_text(size=17))
	plot(p1)
	plot(p2)	
	plot(p1 + ylim(0,4050))
	plot(p2 + ylim(0,4050))
})

dev.off()

pdf(snakemake@output[["plots2"]], width = 16.5, height = 3.7)

	p1 <- deg_plot_data %>% filter(celltype != "All") %>% ggplot(aes(x = age_threshold, y = n_DEGs, color = gender)) +
  		geom_line(linewidth = 1) +
		geom_point(size = 1.5) +
  		scale_color_manual(values = c("Both" = "grey", "Female" = "#E15566", "Male" = "#4981BF")) +
  		facet_wrap(~ celltype, nrow = 1, ncol = 5)+
  		labs(
    			x = "Age (years)",
    			y = "Number of DEGs",
    			color = "Gender"
  		) +
  		theme_linedraw(base_size = 15)+
  		theme(panel.grid = element_blank(), strip.text = element_text(size=17))
		
	p2 <- deg_plot_data %>% filter(celltype != "All") %>% ggplot(aes(x = age_threshold, y = n_DEGs, color = gender)) +
  		geom_line(linewidth = 1) +
		geom_point(size = 1.5) +
  		scale_color_manual(values = c("Both" = "grey", "Female" = "#E15566", "Male" = "#4981BF")) +
  		facet_wrap(~ celltype, nrow = 1, ncol = 5, scales = "free")+
  		labs(
    			x = "Age (years)",
    			y = "Number of DEGs",
    			color = "Gender"
  		) +
  		theme_linedraw(base_size = 15)+
  		theme(panel.grid = element_blank(), strip.text = element_text(size=17), legend.position = "bottom")
	plot(p1)
	plot(p2)

dev.off()

