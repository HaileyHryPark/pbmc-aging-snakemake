library(rio)
library(dplyr)
library(tidyverse)
library(ggpubr)

all_res <- import(snakemake@input[["q"]]) 

pdf(snakemake@output[["plots"]], width = 4.5, height = 3.5)
lapply(as.list(unique(all_res$dataset)), function(ds){

	q <- all_res %>% filter(dataset == ds) %>% select(-dataset)
	q <- q %>% separate(variable, c("celltype", "gene"), remove = FALSE, sep ="\\.")
	q_long <- q %>% 
		pivot_longer(cols = -c(gender, variable, celltype, gene), names_to = "age_threshold", values_to = "qvalue") %>%
		mutate(age_threshold = as.numeric(gsub("X", "", age_threshold)), gender = factor(gender, levels = c("Both","Female","Male")))

	all_combinations <- q_long %>% distinct(gender, celltype, age_threshold)

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
	
	deg_plot_data <- bind_rows(deg_summary, deg_total) %>% mutate(celltype = factor(celltype, levels = c("All","CD4 T","NK","CD8 T","B","other T","Mono")))

	lapply(as.list(unique(deg_plot_data$celltype)), function(ct){	

	p1 <- deg_plot_data %>% filter(celltype == ct) %>% ggplot( aes(x = age_threshold, y = n_DEGs, color = gender)) +
                geom_line(linewidth = 1) +
                geom_point(size = 1.5) +
                scale_color_manual(values = c("Both" = "grey", "Female" = "#E15566", "Male" = "#4981BF")) +
                #facet_wrap(~ celltype)+
                labs(
                        title = paste(ct, ds),
                        x = "Age (years)",
                        y = "Number of DEGs (q < 0.05)",
                        color = "Gender"
                ) +
                theme_test(base_size = 14)+
                theme(panel.grid = element_blank(), plot.title = element_text(hjust=0.5))

        p2 <- deg_plot_data %>% filter(celltype == ct, gender != "Both") %>% ggplot( aes(x = age_threshold, y = n_DEGs, color = gender)) +
                geom_line(linewidth = 1) +
                geom_point(size = 1.5) +
                scale_color_manual(values = c("Female" = "#E15566", "Male" = "#4981BF")) +
                #facet_wrap(~ celltype)+
                labs(
                        title = paste(ct, ds),
                        x = "Age (years)",
                        y = "Number of DEGs (q < 0.05)",
                        color = "Gender"
		) +
                theme_test(base_size = 14)+
                theme(panel.grid = element_blank(), plot.title = element_text(hjust=0.5))
	
	if(ct == "All"){
		plot(p1 + ylim(0,4100))
	        plot(p2 + ylim(0,4100))
	}else{
		plot(p1)
	        plot(p2)
	}
})
})
dev.off()

