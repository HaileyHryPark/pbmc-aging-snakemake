library(rio)
library(dplyr)
library(tidyverse)
library(ggpubr)

all_res <- import(snakemake@input[["res"]]) 

pdf(snakemake@output[["plot2"]], width = 4.5, height = 3.5)
lapply(list(5,10,15,20,25), function(b){

	q <- all_res %>% filter(bucket == b) %>% select(-bucket)
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
                        title = ct,
                        x = "Age (years)",
                        y = paste0("Number of DEGs (Window ", b, ")"),
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
                        title = ct,
                        x = "Age Threshold",
                        y = paste0("Number of DEGs (Window ", b, ")"),
                        color = "Gender"
		) +
                theme_test(base_size = 14)+
                theme(panel.grid = element_blank(), plot.title = element_text(hjust=0.5))
	
if(ct == "All"){
		plot(p1 + ylim(0,4000))
	        plot(p2 + ylim(0,4000))
	}else{
		plot(p1)
	        plot(p2)
	}
})
})
dev.off()

pdf(snakemake@output[["plot1"]], width = 4.5, height = 3.5)

q <- all_res %>% filter(bucket == 20) %>% select(-bucket)
q <- q %>% separate(variable, c("celltype", "gene"), remove = FALSE, sep ="\\.")

lapply(list(0.05,0.01,0.001,0.0001), function(b){

	q_long <- q %>% 
		pivot_longer(cols = -c(gender, variable, celltype, gene), names_to = "age_threshold", values_to = "qvalue") %>%
		mutate(age_threshold = as.numeric(gsub("X", "", age_threshold)), gender = factor(gender, levels = c("Both","Female","Male")))

	all_combinations <- q_long %>% distinct(gender, celltype, age_threshold)

	deg_summary <- q_long %>%
	  filter(qvalue < b) %>%
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
                labs(
                        title = ct,
                        x = "Age (years)",
                        y = paste0("Number of DEGs (q < ", b, ")"),
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
                        title = ct,
                        x = "Age Threshold",
                        y = paste0("Number of DEGs (q < ", b, ")"),
                        color = "Gender"
		) +
                theme_test(base_size = 14)+
                theme(panel.grid = element_blank(), plot.title = element_text(hjust=0.5))
        
	if(ct == "All"){
		plot(p1 + ylim(0,3800))
	        plot(p2 + ylim(0,3800))
	}else{
		plot(p1)
	        plot(p2)
	}
})
})
dev.off()

