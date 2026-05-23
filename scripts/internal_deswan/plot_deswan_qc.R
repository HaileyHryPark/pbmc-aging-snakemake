library(rio)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(svglite)

## q res
all_res <- import(snakemake@input[["q1"]]) %>% mutate(dataset = "all")
ds_res <- import(snakemake@input[["q2"]])

res <- bind_rows(all_res, ds_res)
print(head(res))
print(tail(res))

## sample number by threshold
meta <- import(snakemake@input[["meta"]]) %>% select(rowname, age, sex, dataset, ethnicity)

print(head(meta))

plots <- lapply(as.list(unique(res$dataset)), function(ds){

	q <- res %>% filter(dataset == ds) %>% select(-dataset)
	q <- q %>% separate(variable, c("celltype", "gene"), remove = FALSE, sep ="\\.")
	q_long <- q %>% 
		pivot_longer(cols = -c(gender, variable, celltype, gene), names_to = "age_threshold", values_to = "qvalue") %>%
		mutate(age_threshold = as.numeric(gsub("X", "", age_threshold)), gender = factor(gender, levels = c("Both","Female","Male")))

	all_combinations <- q_long %>% distinct(gender, celltype, age_threshold)

	deg_total <- q_long %>%
	  filter(qvalue < 0.05, !is.na(qvalue)) %>%
	  group_by(gender, age_threshold) %>%
	  summarise(n_DEGs = n(), .groups = "drop")

	## donor count
	meta_ds <- if (ds == "all") meta else meta %>% filter(dataset == ds)

	thresholds <- unique(deg_total$age_threshold)

	donor_counts <- expand.grid(age_threshold = thresholds, gender = unique(deg_total$gender)) %>%
	  rowwise() %>%
	  mutate(
            n_donors = {
            lower <- age_threshold - 10
            upper <- age_threshold + 10

            df <- meta_ds %>% filter(age >= lower, age <= upper)

            if (gender == "Female") df <- df %>% filter(sex == "female")
            if (gender == "Male")   df <- df %>% filter(sex == "male")

            n_distinct(df$rowname)
            }
          ) %>%
          ungroup()

	plot_df <- deg_total %>%
		left_join(donor_counts, by = c("age_threshold", "gender"))
	
	all <- ggplot(plot_df, aes(n_donors, n_DEGs, color = gender)) +
	    geom_point(size = 1) +
	    #geom_smooth(method = "lm", se = FALSE) +
	    scale_color_manual(values = c("Both"="grey","Female"="#E15566","Male"="#4981BF")) +
	    labs(title = ds, x = "Number of donors", y = "Number of DEGs") +
	    theme_test() +
	    theme(panel.grid = element_blank(), plot.title = element_text(hjust=0.5))

	female <- ggplot(plot_df %>% filter(gender == "Female"), aes(n_donors, n_DEGs, color = gender, label = age_threshold)) +
	    geom_point(size = 1) +
	    geom_text(hjust=-1,vjust=1,color="black") + 
	    scale_color_manual(values = c("Both"="grey","Female"="#E15566","Male"="#4981BF")) +
	    labs(title = ds, x = "Number of donors", y = "Number of DEGs") +
	    theme_test() +
	    theme(panel.grid = element_blank(), plot.title = element_text(hjust=0.5))

	male <- ggplot(plot_df %>% filter(gender == "Male"), aes(n_donors, n_DEGs, color = gender, label = age_threshold)) +
	    geom_point(size = 1) +
	    geom_text(hjust=-1,vjust=1,color="black") + 
	    scale_color_manual(values = c("Both"="grey","Female"="#E15566","Male"="#4981BF")) +
	    labs(title = ds, x = "Number of donors", y = "Number of DEGs") +
	    theme_test() +
	    theme(panel.grid = element_blank(), plot.title = element_text(hjust=0.5))

	both <- ggplot(plot_df %>% filter(gender == "Both"), aes(n_donors, n_DEGs, color = gender, label = age_threshold)) +
	    geom_point(size = 1) +
	    geom_text(hjust=-1,vjust=1,color="black") + 
	    scale_color_manual(values = c("Both"="grey","Female"="#E15566","Male"="#4981BF")) +
	    labs(title = ds, x = "Number of donors", y = "Number of DEGs") +
	    theme_test() +
	    theme(panel.grid = element_blank(), plot.title = element_text(hjust=0.5))

	pl <- ggarrange(all, female, male, both, ncol = 4, nrow = 1)
})

p <- ggarrange(plotlist = plots, ncol = 1, nrow = 5)
ggsave(snakemake@output[["plots"]], p, width = 20, height = 18)
