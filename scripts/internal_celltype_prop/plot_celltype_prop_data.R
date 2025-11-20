library(rio)
library(dplyr)
library(tidyverse)
library(purrr)
library(ggpubr)


## Functions
# Function to compute Spearman correlation
get_corr <- function(data, cell) {
  test <- suppressWarnings(cor.test(data[[cell]], data$age, method = "spearman"))
  tibble(
    celltype = cell,
    rho = test$estimate,
    pval = test$p.value
  )
}


raw <- import(snakemake@input[["raw"]])
data <- import(snakemake@input[["data"]])
celltypes <- colnames(data %>% select(-c(donor_id, age, sex, dataset, ethnicity)))

# Compute correlations for each subgroup
res_both   <- bind_rows(lapply(celltypes, function(c) get_corr(raw, c))) %>%
  rename(rho_both = rho, p_both = pval)

res_female <- bind_rows(lapply(celltypes, function(c) get_corr(raw %>% filter(sex == "female"), c))) %>%
  rename(rho_female = rho, p_female = pval)

res_male   <- bind_rows(lapply(celltypes, function(c) get_corr(raw %>% filter(sex == "male"), c))) %>%
  rename(rho_male = rho, p_male = pval)

# Merge all tables
final_table <- res_both %>%
  left_join(res_female, by = "celltype") %>%
  left_join(res_male,   by = "celltype")

export(final_table, snakemake@output[["table"]])


pdf(snakemake@output[["plot1"]], width = 9, height = 5)
lapply(celltypes, function(f){

	plots <- lapply(list("both","female","male"), function(g){ 
		df <- data
		if(g != "both"){
			df <- df %>% filter(sex == g)
		}

		p1 <- ggscatter(df, x = "age", y = f, alpha = 0.7, size = 0.5, cor.coef = T, cor.method = "spearman") +
			geom_smooth(method = "loess", method.args = list(degree = 1), se = F) +
			ggtitle(paste(f, g)) +
			theme_test() + theme(legend.position = "none")
		p2 <- ggscatter(df, x = "age", y = f, color = "dataset", alpha = 0.7, size = 0.5, cor.coef = T, cor.method = "spearman") +
			ggtitle(paste(f, g)) +
			theme_test() + theme(legend.position = "none")
		p <- ggarrange(p1, p2, ncol = 1, nrow = 2, widths = c(1,1))
		return(p)
	})
	ggarrange(plotlist = plots, ncol = 3, nrow = 1)
})
dev.off()

pdf(snakemake@output[["plot2"]], width = 9, height = 5)
lapply(celltypes, function(f){

	plots <- lapply(list("both","female","male"), function(g){ 
		df <- raw 
		if(g != "both"){
			df <- df %>% filter(sex == g)
		}

		p1 <- ggscatter(df, x = "age", y = f, alpha = 0.7, size = 0.5, cor.coef = T, cor.method = "spearman") +
			geom_smooth(method = "loess", method.args = list(degree = 1), se = F) +
			ggtitle(paste(f, g)) +
			theme_test() + theme(legend.position = "none")
		p2 <- ggscatter(df, x = "age", y = f, color = "dataset", alpha = 0.7, size = 0.5, cor.coef = T, cor.method = "spearman") +
			ggtitle(paste(f, g)) +
			theme_test() + theme(legend.position = "none")
		p <- ggarrange(p1, p2, ncol = 1, nrow = 2, widths = c(1,1))
		return(p)
	})
	ggarrange(plotlist = plots, ncol = 3, nrow = 1)
})
dev.off()
