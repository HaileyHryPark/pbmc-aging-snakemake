library(rio)
library(dplyr)
library(tidyverse)
library(purrr)
library(ggpubr)


## Functions
# Function to compute Spearman correlation
run_corr <- function(data, celltype_cols) {
  results <- lapply(celltype_cols, function(ct) {
    test <- cor.test(data[[ct]], data$age, method = "spearman")
    data.frame(
      celltype = ct,
      rho = test$estimate,
      pval = test$p.value
    )
  })
  res <- do.call(rbind, results)
  res$fdr <- p.adjust(res$pval, method = "fdr")
  return(res)
}

raw <- import(snakemake@input[["raw"]]) %>% mutate(`CD4/CD8` = `CD4_T`/`CD8_T`)
data <- import(snakemake@input[["data"]]) %>% mutate(`CD4/CD8` = `CD4_T`/`CD8_T`)
celltypes <- colnames(data %>% select(-c(donor_id, age, sex, dataset, ethnicity)))

# Compute correlations
res_both   <- run_corr(raw, celltypes)
res_female <- run_corr(raw[raw$sex == "female", ], celltypes)
res_male   <- run_corr(raw[raw$sex == "male", ], celltypes)


# Merge all tables
final_table <- Reduce(
  function(x, y) merge(x, y, by = "celltype", suffixes = c("", "")),
  list(
    res_both   |> rename(rho_both = rho, p_both = pval, fdr_both = fdr),
    res_female |> rename(rho_female = rho, p_female = pval, fdr_female = fdr),
    res_male   |> rename(rho_male = rho, p_male = pval, fdr_male = fdr)
  )
)

export(final_table, snakemake@output[["table"]])

celltypes <- celltypes[celltypes != "CD4/CD8"]

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

pdf(snakemake@output[["plot3"]], width = 9, height = 2.5)
lapply(celltypes, function(f){

	plots <- lapply(list("both","female","male"), function(g){ 
		df <- raw %>% mutate(agegroup = case_when(
							age < 40 ~ "<40",
							age >= 40 & age < 60 ~ "40-60",
							age >= 60 & age < 80 ~ "60-80",
							age >= 80 ~ ">80")) %>%
				mutate(agegroup = factor(agegroup, levels = c("<40","40-60","60-80",">80")))

		fill_col <- "grey"
		if(g != "both"){
			df <- df %>% filter(sex == g)
			fill_col <- ifelse(g == "female", "#E15566", "#4981BF")
		}

		p1 <- ggplot(df, aes(x = agegroup, y = .data[[f]])) +
			geom_violin(fill = fill_col, width = 0.7) +
			geom_boxplot(fill = "white", width = 0.2) +
			stat_compare_means(paired = FALSE, comparisons = list(c("60-80", ">80")), label = "p", tip.length = 0, bracket.size = 0.7, vjust = -0.4) +
			xlab("") +
			scale_y_continuous(expand = expansion(mult = c(0,.15)))+
			theme_classic(base_size = 15) +
			theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1, size = 15))

		return(p1)
	})
	ggarrange(plotlist = plots, ncol = 3, nrow = 1)
})
dev.off()
