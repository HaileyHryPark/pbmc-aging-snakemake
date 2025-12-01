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
