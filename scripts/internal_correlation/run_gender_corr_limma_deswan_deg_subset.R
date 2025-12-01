library(rio)
library(dplyr)
library(tidyverse)
library(purrr)

## Functions


limma <- import(snakemake@input[["limma"]])
limma_b <- limma %>% filter(age < 60) 
limma_a <- limma %>% filter(age >= 60) 

limma_list <- list(limma, limma_b, limma_a)
names(limma_list) <- c("All", "Before", "After")

features <- colnames(limma %>% select(-c(rowname,age,sex,dataset,ethnicity)))

cor_results <- lapply(features, function(f) {
res <- lapply(as.list(names(limma_list)), function(l){
  x <- limma_list[[l]][[f]]
  y <- limma_list[[l]]$age
  
  # run correlation
  ct1 <- suppressWarnings(cor.test(x, y, method = "pearson"))
  ct2 <- suppressWarnings(cor.test(x, y, method = "spearman"))
  
  df <- data.frame(
    feature   = f,
    r         = ct1$estimate,
    p_value_p = ct1$p.value,
    rho       = ct2$estimate,
    p_value_s = ct2$p.value
  )

  if(l == "Before"){
    colnames(df) <- paste(colnames(df), "b", sep = ".")
    df <- df %>% rename(feature = feature.b)
  }else if(l == "After"){
    colnames(df) <- paste(colnames(df), "a", sep = ".")
    df <- df %>% rename(feature = feature.a)
  }
  return(df)
})
return(Reduce(function(x,y) merge(x, y, by = "feature", all = TRUE), res))
})

# Bind results
cor_limma <- bind_rows(cor_results)

# Adjust p-values
cor_limma$p_adj_p <- p.adjust(cor_limma$p_value_p, method = "BH")
cor_limma$p_adj_s <- p.adjust(cor_limma$p_value_s, method = "BH")

# Final output
export(cor_limma, snakemake@output[["cor"]])
