library(rio)
library(limma)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggrepel)


# Read data
data <- import(snakemake@input[["mvalue"]])
metadata <- import(snakemake@input[["metadata"]]) %>% filter(!is.na(Age), Gender != "")
print(length(unique(metadata$Basename)))
print(length(unique(metadata$ID)))
print(metadata %>% filter(is.na(Age)))
print(metadata %>% filter(is.na(AgeGroup)))
print(metadata %>% filter(Gender == ""))

mat <- data %>% column_to_rownames("probe_prefix") %>% as.matrix()

print(head(mat))
print(dim(mat))
print(dim(metadata))

# Subset metadata to match NPX matrix
mat <- mat[, metadata$Basename]
print(dim(mat))

# Design matrix
metadata$AgeGroup <- factor(metadata$AgeGroup, levels = c("<40","40-60",">60"))
print(table(metadata$AgeGroup))
design1 <- model.matrix(~ Age * Gender, data = metadata)
design2 <- model.matrix(~ AgeGroup * Gender, data = metadata)

dlist <- list(design1, design2)
names(dlist) <- c("age", "agegroup")

res <- lapply(as.list(names(dlist)), function(d){

	design <- dlist[[d]]
	print(design)
	print(dim(design))
	
	# Fit linear model
	fit <- lmFit(mat, design)
	
	fit <- eBayes(fit)

	res2 <- lapply(as.list(colnames(design)[2:4]), function(c){
		r <- topTable(fit, coef = c, number = Inf) %>% as.data.frame() %>% 
			mutate(design = d, coef = c) %>%
			rownames_to_column("probeID")
		return(r)
	})
	
	return(bind_rows(res2))

})

export(res[[1]], snakemake@output[["res1"]])
export(res[[2]], snakemake@output[["res2"]])

