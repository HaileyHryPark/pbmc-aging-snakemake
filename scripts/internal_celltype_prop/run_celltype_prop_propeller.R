library(rio)
library(dplyr)
library(tidyverse)
library(purrr)
library(ggpubr)
library(speckle)
library(limma)

dataset_cols <- c("onek1k" = "#1b9e77", "aida" = "#d95f02", "perez" = "#7570b3", "marina" = "#e7298a")

## Functions
run_propeller <- function(df, cluster_col){

donor_meta <- df %>% select(donor_id, dataset, age, agegroup, sex) %>% distinct()

prop.logit <- getTransformedProps(clusters = df[[cluster_col]],
				sample = df$donor_id,
				transform = "logit")
print(dim(prop.logit$TransformedProps))

donor_meta <- donor_meta[match(colnames(prop.logit$Counts), donor_meta$donor_id),]

## Anova
design.anova <- model.matrix(~agegroup+dataset, data = donor_meta)
print(head(design.anova))
print(dim(design.anova))

fit.anova <- propeller.anova(prop.list = prop.logit, 
				design = design.anova, 
				coef = 4:6, 
				robust = T, trend = F, sort = T)
fit.anova$celltype <- rownames(fit.anova)

## modeling age as continuous variable
design.age <- model.matrix(~age+dataset, data = donor_meta)
print(head(design.age))

fit.age <- lmFit(prop.logit$TransformedProps, design.age)
fit.age <- eBayes(fit.age, robust = TRUE)
fit.age.tbl <- topTable(fit.age, coef = "age", number = Inf)
fit.age.tbl$celltype <- rownames(fit.age.tbl)


prop_df <- as.data.frame(t(prop.logit$Proportions))
print(head(prop_df))

prop_df$donor_id <- prop_df$sample 
prop_df$celltype <- prop_df$clusters 

prop_df <- prop_df %>% left_join(donor_meta, by = "donor_id")

#prop_long <- prop_df %>%
#    pivot_longer(
#      cols = -c(donor_id, sample, dataset, age, agegroup, sex),
#      names_to = "celltype",
#      values_to = "proportion"
#)
print(head(prop_df))

vlnplots <- lapply(unique(prop_df$celltype), function(ct){

p <- ggplot(prop_df %>% filter(celltype == ct), aes(x = dataset, y = Freq, fill = dataset)) +
	geom_violin(trim = F, alpha = 0.7) +
	geom_boxplot(fill = "white", width = 0.15) +
	scale_fill_manual(values = dataset_cols) +
	theme_classic(base_size = 14)+
	labs(title = ct, x = "Dataset", y = "Cell type proportion") +
	ylim(0, 1) +
	theme(plot.title = element_text(hjust=0.5), legend.position = "none")

return(p)
})

vlnplots_combined <- ggarrange(plotlist = vlnplots, ncol = 4, nrow = ceiling(length(unique(prop_df$celltype))/4))

return(list(anova = fit.anova, age = fit.age.tbl, plots = vlnplots_combined))

}

### Main

data <- import(snakemake@input[["data"]]) 

if(snakemake@params[["gender"]] != "both"){
	data <- data %>% filter(sex == snakemake@params[["gender"]])
}

data <- data %>% mutate(agegroup = ifelse(age < 40, "Young", ifelse(age > 60, "Old", "Middleaged")))
data$agegroup <- factor(data$agegroup, levels = c("Young", "Middleaged", "Old"))

l1 <- run_propeller(df = data, cluster_col = "predicted.celltype.l1")
l2 <- run_propeller(df = data, cluster_col = "predicted.celltype.l2")

export(l1$anova, snakemake@output[["l1res_anova"]])
export(l1$age, snakemake@output[["l1res_age"]])
ggsave(snakemake@output[["l1_vlnplot"]], l1$plots, width = 16, height = 8)

export(l2$anova, snakemake@output[["l2res_anova"]])
export(l2$age, snakemake@output[["l2res_age"]])
ggsave(snakemake@output[["l2_vlnplot"]], l2$plots, width = 16, height = 24)

