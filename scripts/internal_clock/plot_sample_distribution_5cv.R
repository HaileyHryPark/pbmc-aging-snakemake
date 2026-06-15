library(rio)
library(tidyr)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(svglite)

## Functions

sexcolor <- c("female" = "#E96D79", "male" = "#5995CB")
dataset_cols <- c("onek1k" = "#1b9e77", "aida" = "#d95f02", "perez" = "#7570b3", "marina" = "#e7298a")

plotFoldDistributionBoth <- function(df){

p <- lapply(as.list(unique(df$outer_fold)), function(f){

	m <- df %>% filter(outer_fold == f, split == "train")

	p1 <- ggplot(m, aes(x = age)) +
        	geom_histogram(aes(fill = sex), color = "white", binwidth = 5, boundary = 50) +
        	scale_fill_manual(values = sexcolor, name = "Sex") +
        	xlab("Age")+
        	ylab("Number of samples") +
		ggtitle(f) +
		ylim(0,150) +
		xlim(18,100) +
        	theme_classic()+
		theme(panel.grid = element_blank(), strip.text = element_text(size=17))
	p2 <- ggplot(m, aes(x = age)) +
        	geom_histogram(aes(fill = dataset), color = "white", binwidth = 5, boundary = 50) +
        	scale_fill_manual(values = dataset_cols, name = "Dataset") +
        	xlab("Age")+
        	ylab("Number of samples") +
		ggtitle(f) +
		ylim(0,150) +
		xlim(18,100) +
        	theme_classic()+
		theme(panel.grid = element_blank(), strip.text = element_text(size=17))
	return(ggarrange(p1, p2, ncol = 1, nrow = 2))
})

return(ggarrange(plotlist=p, ncol = length(unique(df$outer_fold)), nrow = 1))

}

plotFoldDistributionSex <- function(df){

p <- lapply(as.list(unique(df$fold)), function(f){

	m <- df %>% filter(fold != f)

	p1 <- ggplot(m, aes(x = actual_age)) +
        	geom_histogram(aes(fill = sex), color = "white", binwidth = 5, boundary = 50) +
        	scale_fill_manual(values = sexcolor, name = "Sex") +
        	xlab("Age")+
        	ylab("Number of samples") +
		ggtitle(f) +
		ylim(0,150) +
		xlim(18,100) +
        	theme_classic()+
		theme(panel.grid = element_blank(), strip.text = element_text(size=17))
	p2 <- ggplot(m, aes(x = actual_age)) +
        	geom_histogram(aes(fill = dataset), color = "white", binwidth = 5, boundary = 50) +
        	scale_fill_manual(values = dataset_cols, name = "Dataset") +
        	xlab("Age")+
        	ylab("Number of samples") +
		ggtitle(f) +
		ylim(0,150) +
		xlim(18,100) +
        	theme_classic()+
		theme(panel.grid = element_blank(), strip.text = element_text(size=17))
	return(ggarrange(p1, p2, ncol = 1, nrow = 2))
})

return(ggarrange(plotlist=p, ncol = length(unique(df$fold)), nrow = 1))

}

metab <- import(snakemake@input[["datab"]])
metaf <- import(snakemake@input[["dataf"]])
metam <- import(snakemake@input[["datam"]])

ggsave(snakemake@output[["plotb"]], plotFoldDistributionBoth(metab), width = 17.5, height = 4.5)
ggsave(snakemake@output[["plotf"]], plotFoldDistributionSex(metaf), width = 17.5, height = 4.5)
ggsave(snakemake@output[["plotm"]], plotFoldDistributionSex(metam), width = 17.5, height = 4.5)


sumb <- metab %>% group_by(outer_fold, split) %>% summarize(n = n(), mean_age = mean(age)) %>% mutate(sex = "both")
print(head(sumb))
sumf_train <- lapply(sort(unique(metaf$fold)), function(f) {

  test_ids <- metaf %>%
    filter(fold == f) %>%
    pull(donor_id)

  metaf %>%
    filter(!donor_id %in% test_ids) %>%
    summarize(
      fold = f,
      n = n(),
      mean_age = mean(actual_age),
    )

}) %>%
  bind_rows()
summ_train <- lapply(sort(unique(metam$fold)), function(f) {

  test_ids <- metam %>%
    filter(fold == f) %>%
    pull(donor_id)

  metam %>%
    filter(!donor_id %in% test_ids) %>%
    summarize(
      fold = f,
      n = n(),
      mean_age = mean(actual_age),
    )

}) %>%
  bind_rows()

sumf <- bind_rows(metaf %>% group_by(fold) %>% summarize(n = n(), mean_age = mean(actual_age)) %>% mutate(split = "test"), 
		sumf_train %>% mutate(split = "train")) %>% select(outer_fold = fold, split, n, mean_age) %>% mutate(sex = "female")
summ <- bind_rows(metam %>% group_by(fold) %>% summarize(n = n(), mean_age = mean(actual_age)) %>% mutate(split = "test"), 
		summ_train %>% mutate(split = "train")) %>% select(outer_fold = fold, split, n, mean_age) %>% mutate(sex = "male")


export(bind_rows(sumb, sumf, summ), snakemake@output[["res"]])
