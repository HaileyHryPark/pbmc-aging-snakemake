library(rio)
library(tidyverse)
library(ggpubr)
library(svglite)
library(Metrics)

cvb <- import(snakemake@input[["cvb"]]) %>% filter(!dataset %in% c("immage","soundlife")) %>% mutate(actual_age = as.integer(actual_age)) %>% filter(!is.na(actual_age)) %>% group_by(fold) %>% summarise(RMSE = rmse(actual_age, predicted_age), MAE = mae(actual_age, predicted_age), .groups = "drop") %>% select(fold, RMSE, MAE) %>% mutate(sex = "Both-sex", group = "5-fold CV", fold = as.character(fold))
cvf <- import(snakemake@input[["cvf"]]) %>% filter(!dataset %in% c("immage","soundlife")) %>% mutate(actual_age = as.integer(actual_age)) %>% filter(!is.na(actual_age)) %>% group_by(fold) %>% summarise(RMSE = rmse(actual_age, predicted_age), MAE = mae(actual_age, predicted_age), .groups = "drop") %>% select(fold, RMSE, MAE) %>% mutate(sex = "Female", group = "5-fold CV", fold = as.character(fold))
cvm <- import(snakemake@input[["cvm"]]) %>% filter(!dataset %in% c("immage","soundlife")) %>% mutate(actual_age = as.integer(actual_age)) %>% filter(!is.na(actual_age)) %>% group_by(fold) %>% summarise(RMSE = rmse(actual_age, predicted_age), MAE = mae(actual_age, predicted_age), .groups = "drop") %>% select(fold, RMSE, MAE) %>% mutate(sex = "Male", group = "5-fold CV", fold = as.character(fold))
print(head(cvm))

lodob <- import(snakemake@input[["lodob"]]) %>% select(fold = held_out_train_domain, RMSE = rmse_external, MAE = mae_external) %>% mutate(sex = "Both-sex", group = "Leave-One-Dataset-Out")
lodof <- import(snakemake@input[["lodof"]]) %>% select(fold = held_out_train_domain, RMSE = rmse_external, MAE = mae_external) %>% mutate(sex = "Female", group = "Leave-One-Dataset-Out")
lodom <- import(snakemake@input[["lodom"]]) %>% select(fold = held_out_train_domain, RMSE = rmse_external, MAE = mae_external) %>% mutate(sex = "Male", group = "Leave-One-Dataset-Out")
print(head(lodom))

data <- bind_rows(cvb, cvf, cvm, lodob, lodof, lodom)
export(data, snakemake@output[["res"]])

p <- ggplot(data, aes(x = group, y = RMSE)) +
        geom_boxplot(aes(color = group), width = 0.6) +
        scale_color_manual(values = c("5-fold CV" = "grey", "Leave-One-Dataset-Out" = "red")) +
        geom_point(position = position_jitter(width = 0.2), size = 1, alpha = 0.8, color = "black") +
	facet_wrap(~sex) +
        ylim(12, 25)+
        xlab("")+
        theme_linedraw(base_size = 15)+
        theme(panel.grid = element_blank(), legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(snakemake@output[["plot"]], p, width = 5, height = 3.5)

