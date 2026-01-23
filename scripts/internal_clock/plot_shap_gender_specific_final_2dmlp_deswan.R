library(rio)
library(tidyverse)
library(ggpubr)
library(ggfortify)

gender_cols <- c("female" = "#E96D79", "male" = "#5995CB")

## SHAP analysis

# GAF annoatation
cluster_b <- import(snakemake@input[["clust_b"]]) %>% filter(final_cluster != "")
cluster_f <- import(snakemake@input[["clust_f"]]) %>% filter(final_cluster != "")
cluster_m <- import(snakemake@input[["clust_m"]]) %>% filter(final_cluster != "")

# Preds
b_pred <- import(snakemake@input[["pred_b"]])
f_pred <- import(snakemake@input[["pred_f"]])
m_pred <- import(snakemake@input[["pred_m"]])

# Internal
bi_shap <- import(snakemake@input[["shap_b"]])
bi_pred <- b_pred %>% filter(cohort == "internal")
bi_shap_sub <- bi_pred %>% left_join(bi_shap %>% select(-c(dataset, sex, actual_age, fold, background_n)), by = "donor_id") %>% 
  mutate(age_group = ifelse(actual_age > 60, ">60", ifelse(actual_age < 40, "<40", "40-60")))

fi_shap <- import(snakemake@input[["shap_f"]])
fi_pred <- f_pred %>% filter(cohort == "internal")
fi_shap_sub <- fi_pred %>% left_join(fi_shap %>% select(-c(dataset, sex, actual_age, fold, background_n)), by = "donor_id") %>% 
  mutate(age_diff = predicted_age - actual_age, age_group = ifelse(actual_age > 60, ">60", ifelse(actual_age < 40, "<40", "40-60"))) 

mi_shap <- import(snakemake@input[["shap_m"]])
mi_pred <- m_pred %>% filter(cohort == "internal")
mi_shap_sub <- mi_pred %>% left_join(mi_shap %>% select(-c(dataset, sex, actual_age, fold, background_n)), by = "donor_id") %>% 
  mutate(age_diff = predicted_age - actual_age, age_group = ifelse(actual_age > 60, ">60", ifelse(actual_age < 40, "<40", "40-60")))
  
meta_cols <- c("donor_id", "sample_id", "dataset", "disease", "sex", "actual_age", "predicted_age", "age_group", "age_diff", "c_predicted_age", "c_age_diff", "loess_fit", "fold", "cohort", "type", "age_accel1", "c_age_diff_z","age_accel2")
features_f <- setdiff(colnames(fi_shap_sub), meta_cols)
features_m <- setdiff(colnames(mi_shap_sub), meta_cols)

fi_shap_sub_long <- fi_shap_sub %>% pivot_longer(cols = all_of(features_f), names_to = "feature", values_to = "shap_value") %>% 
  inner_join(cluster_f, by = "feature")
fi_shap_corrected <- fi_shap_sub_long %>% 
  group_by(dataset, feature) %>% 
  mutate(batch_mean = mean(shap_value), shap_corrected = shap_value - batch_mean) %>% 
  ungroup()
mi_shap_sub_long <- mi_shap_sub %>% pivot_longer(cols = all_of(features_m), names_to = "feature", values_to = "shap_value") %>% 
  inner_join(cluster_m, by = "feature")
mi_shap_corrected <- mi_shap_sub_long %>% 
  group_by(dataset, feature) %>% 
  mutate(batch_mean = mean(shap_value), shap_corrected = shap_value - batch_mean) %>% 
  ungroup()

## PCA plots
f_shap_pca <- prcomp(fi_shap_sub %>% select(all_of(features_f)))
m_shap_pca <- prcomp(mi_shap_sub %>% select(all_of(features_m)))

pca_bf_c <- ggarrange(autoplot(f_shap_pca, data = fi_shap_sub, colour = 'dataset', x = 1, y = 2) + theme_test(base_size = 15),
                      autoplot(f_shap_pca, data = fi_shap_sub, colour = 'predicted_age', x = 1, y = 2) + theme_test(base_size = 15),
                      autoplot(m_shap_pca, data = mi_shap_sub, colour = 'dataset', x = 1, y = 2) + theme_test(base_size = 15),
                      autoplot(m_shap_pca, data = mi_shap_sub, colour = 'predicted_age', x = 1, y = 2) + theme_test(base_size = 15),
                      ncol = 2, nrow = 2)
ggsave(snakemake@output[["pcaplot"]], plot = pca_bf_c, width = 10, height = 6)

fi_shap_corrected_wide <- fi_shap_corrected %>% pivot_wider(names_from = feature, id_cols = sample_id, values_from = shap_corrected)
mi_shap_corrected_wide <- mi_shap_corrected %>% pivot_wider(names_from = feature, id_cols = sample_id, values_from = shap_corrected)
identical(fi_shap_corrected_wide$sample_id, fi_shap_sub$sample_id)
f_shap_pca_corrected <- prcomp(fi_shap_corrected_wide %>% select(-sample_id))
m_shap_pca_corrected <- prcomp(mi_shap_corrected_wide %>% select(-sample_id))

pca_af_c <- ggarrange(autoplot(f_shap_pca_corrected, data = fi_shap_sub, colour = 'dataset', x = 1, y = 2) + theme_test(base_size = 15),
                      autoplot(f_shap_pca_corrected, data = fi_shap_sub, colour = 'predicted_age', x = 1, y = 2) + theme_test(base_size = 15),
                      autoplot(m_shap_pca_corrected, data = mi_shap_sub, colour = 'dataset', x = 1, y = 2) + theme_test(base_size = 15),
                      autoplot(m_shap_pca_corrected, data = mi_shap_sub, colour = 'predicted_age', x = 1, y = 2) + theme_test(base_size = 15),
                      ncol = 2, nrow = 2)
ggsave(snakemake@output[["pcaplot_c"]], plot = pca_af_c, width = 10, height = 6)


## SHAP sum by cluster
fi_shap_sub_sum <- fi_shap_sub_long %>% group_by(across(all_of(meta_cols)), final_cluster) %>% 
  summarise(cluster_shap = sum(shap_value, na.rm = T), cluster_shap_abs = abs(sum(shap_value, na.rm = T)), .groups = "drop") %>% 
  mutate(final_cluster = factor(final_cluster, levels = c("Early\nincrease", "Early\ndecrease", "Continuous\ndecrease", "Continuous\nincrease", "Late\nincrease")), 
         age_group = factor(age_group, levels = c("<40", "40-60", ">60")))
mi_shap_sub_sum <- mi_shap_sub_long %>% group_by(across(all_of(meta_cols)), final_cluster) %>% 
  summarise(cluster_shap = sum(shap_value, na.rm = T), cluster_shap_abs = abs(sum(shap_value)), .groups = "drop") %>% 
  mutate(final_cluster = factor(final_cluster, levels = c("Early\nincrease", "Early\ndecrease", "Continuous\ndecrease", "Irregular\nchange")), 
         age_group = factor(age_group, levels = c("<40", "40-60", ">60")))

fi_shap_corrected_sum <- fi_shap_corrected %>% group_by(across(all_of(meta_cols)), final_cluster) %>% 
  summarise(cluster_shap = sum(shap_corrected, na.rm = T), cluster_shap_abs = abs(sum(shap_corrected, na.rm = T)), .groups = "drop") %>% 
  mutate(final_cluster = factor(final_cluster, levels = c("Early\nincrease", "Early\ndecrease", "Continuous\ndecrease", "Continuous\nincrease", "Late\nincrease")), 
         age_group = factor(age_group, levels = c("<40", "40-60", ">60")))
mi_shap_corrected_sum <- mi_shap_corrected %>% group_by(across(all_of(meta_cols)), final_cluster) %>% 
  summarise(cluster_shap = sum(shap_corrected, na.rm = T), cluster_shap_abs = abs(sum(shap_corrected)), .groups = "drop") %>% 
  mutate(final_cluster = factor(final_cluster, levels = c("Early\nincrease", "Early\ndecrease", "Continuous\ndecrease", "Irregular\nchange")), 
         age_group = factor(age_group, levels = c("<40", "40-60", ">60")))

## Correlation with age
fi_corrplots <- ggscatter(fi_shap_corrected_sum, x = "predicted_age", y = "cluster_shap", color = gender_cols["female"], 
                          alpha = 0.7, size = 0.5, cor.coef = T, add = "reg.line", cor.method = "pearson") +
  facet_wrap(~final_cluster, ncol = 5)+
  theme_linedraw(base_size = 11)+
  theme(panel.grid = element_blank())+
  labs(x = "predicted age", y = "summed SHAP")
ggsave(snakemake@output[["corrplot_f"]], fi_corrplots, width = 10, height = 2.5)

mi_corrplots <- ggscatter(mi_shap_corrected_sum, x = "predicted_age", y = "cluster_shap", color = gender_cols["male"], 
                          alpha = 0.7, size = 0.5, cor.coef = T, add = "reg.line", cor.method = "pearson") +
  facet_wrap(~final_cluster, ncol = 4)+
  theme_linedraw(base_size = 11)+
  theme(panel.grid = element_blank())+
  labs(x = "predicted age", y = "summed SHAP")
ggsave(snakemake@output[["corrplot_m"]], mi_corrplots, width = 8, height = 2.5)


## Corrected SHAP by cluster
f_clust_shap <- ggviolin(fi_shap_corrected_sum, x = "final_cluster", y = "cluster_shap_abs", fill = gender_cols["female"], width = 0.6) +
  geom_boxplot(width = 0.1, fill = "white") +
  stat_compare_means(label = "p", comparisons = list(c("Continuous\nincrease","Late\nincrease"),c("Continuous\nincrease","Continuous\ndecrease"),c("Continuous\nincrease","Early\ndecrease"),c("Continuous\nincrease","Early\nincrease")),
                     step.increase = 0.2, tip.length = 0, bracket.size = 0.5, vjust = -0.4) +
  theme_linedraw(base_size = 15) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 15), 
        panel.grid = element_blank(),
        legend.position = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15))) +
  labs(x = "", y = "summed absolute SHAP")
ggsave(snakemake@output[["vlnplot_f"]], f_clust_shap, width = 4, height = 5)

m_clust_shap <- ggviolin(mi_shap_corrected_sum, x = "final_cluster", y = "cluster_shap_abs", fill = gender_cols["male"], width = 0.6) +
  geom_boxplot(width = 0.1, fill = "white") +
  stat_compare_means(label = "p", comparisons = list(c("Early\nincrease","Early\ndecrease"),c("Early\nincrease","Continuous\ndecrease"),c("Early\nincrease","Irregular\nchange")),
                     step.increase = 0.2, tip.length = 0, bracket.size = 0.5, vjust = -0.4) +
  theme_linedraw(base_size = 15) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 15), 
        panel.grid = element_blank(),
        legend.position = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15))) +
  labs(x = "", y = "summed absolute SHAP")
ggsave(snakemake@output[["vlnplot_m"]], m_clust_shap, width = 3.2, height = 5)


## Corrected SHAP by agegroup
f_ag_shap <- ggviolin(fi_shap_corrected_sum, x = "age_group", y = "cluster_shap", fill = gender_cols["female"], width = 0.6, alpha = "age_group") +
  geom_boxplot(width = 0.2, fill = "white") +
  facet_grid(~final_cluster, scales = "free_y") +
  stat_compare_means(label = "p", comparisons = list(c("<40","40-60"),c("40-60",">60"),c("<40",">60")),
                     step.increase = 0.2, tip.length = 0, bracket.size = 0.5, vjust = -0.4) +
  theme_linedraw(base_size = 13) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), 
        panel.grid = element_blank(),
        legend.position = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2))) +
  labs(x = "", y = "summed SHAP")
ggsave(snakemake@output[["vlnplot2_f"]], f_ag_shap, width = 10, height = 3.5)

m_ag_shap <- ggviolin(mi_shap_corrected_sum, x = "age_group", y = "cluster_shap", fill = gender_cols["male"], width = 0.6, alpha = "age_group") +
  geom_boxplot(width = 0.2, fill = "white") +
  facet_grid(~final_cluster, scales = "free_y") +
  stat_compare_means(label = "p", comparisons = list(c("<40","40-60"),c("40-60",">60"),c("<40",">60")),
                     step.increase = 0.2, tip.length = 0, bracket.size = 0.5, vjust = -0.4) +
  theme_linedraw(base_size = 13) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), 
        panel.grid = element_blank(),
        legend.position = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2))) +
  labs(x = "", y = "summed SHAP")
ggsave(snakemake@output[["vlnplot2_m"]], m_ag_shap, width = 8, height = 3.5)

## Corrected SHAP by agegroup and acceleration
f_accel <- ggviolin(fi_shap_corrected_sum %>% mutate(age_accel1 = factor(age_accel1, levels = c("Decelerated", "Accelerated"))), 
         x = "age_accel1", y = "cluster_shap", fill = "age_accel1", width = 0.6) +
  geom_boxplot(width = 0.15, fill = "white") +
  facet_grid(age_group~final_cluster, scales = "free", switch = "y") +
  stat_compare_means(label = "p", comparisons = list(c("Accelerated","Decelerated")), tip.length = 0, bracket.size = 0.5, vjust = -0.4) +
  scale_fill_manual(values = c("Decelerated" = "#E96D7930", "Accelerated" = "darkred"), name = "Age acceleration", labels = c("Slower aging", "Accelerated aging")) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2)), position = "right") +
  theme_linedraw(base_size = 15) +
  theme(axis.text.x = element_blank(), legend.position = "bottom", 
        panel.grid.major = element_blank()) +
  labs(x = "", y = "summed SHAP")
ggsave(snakemake@output[["vlnplot3_f"]], f_accel, width = 8, height = 6)

m_accel <- ggviolin(mi_shap_corrected_sum %>% mutate(age_accel1 = factor(age_accel1, levels = c("Decelerated", "Accelerated"))), 
                    x = "age_accel1", y = "cluster_shap", fill = "age_accel1", width = 0.6) +
  geom_boxplot(width = 0.15, fill = "white") +
  facet_grid(age_group~final_cluster, scales = "free", switch = "y") +
  stat_compare_means(label = "p", comparisons = list(c("Accelerated","Decelerated")), tip.length = 0, bracket.size = 0.5, vjust = -0.4) +
  scale_fill_manual(values = c("Decelerated" = "#5995CB30", "Accelerated" = "navy"), name = "Age acceleration", labels = c("Slower aging", "Accelerated aging")) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2)), position = "right") +
  theme_linedraw(base_size = 15) +
  theme(axis.text.x = element_blank(), legend.position = "bottom", 
        panel.grid.major = element_blank()) +
  labs(x = "", y = "summed SHAP")
ggsave(snakemake@output[["vlnplot3_m"]], m_accel, width = 6.4, height = 6)
