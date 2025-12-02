library(rio)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(ggplot2)
library(colorspace)
library(ggVennDiagram)

celltype_colors <- c(
  "CD4 T" = "#D2533B",
  "CD8 T"    = "#E6974D",
  "NK"= "#73AF68",
  "B"   = "#79629E",
  "Mono" = "#5B83BF"
)
celltype_col_df <-  data.frame(celltype = names(celltype_colors), color = celltype_colors)

cluster_col = qualitative_hcl(6, palette = "Set 2")
cluster_level = c("Early\nincrease", "Early\ndecrease", "Continuous\ndecrease", "Irregular\nchange", "Late\nincrease", "Continuous\nincrease")
cluster_col_df <- data.frame(cluster = cluster_level, color = cluster_col)

## Import cor table
cor <- import(snakemake@input[["annot_cor"]])

cor_sub_f1 <- cor %>% filter(type == "female") %>% filter(rho <= quantile(rho, 0.10, na.rm = T))
cor_sub_f2 <- cor %>% filter(type == "female") %>% filter(rho >= quantile(rho, 0.90, na.rm = T))

cor_sub_m1 <- cor %>% filter(type == "male") %>% filter(rho <= quantile(rho, 0.10, na.rm = T))
cor_sub_m2 <- cor %>% filter(type == "male") %>% filter(rho >= quantile(rho, 0.90, na.rm = T))

## Import fit and spans
mat_f <- import(snakemake@input[["fit_res_f"]]) %>% filter(feature %in% c(cor_sub_f1$feature, cor_sub_f2$feature))
span_f <- import(snakemake@input[["span_res_f"]]) %>% filter(feature %in% c(cor_sub_f1$feature, cor_sub_f2$feature))
span_f <- span_f %>% group_by(feature) %>% slice_min(rmse, n = 1, with_ties = FALSE)

mat_m <- import(snakemake@input[["fit_res_m"]]) %>% filter(feature %in% c(cor_sub_m1$feature, cor_sub_m2$feature))
span_m <- import(snakemake@input[["span_res_m"]]) %>% filter(feature %in% c(cor_sub_m1$feature, cor_sub_m2$feature))
span_m <- span_m %>% group_by(feature) %>% slice_min(rmse, n = 1, with_ties = FALSE)

mat_all <- bind_rows(mat_f %>% mutate(gender = "female"), mat_m %>% mutate(gender = "male"))

mat_f <- mat_f %>%
  left_join(span_f %>% select(feature, span), by = "feature")
mat_m <- mat_m %>%
  left_join(span_m %>% select(feature, span), by = "feature")

## loess
pdf(snakemake@output[["loess"]], width = 5, height = 2.5)
fp <- lapply(list(cor_sub_f1$feature, cor_sub_f2$feature), function(sub){

p <- ggplot(mat_f %>% filter(feature %in% sub), aes(x=age, y=fitted, group=feature)) +
        geom_line(color = "#E15566", size = 0.3, alpha = 0.5) +
        geom_smooth(aes(group = 1), method = "gam", formula = y ~ s(x, k = 10),
              color = "black", size = 2, linetype = "solid") +
        theme_linedraw(base_size = 15)+
        theme(panel.grid = element_blank(), strip.text = element_text(size=17), legend.position = "bottom")+
        ylim(-1.7,1.7)+
        xlim(19,93)+
        labs(x = "Age", y = "Scaled expression")
return(p)

})

mp <- lapply(list(cor_sub_m1$feature, cor_sub_m2$feature), function(sub){

p <- ggplot(mat_m %>% filter(feature %in% sub), aes(x=age, y=fitted, group=feature)) +
        geom_line(color = "#E15566", size = 0.3, alpha = 0.5) +
        geom_smooth(aes(group = 1), method = "gam", formula = y ~ s(x, k = 10),
              color = "black", size = 2, linetype = "solid") +
        theme_linedraw(base_size = 15)+
        theme(panel.grid = element_blank(), strip.text = element_text(size=17), legend.position = "bottom")+
        ylim(-1.7,1.7)+
        xlim(19,93)+
        labs(x = "Age", y = "Scaled expression")
return(p)

})

ggarrange(plotlist = mp, ncol = 2, nrow = 1)
dev.off()

## Venn diagram
pdf(snakemake@output[["venn"]], height = 3.5, width = 4)

venn_pos <- list(Female = cor_sub_f2 %>% pull(feature) %>% unique(), 
                Male = cor_sub_m2 %>% pull(feature) %>% unique())
plot(ggVennDiagram(venn_pos, label_alpha = 0, edge_size = 1.5, set_size = 6, label_size = 5, set_color=c("Female" = "#E15566", "Male" = "#4981BF"), label = "both", fill = "white") + labs(title = "Top percentile (Positive rho)") + theme(legend.position = "none") + scale_fill_gradient(low = "white", high = "white"))

venn_neg <- list(Female = cor_sub_f1 %>% pull(feature) %>% unique(), 
                Male = cor_sub_m1 %>% pull(feature) %>% unique())
plot(ggVennDiagram(venn_neg, label_alpha = 0, edge_size = 1.5, set_size = 6, label_size = 5, set_color=c("Female" = "#E15566", "Male" = "#4981BF"), label = "both", fill = "white") + labs(title = "Top percentile (Negative rho)") + theme(legend.position = "none") + scale_fill_gradient(low = "white", high = "white"))

dev.off()

## Pie chart
pdf(snakemake@output[["pie"]], height = 5, width = 5)

sub_list <- list(cor_sub_f1, cor_sub_f2, cor_sub_m1, cor_sub_m2)
names(sub_list) <- c("Female neg", "Female pos", "Male neg", "Male pos")

lapply(as.list(names(sub_list)), function(s){

  df <- sub_list[[s]]

  bycluster <- as.data.frame(table(df$final_cluster)) %>% arrange(Freq)
  bycluster <- bycluster %>% left_join(cluster_col_df, by = join_by(Var1 == cluster))

  pie(bycluster$Freq, labels = rep("", nrow(bycluster)), col = bycluster$color, main = s)
  pie(bycluster$Freq, labels = bycluster$Var1, col = bycluster$color, main = s)
  mtext(s, side = 3, line = 1, outer = F)

  bycelltype <- as.data.frame(table(df$celltype)) %>% arrange(Freq)
  bycelltype <- bycelltype %>% left_join(celltype_col_df, by = join_by(Var1 == celltype))

  pie(bycelltype$Freq, labels = rep("", nrow(bycelltype)), col = bycelltype$color, main = s)
  pie(bycelltype$Freq, labels = bycelltype$Var1, col = bycelltype$color, main = s)
  mtext(s, side = 3, line = 1, outer = F)
})

dev.off()
