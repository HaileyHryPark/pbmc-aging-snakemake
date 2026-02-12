library(rio)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(ggplot2)
library(colorspace)
library(ggVennDiagram)
library(svglite)

celltype_colors <- c(
  "CD4 T" = "#D2533B",
  "CD8 T"    = "#E6974D",
  "B"   = "#79629E",
  "NK"= "#73AF68",
  "Mono" = "#5B83BF"
)
celltype_col_df <-  data.frame(celltype = names(celltype_colors), color = celltype_colors)

cluster_col = qualitative_hcl(7, palette = "Set 2")
cluster_level = c("Early\nincrease", "Early\ndecrease", "Continuous\ndecrease", "Early\nfluctuation", "Inverted\nUshape", "Continuous\nincrease", "Late\nincrease")
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

## Bar chart
cor_sub <- bind_rows(cor_sub_f2 %>% mutate(top = "Female top positive correlation"), 
                     cor_sub_f1 %>% mutate(top = "Female top negative correlation"), 
                     cor_sub_m2 %>% mutate(top = "Male top positive correlation"), 
                     cor_sub_m1 %>% mutate(top = "Male top negative correlation")) %>% 
  mutate(top = factor(top, levels = rev(c("Female top positive correlation", "Male top positive correlation", 
                                      "Female top negative correlation", "Male top negative correlation"))))

cor_sub_df <- as.data.frame(table(cor_sub$celltype,cor_sub$final_cluster,cor_sub$top))

cluster_col <- c(cluster_col_df$color)
names(cluster_col) <- cluster_col_df$cluster

p1 <- ggplot(cor_sub_df %>% mutate(Var2 = factor(Var2, levels = cluster_level)), aes(fill = Var2, y = Freq, x = Var3)) + 
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = cluster_col, name = "Cluster") +
  theme_linedraw(base_size = 15)+
  theme(legend.position = "bottom", axis.text.y = element_text(size = 15)) +
  coord_flip()+
  labs(y = "Percentage", x = "")

p2 <- ggplot(cor_sub_df %>% mutate(Var1 = factor(Var1, levels = c("CD4 T", "CD8 T", "B", "NK", "Mono"))), aes(fill = Var1, y = Freq, x = Var3)) + 
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = celltype_colors, name = "Cell type") +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  theme_linedraw(base_size = 15)+
  theme(legend.position = "bottom", axis.text.y = element_blank()) +
  coord_flip()+
  labs(y = "Percentage", x = "")

p <- ggarrange(p1, p2, widths = c(8,5), ncol = 2, nrow = 1)
ggsave(snakemake@output[["bar"]], p, width = 13, height = 4)


