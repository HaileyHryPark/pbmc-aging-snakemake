library(dplyr)
library(tidyverse)
library(ggpubr)
library(ggplot2)
library(colorspace)
library(circlize)
library(rio)

cluster_level = c("Early\nincrease", "Early\ndecrease", "Continuous\ndecrease", "Irregular\nchange", "Late\nincrease", "Continuous\nincrease", "Inverted\nU-shape")
celltype_level = c("CD4 T", "CD8 T", "NK", "B", "Mono")

### Functions
PlotFARes <-function(df, database){

fa_top <- df %>% filter(db == database, qvalue < 0.05) %>% group_by(cluster, gender) %>% 
  arrange(qvalue) %>% slice_head(n = 10) %>% ungroup() %>% pull(term) %>% unique()
print(head(fa_top))

df_top <- df %>% filter(term %in% fa_top, fa_celltype != "All celltype") %>% 
  mutate(gender = factor(gender, levels = c("both","female","male")), 
         cluster = factor(cluster, levels = cluster_level),
         celltype = factor(fa_celltype, levels = celltype_level)) %>% 
  select(Description, Gender = gender, Celltype = celltype, Cluster = cluster, qvalue) %>% 
  filter(qvalue < 0.05)
print(head(df_top))

term_order <- df_top %>% group_by(Description) %>% summarise(Count = n()) %>% arrange(desc(Count)) %>% pull(Description)

celltypes <- unique(df_top$Celltype)
print(celltypes)

# Create alternating background rectangles
bg <- data.frame(
  xmin = seq(0.5, length(celltypes)-0.5, by = 1),
  xmax = seq(1.5, length(celltypes)+0.5, by = 1),
  ymin = -Inf,
  ymax = Inf,
  fill = rep(c("white", "grey85"), length.out = length(celltypes))
)
print(head(bg))

p <- ggplot(df_top %>% mutate(Description = factor(Description, levels = rev(term_order))), 
       aes(x = Celltype, y = Description, color = Gender, group = Gender)) +
  #geom_rect(data = bg, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
  #          fill = bg$fill, color = NA, inherit.aes = FALSE, alpha = 0.5) +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  facet_grid(~Cluster)+
  scale_color_manual(values = c("both" = "grey", "female" = "#E15566", "male" = "#4981BF")) +
  theme_linedraw()+
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    panel.grid.major = element_line(color = "grey"),
    panel.grid.minor = element_blank()
  ) +
  labs(x = "", y = "Pathways", color = "Gender")

plot(p)

}

both <- import(snakemake@input[["both1"]]) %>% mutate(gender = "both")
both_nr <- import(snakemake@input[["both2"]]) %>% mutate(gender = "both")
female <- import(snakemake@input[["female1"]]) %>% mutate(gender = "female")
female_nr <- import(snakemake@input[["female2"]]) %>% mutate(gender = "female")
male <- import(snakemake@input[["male1"]]) %>% mutate(gender = "male")
male_nr <- import(snakemake@input[["male2"]]) %>% mutate(gender = "male")

fa_res <- bind_rows(both, female, male)
fa_res_nr <- bind_rows(both_nr, female_nr, male_nr)

pdf(snakemake@output[["plot1"]], width = 18, height = 10)

PlotFARes(fa_res, "GO")
#PlotFARes(fa_res, "KEGG")
PlotFARes(fa_res, "WP")
PlotFARes(fa_res, "Reactome")

dev.off()

pdf(snakemake@output[["plot2"]], width = 27, height = 10)

PlotFARes(fa_res_nr, "GO")
#PlotFARes(fa_res_nr, "KEGG")
PlotFARes(fa_res_nr, "WP")
PlotFARes(fa_res_nr, "Reactome")

dev.off()

