library(dplyr)
library(tidyverse)
library(ggpubr)
library(ggplot2)
library(colorspace)
library(circlize)
library(rio)

### Functions
fa <- import(snakemake@input[["table"]])

fa_top <- fa %>% filter(cluster != "", fa_celltype != "All celltype", type %in% c("both","female","male")) %>% 
	group_by(fa_celltype, cluster, type) %>% 
	arrange(qvalue) %>% slice_head(n = 5) %>% ungroup() %>% pull(term) %>% unique()

fa_res_top <- fa %>% filter(term %in% fa_top,  cluster != "", fa_celltype != "All celltype", type %in% c("both","female","male")) %>% 
  mutate(type = factor(type, levels = c("both","female","male")), 
         cluster = factor(cluster, levels = c("Early\nincrease", "Early\ndecrease", "Continuous\ndecrease", "Early\nfluctuation","Late\nincrease", "Continuous\nincrease")),
         fa_celltype = factor(fa_celltype, levels = c("CD4 T", "CD8 T", "NK", "B", "Mono"))) %>% 
  select(Description, Gender = type, Celltype = fa_celltype, Cluster = cluster, qvalue, geneID) %>% 
  filter(qvalue < 0.01)

term_order <- fa_res_top %>% 
  group_by(Description) %>% summarise(Count = n()) %>% arrange(desc(Count)) %>% pull(Description)

celltypes <- unique(fa_res_top$Celltype)

bg <- data.frame(
  xmin = seq(0.5, length(celltypes)-0.5, by = 1),
  xmax = seq(1.5, length(celltypes)+0.5, by = 1),
  ymin = -Inf,
  ymax = Inf,
  fill = rep(c("white", "grey85"), length.out = length(celltypes))
)

pdf(snakemake@output[["plot"]], width = 12, height =5)
ggplot(fa_res_top %>%  mutate(Description = factor(Description, levels = rev(term_order))), 
       aes(x = Celltype, y = Description, color = Gender, group = Gender)) +
  geom_rect(data = bg, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
            inherit.aes = FALSE, alpha = 0.5) +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  facet_grid(~Cluster)+
  scale_fill_identity() + 
  scale_color_manual(values = c("both" = "grey","female" = "#E15566", "male" = "#4981BF")) +
  theme_linedraw()+
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    panel.grid.major = element_line(color = "grey30"),
    panel.grid.minor = element_blank()
  ) +
  labs(x = "", y = "Pathways", color = "Gender")
dev.off()
