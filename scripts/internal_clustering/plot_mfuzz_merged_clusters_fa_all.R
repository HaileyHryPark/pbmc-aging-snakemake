library(dplyr)
library(tidyverse)
library(ggpubr)
library(ggplot2)
library(colorspace)
library(circlize)
library(rio)

### Functions
fa <- import(snakemake@input[["table"]])
annot <- import(snakemake@input[["annot"]])

fa_top <- fa %>% filter(db == "GO", cluster != "", fa_celltype != "All celltype", type %in% c("both","female","male")) %>% 
	group_by(fa_celltype, cluster, type) %>% 
	arrange(qvalue) %>% slice_head(n = 5) %>% ungroup() %>% pull(term) %>% unique()

fa_res_top <- fa %>% filter(term %in% fa_top,  cluster != "", fa_celltype != "All celltype", type %in% c("both","female","male")) %>% 
  mutate(type = factor(type, levels = c("both","female","male")), 
         cluster = factor(cluster, levels = c("Early\nincrease", "Early\ndecrease", "Continuous\ndecrease", "Early\nfluctuation","Late\nincrease", "Continuous\nincrease")),
         fa_celltype = factor(fa_celltype, levels = c("CD4 T", "CD8 T", "NK", "B", "Mono"))) %>% 
  select(Description, ID, Gender = type, Celltype = fa_celltype, Cluster = cluster, qvalue, gene_name) %>% 
  filter(qvalue < 0.01)

gobp_fa_res_top <- merge(annot, fa_res_top, by.x = "term", by.y = "Description") %>% rename(Category = category, Description = term)

gobp_fa_res_top <- gobp_fa_res_top %>% mutate(Category = factor(Category, levels = c("RNA biosynthetic process", "Translation", "Proteostasis", "Chromosome organization", "OXPHOS/Energy metabolism", "Mitochondria","Immune response", "Immune cell differentiation/activation", "Antigen processing and presentation", "Signaling", "Apoptosis", "Actin fiber organization", "Cellular transport", "Others")))
export(gobp_fa_res_top, snakemake@output[["annotgo"]])

term_order2 <- gobp_fa_res_top %>% 
  group_by(Description, Category) %>% summarise(Count = n()) %>% arrange(Category, desc(Count)) %>% pull(Description)

data_to_plot <- fa %>% filter(Description %in% term_order2, cluster != "", fa_celltype != "All celltype", type %in% c("both","female","male")) %>%
  mutate(type = factor(type, levels = c("both","female","male")), 
         cluster = factor(cluster, levels = c("Early\nincrease", "Early\ndecrease", "Continuous\ndecrease", "Early\nfluctuation","Late\nincrease", "Continuous\nincrease")),
         fa_celltype = factor(fa_celltype, levels = c("CD4 T", "CD8 T", "NK", "B", "Mono"))) %>% 
  select(Description, Gender = type, Celltype = fa_celltype, Cluster = cluster, qvalue) %>% 
  filter(qvalue < 0.01)


celltypes <- unique(data_to_plot$Celltype)

bg <- data.frame(
  xmin = seq(0.5, length(celltypes)-0.5, by = 1),
  xmax = seq(1.5, length(celltypes)+0.5, by = 1),
  ymin = -Inf,
  ymax = Inf,
  fill = rep(c("white", "grey85"), length.out = length(celltypes))
)

pdf(snakemake@output[["go"]], width = 16, height =15)
ggplot(data_to_plot %>%  mutate(Description = factor(Description, levels = rev(term_order2))), 
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
