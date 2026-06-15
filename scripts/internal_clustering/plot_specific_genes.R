library(rio)
library(ggpubr)
library(tidyverse)
library(svglite)

celltype_cols <- c("CD4 T" = "#D2533B", "CD8 T" = "#E6974D", "B" = "#79629E", "NK" = "#73AF68", "Mono" = "#5B83BF")

## Functions
plotLimmaScatterFig <- function(limmab, spanf, spanm, f){
  spanff <- spanf %>% filter(feature == f) %>% slice_min(rmse) %>% pull(span)
  spanmf <- spanm %>% filter(feature == f) %>% slice_min(rmse) %>% pull(span)
  
  d1 <- limmab %>% select(all_of(f), age, sex, dataset)
  
  p1 <- ggplot(d1, aes(x = age, y = .data[[f]])) + 
    geom_point(aes(color = sex), alpha = 0.2, size = 0.2) +
    geom_smooth(aes(color = sex), method = "loess") +
    scale_color_manual(values = c("female" = "#E15566", "male" = "#4981BF")) +
    theme_linedraw(base_size = 15) + theme(legend.position = "none", panel.grid = element_blank())
  
  return(p1)
}

## Main
limma <- import(snakemake@input[["limma"]])
limma_f <- import(snakemake@input[["limma_f"]])
span_f <- import(snakemake@input[["span_f"]])
span_m <- import(snakemake@input[["span_m"]])

## Continuous increase, late increase, inverted U-shape genes
p1 <- ggplot(limma %>% select(`CD8 T.HLA-C`, age, sex, dataset), aes(x = age, y = `CD8 T.HLA-C`)) +
  geom_point(aes(color = sex), alpha = 0.2, size = 0.2) +
  geom_smooth(aes(color = sex), method = "loess") +
  scale_color_manual(values = c("female" = "#E15566", "male" = "#4981BF"), labels = c("Female", "Male"), name = "Gender") +
  theme_linedraw(base_size = 15) + theme(legend.position = "none", panel.grid = element_blank())

ps <- lapply(list("CD8 T.SIGIRR", "CD8 T.MYO1G", "CD8 T.CD74", "CD8 T.ITGB2", "CD8 T.CTSS", "CD4 T.ATP6AP2", "CD4 T.AP1S2", "CD4 T.ANP32E", "CD4 T.ARRDC3", "CD4 T.CAPZA1", "CD4 T.CLN5", "CD4 T.TRAM1", "NK.TPT1", "CD4 T.TPT1", "NK.RACK1", "CD4 T.RACK1", "CD4 T.LAPTM5"), function(f){
  return(plotLimmaScatterFig(limma, span_f, span_m, f))
})

allp <- ggarrange(plotlist = c(p1, ps), nrow = 4, ncol = 5)
ggsave(snakemake@output[["f_genes"]], allp, width = 16, height = 12)


## Inverted Ushape genes main
p2 <- ggplot(limma %>% select(`B.HLA-DMA`, age, sex, dataset), aes(x = age, y = `B.HLA-DMA`)) +
  geom_smooth(aes(color = sex), method = "loess") +
  scale_color_manual(values = c("female" = "#E15566", "male" = "#4981BF"), labels = c("Female", "Male"), name = "Gender") +
  ylim(-2,2)+
  theme_linedraw(base_size = 15) + theme(legend.position = "none", panel.grid = element_blank())

p3 <- ggplot(limma %>% select(`B.HLA-DMB`, age, sex, dataset), aes(x = age, y = `B.HLA-DMB`)) +
  geom_smooth(aes(color = sex), method = "loess") +
  scale_color_manual(values = c("female" = "#E15566", "male" = "#4981BF"), labels = c("Female", "Male"), name = "Gender") +
  ylim(-2,2)+
  theme_linedraw(base_size = 15) + theme(legend.position = "none", panel.grid = element_blank())

p4 <- ggplot(limma %>% select(`B.HLA-DRA`, age, sex, dataset), aes(x = age, y = `B.HLA-DRA`)) +
  geom_smooth(aes(color = sex), method = "loess") +
  scale_color_manual(values = c("female" = "#E15566", "male" = "#4981BF"), labels = c("Female", "Male"), name = "Gender") +
  ylim(-2,2)+
  theme_linedraw(base_size = 15) + theme(legend.position = "none", panel.grid = element_blank())

p5 <- ggplot(limma_f %>% select(`B.HLA-DRB1`, age, sex, dataset), aes(x = age, y = `B.HLA-DRB1`)) +
  geom_smooth(aes(color = sex), method = "loess") +
  scale_color_manual(values = c("female" = "#E15566", "male" = "#4981BF"), labels = c("Female", "Male"), name = "Gender") +
  ylim(-2,2)+
  theme_linedraw(base_size = 15) + theme(legend.position = "none", panel.grid = element_blank())

allp2 <- ggarrange(p2, p3, p4, p5, nrow = 2, ncol = 2)
ggsave(snakemake@output[["fiu_genes"]], allp2, width = 4, height = 3.9)

## Late increase iron genes
iron_df <- limma %>% pivot_longer(cols = ends_with("FTL"), names_to = "celltype", values_to = "FTL_expr") %>% 
  mutate(celltype = gsub("\\.FTL$", "", celltype), sex = factor(sex, levels = c("female", "male"), labels = c("Female", "Male"))) 
iron_df2 <- limma %>% pivot_longer(cols = ends_with("FTH1"), names_to = "celltype", values_to = "FTH1_expr") %>% 
  mutate(celltype = gsub("\\.FTH1$", "", celltype), sex = factor(sex, levels = c("female", "male"), labels = c("Female", "Male"))) 

ftl <- ggplot(iron_df, aes(x = age, y = FTL_expr, color = celltype)) +
  facet_wrap(~sex) +
  geom_point(color = "white", alpha = 0) +
  geom_smooth(method = "loess", se = FALSE, span = 0.8, linewidth = 1.2) +
  theme_linedraw(base_size = 15) +
  theme(panel.grid.major = element_blank(),legend.position = "bottom") + 
  scale_color_manual(values = celltype_cols) +
  ylim(-2,2)+
  labs(
    x = "Age",
    y = "Scaled FTL Expression",
    color = "Cell type"
  )

fth1 <- ggplot(iron_df2, aes(x = age, y = FTH1_expr, color = celltype)) +
  facet_wrap(~sex) +
  geom_point(color = "white", alpha = 0) +
  geom_smooth(method = "loess", se = FALSE, span = 0.8, linewidth = 1.2) +
  theme_linedraw(base_size = 15) +
  theme(panel.grid.major = element_blank(),legend.position = "bottom") + 
  scale_color_manual(values = celltype_cols) +
  ylim(-2,2)+
  labs(
    x = "Age",
    y = "Scaled FTH1 Expression",
    color = "Cell type"
  )

ggsave(snakemake@output[["ftl"]], ftl, width = 5, height = 4)
ggsave(snakemake@output[["fth1"]], fth1, width = 5, height = 4)

## FCIMEI genes related to dnam
ps2 <- lapply(list("CD4 T.CLIC1", "CD8 T.CLIC1"), function(f){
	return(plotLimmaScatterFig(limma, span_f, span_m, f))
})
allp2 <- ggarrange(plotlist = ps2, nrow = 2, ncol = 1)
ggsave(snakemake@output[["fcimei_genes"]], allp2, width = 3, height = 5.8)

