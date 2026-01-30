library(rio)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(ggplot2)
library(ggalluvial)
library(colorspace)
library(ggVennDiagram)


cluster_col = qualitative_hcl(6, palette = "Set 2")
celltype_colors <- c(
  "CD4 T" = "#D2533B",
  "CD8 T"    = "#E6974D",
  "NK"= "#73AF68",
  "B"   = "#79629E",
  "Mono" = "#5B83BF"
)
cluster_level = c("Early\nincrease", "Early\ndecrease", "Continuous\ndecrease", "Early\nfluctuation", "Late\nincrease", "Continuous\nincrease")
cluster_col_df <- data.frame(cluster = cluster_level, color = cluster_col)

both <- import(snakemake@input[["both"]]) %>% mutate(gender = "both")
female <- import(snakemake@input[["female"]]) %>% mutate(gender = "female")
male <- import(snakemake@input[["male"]]) %>% mutate(gender = "male")

df <- bind_rows(both, female, male)
print(head(df))

df <- df %>% filter(!is.na(final_cluster)) %>% mutate(final_cluster = factor(final_cluster, levels = cluster_level), gender = factor(gender, levels = c("male","both","female")))

freqdf <- df %>% group_by(gender, feature, final_cluster) %>% summarise(Freq = n(), .groups = "drop")
print(head(freqdf))

pdf(snakemake@output[["flowplot"]], height = 7, width = 5)

## Male -> Both -> Female Flow plot
ggplot(freqdf, aes(x = gender, stratum = final_cluster, alluvium = feature, y = Freq, fill = final_cluster)) + 
	scale_x_discrete(expand = c(.1, .1)) +
	geom_flow() +
	geom_stratum(alpha = .5) +
	scale_fill_manual(values = cluster_col)+
	theme_classic(base_size = 15)+
	theme(legend.position = "right", axis.line = element_blank()) +
	xlab("") + ylab("")

## Male -> Female Flow plot
ggplot(freqdf %>% filter(gender != "both"), aes(x = gender, stratum = final_cluster, alluvium = feature, y = Freq, fill = final_cluster)) + 
	scale_x_discrete(expand = c(.1, .1)) +
	geom_flow() +
	geom_stratum(alpha = .5) +
	scale_fill_manual(values = cluster_col)+
	theme_classic(base_size = 15)+
	theme(legend.position = "right", axis.line = element_blank()) +
	xlab("") + ylab("")

dev.off()

## Venn diagram
pdf(snakemake@output[["venn"]], height = 7, width = 7)
venns <- lapply(as.list(unique(df$final_cluster)), function(cl){

cldf <- df %>% filter(final_cluster == cl)
  
venn_list <- list(Both = cldf %>% filter(gender == "both") %>% pull(feature) %>% unique(), 
		Female = cldf %>% filter(gender == "female") %>% pull(feature) %>% unique(),
		Male = cldf %>% filter(gender == "male") %>% pull(feature) %>% unique())

plot(ggVennDiagram(venn_list, label_alpha = 0, edge_size = 0.5, set_size = 6) + ggtitle(cl))

})
dev.off()

## Pie chart
pdf(snakemake@output[["pie"]], height = 5, width = 5)

lapply(list("both", "female", "male"), function(g){
lapply(as.list(unique(df$celltype)), function(ct){

ctdf <- df %>% filter(celltype == ct, gender == g)

  bycluster <- as.data.frame(table(ctdf$final_cluster)) %>% arrange(Freq)
  bycluster <- bycluster %>% left_join(cluster_col_df, by = join_by(Var1 == cluster))

  pie(bycluster$Freq, labels = rep("", nrow(bycluster)), col = bycluster$color, main = paste(g, ct))
  pie(bycluster$Freq, labels = bycluster$Var1, col = bycluster$color, main = paste(g, ct))
  mtext(ct, side = 3, line = 1, outer = F)

})
})

dev.off()
