library(dplyr)
library(tidyverse)
library(ggpubr)
library(ggplot2)
library(colorspace)
library(circlize)
library(rio)

celltype_colors <- c(
  "CD4 T" = "#D2533B",
  "CD8 T"    = "#E6974D",
  "NK"= "#73AF68",
  "B"   = "#79629E",
  "Mono" = "#5B83BF"
)

### Functions
fa <- import(snakemake@input[["table"]]) %>% filter(qvalue < 0.05, fa_celltype != "", db == "GO")

pdf(snakemake@output[["plot"]], width = 10, height = 5)
lapply(as.list(unique(fa$cluster)), function(subset){

lapply(as.list(names(celltype_colors)), function(ct){

	sub_fa <- fa %>% filter(cluster == subset, fa_celltype == ct) %>% arrange(qvalue) %>% slice_head(n=5)
	p <- ggplot(sub_fa, aes(y = reorder(Description, -log10(qvalue)))) +
	  geom_col(aes(x = -log10(qvalue)), fill = celltype_colors[ct]) +
	  geom_text(
	    aes(x = 0, label = paste0(Description, "\n")),
	    hjust = 0,
	    nudge_x = max(-log10(sub_fa$qvalue)) / 40,
	    size = 6.5,
	    lineheight = 0.95
	  ) +
	  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
	  labs(
	    x = "-log10(qvalue)",
	    y = "",
	    title = paste(subset, ct)
	  ) +
	  theme_linedraw(base_size = 18) +
	  theme(
	    panel.grid=element_blank(),
	    axis.text.y = element_blank(),  # hide original y-axis text
	    axis.ticks.y = element_blank()
	  )
	plot(p)	
})

})
dev.off()
