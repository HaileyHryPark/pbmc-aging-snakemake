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
  "NK"= "#73AF68",
  "B"   = "#79629E",
  "Mono" = "#5B83BF"
)

both <- import(snakemake@input[["degb"]]) %>% mutate(gender = "both")
female <- import(snakemake@input[["degf"]]) %>% mutate(gender = "female")
male <- import(snakemake@input[["degm"]]) %>% mutate(gender = "male")

df <- bind_rows(both, female, male)
print(head(df))

## Venn diagram
pdf(snakemake@output[["venn"]], height = 3.5, width = 4)

all_venn_list <- list(Both = df %>% filter(gender == "both") %>% pull(variable) %>% unique(), 
		Female = df %>% filter(gender == "female") %>% pull(variable) %>% unique(),
		Male = df %>% filter(gender == "male") %>% pull(variable) %>% unique())
plot(ggVennDiagram(all_venn_list, label_alpha = 0, edge_size = 1.5, set_size = 6, label_size = 5, set_color=c("Both" = "grey", "Female" = "#E15566", "Male" = "#4981BF"), label = "both", fill = "white") + labs(title = "All") + theme(legend.position = "none") + scale_fill_gradient(low = "white", high = "white"))

venns <- lapply(as.list(unique(df$celltype)), function(ct){

ctdf <- df %>% filter(celltype == ct)
  
venn_list <- list(Both = ctdf %>% filter(gender == "both") %>% pull(variable) %>% unique(), 
		Female = ctdf %>% filter(gender == "female") %>% pull(variable) %>% unique(),
		Male = ctdf %>% filter(gender == "male") %>% pull(variable) %>% unique())

plot(ggVennDiagram(venn_list, label_alpha = 0, edge_size = 1.5, set_size = 6, label_size = 5, set_color=c("Both" = "grey", "Female" = "#E15566", "Male" = "#4981BF"), label = "both", fill = "white") + labs(title = ct) + theme(legend.position = "none") + scale_fill_gradient(low = "white", high = "white"))

})
dev.off()

