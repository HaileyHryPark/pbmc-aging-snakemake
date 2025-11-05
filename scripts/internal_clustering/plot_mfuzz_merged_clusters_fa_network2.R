library(rio)
library(tidyr)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(circlize)
library(reshape2)

library(igraph)
library(GGally)
library(network)
library(sna)
library(intergraph)
library(ggnetwork)


celltype_level <- c(
  "CD4 T",
  "CD8 T",
  "NK",
  "B",
  "Mono"
)
cluster_level = c("Early\nincrease", "Early\ndecrease", "Continuous\ndecrease", "Irregular\nchange", "Late\nincrease", "Continuous\nincrease")

### Functions
calculateJaccardIndex <- function(x,y){

	x <- unlist(strsplit(x,"/"))
	y <- unlist(strsplit(y,"/"))
	return(length(intersect(x,y))/length(unique(c(x,y))))
}

PlotNetworkFA <- function(df){

node <- df %>% mutate(node_name = term)
node$id <- as.character(1:nrow(node))
node <- node %>% select(id, node_name, cluster, type, Category, color, qvalue, Count, geneID)
print(head(node))

net_df <- lapply(list("female","male"), function(t){

node2 <- node %>% filter(type == t)
		
n <- nrow(node2)
w <- matrix(NA, nrow = n, ncol = n)
colnames(w) <- rownames(w) <- node2$id
	
for(i in seq_len(n-1)){
  for(j in (i+1):n){
    w[i,j] <- calculateJaccardIndex(node2[i, "geneID"], node2[j, "geneID"])
  }
}
		
wd <- melt(w)
wd <- wd[wd[,1] != wd[,2],]
wd <- wd[!is.na(wd[,3]),]
colnames(wd) <- c("node1","node2","Jaccard")
wd <- wd[wd$Jaccard > 0.25,]
wd$type = t
print(head(wd))

return(wd)

})

net <- bind_rows(net_df)
print(head(net))

#fa.net <- network(net[, 1:2], directed = FALSE, vertices = data.frame(id = node$id))
#print(fa.net)

fa.net <- network.initialize(nrow(node), directed = FALSE)
network.vertex.names(fa.net) <- node$id 

# Now add edges
add.edges(
  fa.net,
  tail = match(net$node1, node$id),
  head = match(net$node2, node$id)
)

vnames <- network.vertex.names(fa.net)
matched_nodes <- node %>% filter(id %in% vnames)

fa.net %v% "type" <- matched_nodes$type[match(vnames, matched_nodes$id)]
fa.net %v% "cluster" <- matched_nodes$cluster[match(vnames, matched_nodes$id)]
fa.net %v% "color" <- matched_nodes$color[match(vnames, matched_nodes$id)]
fa.net %v% "count" <- matched_nodes$Count[match(vnames, matched_nodes$id)]
fa.net %v% "Category" <- matched_nodes$Category[match(vnames, matched_nodes$id)]

set.seed(123)
ggplot(ggnetwork(fa.net), 
		aes(x, y, xend = xend, yend = yend)) +
	geom_edges(color = "grey21", alpha = 0.7) +
	geom_nodes(aes(fill = Category, size = count), shape = 21, color = "grey21", stroke = 0.5) + 
	#scale_color_identity()+
	scale_fill_manual(name = NULL, values = c(
					"RNA biosynthetic process" = "#E64B35", 
					"Translation" = "#4DBBD5", 
					"Proteostasis" = "#0073C2", 
					"Chromosome organization" = "#F39B7F", 
					"OXPHOS/Energy metabolism" = "#EFC000", 
					"Mitochondria" = "#B09C85",
					"Immune response" = "#00A087", 
					"Immune cell differentiation/activation" = "#91D1C2", 
					"Antigen processing and presentation" = "#8491B4", 
					"Signaling" = "#DC0000", 
					"Apoptosis" = "#A73030", 
					"Actin fiber organization" = "#7E6148", 
					"Cellular transport" = "#CD534C", 
					"Others" = NA)) +
	facet_grid(type~factor(cluster, levels = c("Early\nincrease", "Early\ndecrease", "Continuous\ndecrease", "Irregular\nchange","Late\nincrease", "Continuous\nincrease")), labeller = "label_value", switch= "y")+
	labs(x = "", y = "", color = "Gender")+
	theme_linedraw(base_size = 15)+
	theme(panel.grid = element_blank(), legend.position = "bottom", axis.text = element_blank(), axis.ticks = element_blank())

}

## Main
fa_res <- import(snakemake@input[["table"]])
annot <- import(snakemake@input[["annot"]])
annot <- annot %>% mutate(color = case_match(category,
					"RNA biosynthetic process" ~ "#E64B35", 
					"Translation" ~ "#4DBBD5", 
					"Proteostasis" ~ "#0073C2", 
					"Chromosome organization" ~ "#F39B7F", 
					"OXPHOS/Energy metabolism" ~ "#EFC000", 
					"Mitochondria" ~ "#B09C85",
					"Immune response" ~ "#00A087", 
					"Immune cell differentiation/activation" ~ "#91D1C2", 
					"Antigen processing and presentation" ~ "#8491B4", 
					"Signaling" ~ "#DC0000", 
					"Apoptosis" ~ "#A73030", 
					"Actin fiber organization" ~ "#7E6148", 
					"Cellular transport" ~ "#CD534C", 
					"Others" ~ "grey",
					.default = "grey"))


fa_res <- fa_res %>% filter(qvalue < 0.01, cluster != "", fa_celltype != "All celltype", type %in% c("female","male"))
gobp_fa_res <- merge(fa_res, annot, by.y = "term", by.x = "Description", all.x = T) %>% rename(Category = category)
gobp_fa_res <- gobp_fa_res %>% mutate(Category = ifelse(is.na(Category), "Others", Category))
print(head(gobp_fa_res))
print(table(gobp_fa_res$cluster))
print(table(gobp_fa_res$fa_celltype))
print(table(gobp_fa_res$cluster, gobp_fa_res$fa_celltype))

pdf(snakemake@output[["plot"]], width = 17.5, height = 7)

PlotNetworkFA(gobp_fa_res %>% filter(fa_celltype == "CD4 T"))
PlotNetworkFA(gobp_fa_res %>% filter(fa_celltype == "CD8 T"))
PlotNetworkFA(gobp_fa_res %>% filter(fa_celltype == "NK"))
PlotNetworkFA(gobp_fa_res %>% filter(fa_celltype == "B"))
PlotNetworkFA(gobp_fa_res %>% filter(fa_celltype == "Mono"))

dev.off()

