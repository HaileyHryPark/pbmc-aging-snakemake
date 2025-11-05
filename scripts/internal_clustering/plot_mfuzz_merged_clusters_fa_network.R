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

lapply(as.list(unique(df$cluster)), function(clust){
	ps <- lapply(as.list(c("both","female","male")), function(g){
	plots <- lapply(as.list(c("All celltype", celltype_level)), function(ct){
		
		node <- df %>% filter(cluster == clust, fa_celltype == ct, type == g) %>% mutate(node_name = term)
		
		if(nrow(node) == 0){
			cat("No node in ", clust, " ", ct, "\n")
			return(list(ggplot() + geom_blank(), ggplot() + geom_blank()))
		}
		node$id <- 1:nrow(node)
		node <- node %>% select(id, node_name, qvalue, Count, geneID)
		
		n <- nrow(node)
		w <- matrix(NA, nrow = n, ncol = n)
		colnames(w) <- rownames(w) <- node$id
	
		for(i in seq_len(n-1)){
		  for(j in (i+1):n){
		    w[i,j] <- calculateJaccardIndex(node[i, "geneID"], node[j, "geneID"])
		  }
		}
		
		wd <- melt(w)
		wd <- wd[wd[,1] != wd[,2],]
		wd <- wd[!is.na(wd[,3]),]
		colnames(wd) <- c("node1","node2","Jaccard")
		wd <- wd[wd$Jaccard > 0.25,]

		net <- graph_from_data_frame(d = wd, vertices = node, directed = F)
		node$network_cluster <- igraph::components(net)$membership

		set.seed(123)
		netp1 <- ggnet2(net, mode = "kamadakawai", size = "Count", edge.size = sqrt(E(net)$Jaccard*5), edge.color = "grey", label = V(net)$node_name, label.size = 3, legend.position = "bottom", edge.alpha = 0.5, max_size = 10) + ggtitle(paste(ct, clust, g))

		set.seed(123)
		netp2 <- ggnet2(net, mode = "kamadakawai", size = "Count", edge.size = sqrt(E(net)$Jaccard*5), edge.color = "grey", label = V(net)$node_name, label.size = 0, legend.position = "bottom", edge.alpha = 0.5, max_size = 10) + ggtitle(paste(ct, clust, g))
	
		return(list(netp1, netp2))
	})
	
	return(list(lapply(plots, `[[`, 1), lapply(plots, `[[`, 2)))
	
	})
	both_plots <- ps[[1]]
	female_plots <- ps[[2]]
	male_plots <- ps[[3]]
	plot(ggarrange(plotlist = c(both_plots[[1]],female_plots[[1]],male_plots[[1]]), ncol = 6, nrow = 3))
	plot(ggarrange(plotlist = c(both_plots[[2]],female_plots[[2]],male_plots[[2]]), ncol = 6, nrow = 3))
})

}

## Main
fa_res <- import(snakemake@input[["table"]])
fa_res <- fa_res %>% filter(qvalue < 0.05)
print(table(fa_res$cluster))
print(table(fa_res$fa_celltype))
print(table(fa_res$cluster, fa_res$fa_celltype))

pdf(snakemake@output[["all"]], width = 90, height = 60)
PlotNetworkFA(fa_res)
dev.off()

pdf(snakemake@output[["go"]], width = 90, height = 60)
PlotNetworkFA(fa_res %>% filter(db == "GO"))
dev.off()

pdf(snakemake@output[["wp"]], width = 90, height = 60)
PlotNetworkFA(fa_res %>% filter(db == "WP"))
dev.off()

pdf(snakemake@output[["r"]], width = 90, height = 60)
PlotNetworkFA(fa_res %>% filter(db == "Reactome"))
dev.off()

