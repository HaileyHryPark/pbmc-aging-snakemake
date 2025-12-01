library(rio)
library(tidyr)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(circlize)
library(reshape2)

library(igraph)
library(ggraph)
library(scatterpie)
library(GGally)
library(network)
library(sna)
library(intergraph)
library(ggnetwork)


celltype_cols <- c("CD4 T" = "#D2533B", "CD8 T" = "#E6974D", "NK" = "#73AF68", "B" = "#79629E", "Mono" = "#5B83BF")
celltypes <- factor(names(celltype_cols), levels = names(celltype_cols))

cluster_level = c("Early\nincrease", "Early\ndecrease", "Continuous\ndecrease", "Irregular\nchange", "Late\nincrease", "Continuous\nincrease")

### Functions
jaccard <- function(a, b) {
  length(intersect(a,b)) / length(union(a,b))
}

PlotGSEANetwork <- function(df, title){

	terms <- unique(df$term)
	sig_mat <- matrix(0, nrow = length(terms), ncol = length(celltypes),
                  dimnames = list(terms, celltypes))
	for (i in seq_len(nrow(df))) {
		sig_mat[df$term[i], df$celltype[i]] <- 1
	}

	sig_df <- as.data.frame(sig_mat) %>% rownames_to_column("term")

	# leading edge list per term
	term2genes <- df %>%
	  group_by(term) %>%
	  summarise(genes = list(unique(unlist(gene_name))))
	
	# pairwise Jaccard
	edges <- expand.grid(term1 = terms, term2 = terms, 
	                     stringsAsFactors = FALSE) %>%
	  filter(term1 < term2) %>%
	  left_join(term2genes, by=c("term1"="term")) %>%
	  rename(genes1 = genes) %>%
	  left_join(term2genes, by=c("term2"="term")) %>%
	  rename(genes2 = genes) %>%
	  mutate(jacc = map2_dbl(genes1, genes2, jaccard)) %>%
	  filter(jacc > 0.25)
	
	g <- graph_from_data_frame(d = edges %>% select(term1, term2), 
	                           vertices = sig_df, directed = FALSE)
	
	# Create a ggraph layout object
	set.seed(123)
	ggraph_layout <- ggraph(g, layout = "kk")
	
	# Extract node positions from ggraph
	node_positions <- ggraph_layout$data %>%
	  select(name, x, y) %>%
	  rename(term = name)
	
	# Merge with sig_df
	plot_df <- node_positions %>% left_join(sig_df, by="term")
	
	set.seed(123)
	p1 <- ggraph(g, layout = "kk") +
	  geom_edge_link(color="grey80", alpha=.7, width = 2) +
	  scatterpie::geom_scatterpie(
	    data = plot_df,
	    aes(x = x, y = y, r = 0.18),   # constant node size
	    cols = celltypes,
	    color = NA
	  ) +
	  scale_fill_manual(values = celltype_cols) +
	  theme_void() +
	  theme(legend.position = "right") +
	  labs(title = title) +
	  coord_equal()

	set.seed(123)
	p2 <- ggraph(g, layout = "kk") +
	  geom_edge_link(color="grey80", alpha=.7, width = 2) +
	  scatterpie::geom_scatterpie(
	    data = plot_df,
	    aes(x = x, y = y, r = 0.18),   # constant node size
	    cols = celltypes,
	    color = NA
	  ) +
	  geom_text(
	    data = plot_df,
	    aes(x = x, y = y, label = term),
	    size = 3
	  ) +
	  scale_fill_manual(values = celltype_cols) +
	  theme_void() +
	  theme(legend.position = "right") +
	  labs(title = title) +
	  coord_equal()
	
	p <- ggarrange(p1, p2, ncol = 2, nrow = 1)
	plot(p)
}

## Main
res <- import(snakemake@input[["res"]]) %>%
	filter(db == "Reactome", qvalues < 0.01) %>% 
	rename(term = Description)

pdf(snakemake@output[["plots"]], width = 12, height = 6)
lapply(as.list(unique(res$type)), function(t){

	df_down <- res %>% filter(type == t, NES < 0) %>% mutate(gene_name = strsplit(gene_name, "/"))
	df_up <- res %>% filter(type == t, NES > 0) %>% mutate(gene_name = strsplit(gene_name, "/"))

	PlotGSEANetwork(df_down, paste(t, "Downregulated"))	
	PlotGSEANetwork(df_up, paste(t, "Upregulated"))	

})
dev.off()
