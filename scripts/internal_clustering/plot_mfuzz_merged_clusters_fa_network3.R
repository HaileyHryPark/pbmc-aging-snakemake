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


celltype_cols <- c("CD4 T" = "#D2533B", "CD8 T" = "#E6974D", "B" = "#79629E", "NK" = "#73AF68", "Mono" = "#5B83BF")
celltypes <- factor(names(celltype_cols), levels = names(celltype_cols))

cluster_level = c("Early\nincrease", "Early\ndecrease", "Continuous\ndecrease", "Early\nfluctuation", "Inverted\nUshape", "Continuous\nincrease", "Late\nincrease")

### Functions
jaccard <- function(a, b) {
  length(intersect(a,b)) / length(union(a,b))
}

PlotFANetwork <- function(df, title){

	terms <- unique(df$term)
	sig_mat <- matrix(0, nrow = length(terms), ncol = length(celltypes),
                  dimnames = list(terms, celltypes))
	print(head(sig_mat))
	for (i in seq_len(nrow(df))) {
		sig_mat[df$term[i], df$fa_celltype[i]] <- 1
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

	# Extract top terms per connected component
	# Add qvalue into vertex attributes
	term_q <- df %>% 
	  distinct(term, qvalue)
	
	# Match qvalues to actual graph vertex names
	vertex_q <- term_q$qvalue[match(V(g)$name, term_q$term)]
	
	g <- igraph::set_vertex_attr(g, "qvalue", value = vertex_q)
	
	# Connected components
	comp <- igraph::components(g)
	
	# Build dataframe: term, component, qvalue, ranking
	top_terms_df <- data.frame(
	  term       = names(comp$membership),
	  component  = comp$membership,
	  qvalue     = igraph::vertex_attr(g, "qvalue")
	) %>%
	  group_by(component) %>%
	  arrange(qvalue) %>%
	  mutate(rank = row_number()) %>%
	  filter(rank <= 3) %>%
	  ungroup() %>% 
	  mutate(title = title)
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
	return(top_terms_df)

}

## Main
res <- import(snakemake@input[["table"]]) %>%
	filter(db == snakemake@params[["db"]], qvalue < 0.01, fa_celltype != "All celltype") %>% 
	select(-term) %>%
	rename(term = Description)
print(head(res))

pdf(snakemake@output[["plot"]], width = 22, height = 10)
top_terms <- lapply(as.list(unique(res$type)), function(t){

	res <- lapply(as.list(cluster_level), function(cl){
		df <- res %>% filter(type == t, cluster == cl) %>% mutate(gene_name = strsplit(gene_name, "/"))
		print(head(df))
		if(nrow(df) > 0)
			return(PlotFANetwork(df, paste(t, cl)))
		else
			return(data.frame())
	})
	return(bind_rows(res))

})
dev.off()

top <- bind_rows(top_terms)

# For female inverted Ushape
top_terms_to_plot <- top %>% filter(title == "female Inverted\nUshape", rank == 1) %>% arrange(qvalue) %>% slice_head(n = 5) %>% pull(term)

res_to_plot <- res %>% filter(type == "female", cluster == "Inverted\nUshape", term %in% top_terms_to_plot) %>% arrange(qvalue) %>% mutate(term = factor(term, levels = rev(top_terms_to_plot)))

res_to_plot <- res_to_plot  %>% 
  arrange(term, desc(qvalue))  %>% 
  mutate(group_order = forcats::fct_inorder(interaction(term, fa_celltype)),
         fa_celltype = factor(fa_celltype, levels= c("CD4 T", "CD8 T", "B", "NK", "Mono")))

p <- ggplot(res_to_plot, aes(y = term, group = group_order)) +
  geom_col(aes(x = -log10(qvalue), fill = fa_celltype), position = position_dodge2(width = 0.5, preserve = "single"), width = 0.25) +
  geom_text(
    aes(x = 0, label = term),
    hjust = 0,
    nudge_x = max(-log10(res_to_plot$qvalue)) / 80,
    nudge_y = 0.4,
    size = 6,
    lineheight = 0.95
  ) +
  geom_hline(yintercept = (1:4)+0.6, linewidth = 0.5) +
  scale_fill_manual(values = celltype_cols) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    x = "-log10(qvalue)",
    y = "GOBP Terms",
  ) +
  theme_linedraw(base_size = 18) +
  theme(
    panel.grid=element_blank(),
    axis.text.y = element_blank(),  # hide original y-axis text
    axis.ticks.y = element_blank()
  )
ggsave(snakemake@output[["plot_fiu"]], plot = p, width = 10, height = 5)

# For female continuous increase
top_terms_to_plot <- top %>% filter(title == "female Continuous\nincrease", rank == 1) %>% arrange(qvalue) %>% slice_head(n = 5) %>% pull(term)

res_to_plot <- res %>% filter(type == "female", cluster == "Continuous\nincrease", term %in% top_terms_to_plot) %>% arrange(qvalue) %>% mutate(term = factor(term, levels = rev(top_terms_to_plot)))

res_to_plot <- res_to_plot  %>% 
  arrange(term, desc(qvalue))  %>% 
  mutate(group_order = forcats::fct_inorder(interaction(term, fa_celltype)),
         fa_celltype = factor(fa_celltype, levels= c("CD4 T", "CD8 T", "B", "NK", "Mono")))

p <- ggplot(res_to_plot, aes(y = term, group = group_order)) +
  geom_col(aes(x = -log10(qvalue), fill = fa_celltype), position = position_dodge2(width = 0.5, preserve = "single"), width = 0.5) +
  geom_text(
    aes(x = 0, label = term),
    hjust = 0,
    nudge_x = max(-log10(res_to_plot$qvalue)) / 80,
    nudge_y = 0.4,
    size = 6,
    lineheight = 0.95
  ) +
  geom_hline(yintercept = (1:4)+0.6, linewidth = 0.5) +
  scale_fill_manual(values = celltype_cols) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(
    x = "-log10(qvalue)",
    y = "GOBP Terms",
  ) +
  theme_linedraw(base_size = 18) +
  theme(
    panel.grid=element_blank(),
    axis.text.y = element_blank(),  # hide original y-axis text
    axis.ticks.y = element_blank()
  )
ggsave(snakemake@output[["plot_fci"]], plot = p, width = 10, height = 5)

# For female late increase
top_terms_to_plot <- top %>% filter(title == "female Late\nincrease", rank == 1) %>% arrange(qvalue) %>% slice_head(n = 5) %>% pull(term)

res_to_plot <- res %>% filter(type == "female", cluster == "Late\nincrease", term %in% top_terms_to_plot) %>% arrange(qvalue) %>% mutate(term = factor(term, levels = rev(top_terms_to_plot)))

res_to_plot <- res_to_plot  %>% 
  arrange(term, desc(qvalue))  %>% 
  mutate(group_order = forcats::fct_inorder(interaction(term, fa_celltype)),
         fa_celltype = factor(fa_celltype, levels= c("CD4 T", "CD8 T", "B", "NK", "Mono")))

p <- ggplot(res_to_plot, aes(y = term, group = group_order)) +
  geom_col(aes(x = -log10(qvalue), fill = fa_celltype), position = position_dodge2(width = 0.5, preserve = "single"), width = 0.5) +
  geom_text(
    aes(x = 0, label = term),
    hjust = 0,
    nudge_x = max(-log10(res_to_plot$qvalue)) / 80,
    nudge_y = 0.4,
    size = 6,
    lineheight = 0.95
  ) +
  geom_hline(yintercept = (1:4)+0.6, linewidth = 0.5) +
  scale_fill_manual(values = celltype_cols) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(
    x = "-log10(qvalue)",
    y = "GOBP Terms",
  ) +
  theme_linedraw(base_size = 18) +
  theme(
    panel.grid=element_blank(),
    axis.text.y = element_blank(),  # hide original y-axis text
    axis.ticks.y = element_blank()
  )
ggsave(snakemake@output[["plot_fli"]], plot = p, width = 10, height = 5)

## For male early fluctuation
top_terms_to_plot <- top %>% filter(title == "male Early\nfluctuation", rank == 1) %>% arrange(qvalue) %>% slice_head(n = 5) %>% pull(term)

res_to_plot <- res %>% filter(type == "male", cluster == "Early\nfluctuation", term %in% top_terms_to_plot) %>% 
  mutate(term = factor(term, levels = rev(top_terms_to_plot))) %>%
  arrange(term, desc(qvalue))

res_to_plot <- res_to_plot  %>% 
  mutate(group_order = forcats::fct_inorder(interaction(term, fa_celltype)),
         fa_celltype = factor(fa_celltype, levels= c("CD4 T", "CD8 T", "B", "NK", "Mono")))

p <- ggplot(res_to_plot, aes(y = term, group = group_order)) +
  geom_col(aes(x = -log10(qvalue), fill = fa_celltype), position = position_dodge2(width = 0.5, preserve = "single"), width = 0.5) +
  geom_text(
    aes(x = 0, label = term),
    hjust = 0,
    nudge_x = max(-log10(res_to_plot$qvalue)) / 80,
    nudge_y = 0.4,
    size = 6,
    lineheight = 0.95
  ) +
  geom_hline(yintercept = (1:4)+0.6, linewidth = 0.5) +
  scale_fill_manual(values = celltype_cols) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(
    x = "-log10(qvalue)",
    y = "GOBP Terms",
  ) +
  theme_linedraw(base_size = 18) +
  theme(
    panel.grid=element_blank(),
    axis.text.y = element_blank(),  # hide original y-axis text
    axis.ticks.y = element_blank()
  )
ggsave(snakemake@output[["plot_mef"]], plot = p, width = 10, height = 5)

