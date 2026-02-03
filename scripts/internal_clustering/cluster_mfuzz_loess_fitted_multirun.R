library(Mfuzz)
library(dplyr)
library(tibble)
library(tidyverse)
library(purrr)
library(igraph)
library(rio)

# ----------------------------
# Parameters
# ----------------------------
seeds <- 1:10
cor_threshold <- 0.8
membership_cutoff <- 0.5

# ----------------------------
# Function: graph-based merging
# ----------------------------
merge_clusters_graph <- function(center, cor_threshold = 0.8) {
  cor_mat <- cor(t(center))
  diag(cor_mat) <- 0

  g <- graph_from_adjacency_matrix(
    cor_mat > cor_threshold,
    mode = "undirected",
    diag = FALSE
  )

  comps <- components(g)$membership
  cluster_map <- setNames(
    paste0("Merged_", comps),
    rownames(center)
  )

  merged_centers <- lapply(unique(cluster_map), function(cl) {
    vars <- names(cluster_map)[cluster_map == cl]
    colMeans(center[vars, , drop = FALSE])
  }) %>% bind_rows()

  rownames(merged_centers) <- unique(cluster_map)

  list(
    cluster_map = cluster_map,
    merged_centers = merged_centers
  )
}

# ----------------------------
# Function: graph-based merging
# ----------------------------
FindElbowPoint <- function(values, index = NULL, sign = 1, main = "Elbow Detection",xlab = "Index", ylab = "Value") {

  # Default x-values if not provided
  if (is.null(index)) {
    index <- seq_along(values)
  }

  # Check lengths
  if (length(index) != length(values)) {
    stop("'index' and 'values' must be the same length")
  }

  start <- c(index[1], values[1])
  end <- c(index[length(index)], values[length(values)])

  dist2d <- function(a, b, c) {
    v1 <- b - c
    v2 <- a - b
    abs(det(cbind(v1, v2))) / sqrt(sum(v1 * v1))  # perpendicular distance
  }

  vals <- sapply(
    seq_along(values),
    function(i) {
      sign * dist2d(c(index[i], values[i]), start, end)
    }
  )

  elbow <- which.max(vals)

  # Plot
  plot(index, values, type = "b", pch = 19, col = "steelblue",
       xlab = xlab, ylab = ylab, main = main)

  # Reference line from start to end
  slope <- (end[2] - start[2]) / (end[1] - start[1])
  abline(a = start[2] - slope * start[1], b = slope,
         col = "gray", lty = 2)

  # Elbow marker
  points(index[elbow], values[elbow], col = "red", pch = 19, cex = 1.5)
  text(index[elbow], values[elbow],
       labels = paste("Elbow:", index[elbow]), pos = 3, col = "red")

  return(index[elbow])
}


# ----------------------------
# Load data
# ----------------------------
clust_num <- import(snakemake@input[["table"]])
cnum <- FindElbowPoint(values = clust_num$distance, index = clust_num$k, sign = 1, xlab = "Cluster number", ylab = "Min. centroid distance")
print(cnum)

eset <- table2eset(snakemake@input[["mfuzz_mat"]])
set.seed(123)
m1 <- mestimate(eset)

# ----------------------------
# Run 10 times
# ----------------------------
all_runs <- lapply(seeds, function(seed) {
  set.seed(seed)

  cl <- mfuzz(eset, c = cnum, m = m1)

  membership <- cl$membership
  centers <- t(sapply(1:cnum, function(k) {
    vars <- rownames(membership)[membership[, k] >= membership_cutoff]
    colMeans(exprs(eset)[vars, , drop = FALSE])
  }))

  rownames(centers) <- paste0("Cluster_", seq_len(nrow(centers)))

  merged <- merge_clusters_graph(centers, cor_threshold)

  merged2 <- merge_clusters_graph(merged$merged_centers, cor_threshold)

  list(
    seed = seed,
    n_clusters = length(unique(merged2$cluster_map)),
    merged_centers = merged2$merged_centers,
    cluster_map = merged2$cluster_map
  )
})

# ----------------------------
# Save outputs
# ----------------------------
cluster_counts <- tibble(
  seed = seeds,
  n_clusters = map_int(all_runs, "n_clusters")
)

export(cluster_counts, snakemake@output[["cluster_counts"]])

saveRDS(all_runs, snakemake@output[["all_runs_rds"]])

