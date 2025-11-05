library(rio)
library(Mfuzz)
library(dplyr)
library(tidyverse)
library(corrplot)
library(RColorBrewer)

## Functions to use
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

get_mfuzz_center <- function(data,
                             c,
                             membership_cutoff = 0.5) {
  data <-
    data@assayData$exprs
  
  membership <- c$membership
  
  membership <-
    membership %>%
    as.data.frame() %>%
    purrr::map(function(x) {
      rownames(membership)[which(x >= membership_cutoff)]
    })
  
  centers <-
    membership %>%
    purrr::map(function(x) {
      apply(data[x, , drop = FALSE], 2, mean)
    }) %>%
    dplyr::bind_rows() %>%
    as.data.frame()
  centers
  
}

plot_clusters <- function(eset, clusters, time_labels = NULL, ylim = NULL, main_prefix = "Cluster") {
  exprs <- eset@assayData$exprs
  
  # If not provided, take time labels from eset
  if (is.null(time_labels)) {
    time_labels <- eset@phenoData@data$time
  }
  
  # Ensure ylim is consistent across all clusters
  if (is.null(ylim)) {
    ylim <- range(exprs, na.rm = TRUE)
  }
  
  par(mfrow = c(1, length(clusters)))
  
  for (i in seq_along(clusters)) {
    vars <- clusters[[i]]
    mat <- exprs[vars, , drop = FALSE]
    
    matplot(
      t(mat), type = "l", lty = 1, col = rgb(0, 0, 0, 0.2),
      xaxt = "n", ylim = ylim, main = paste(main_prefix, i),
      xlab = "Age", ylab = "Scaled Expression"
    )
    axis(1, at = match(c(20, 40, 60, 80), time_labels), labels = c(20, 40, 60, 80))
    lines(colMeans(mat), col = "red", lwd = 2)
  }
}

## Determine cluster number
clust_num <- import(snakemake@input[["table"]])
pdf(snakemake@output[["cnumplot"]], width = 5, height = 4)
cnum <- FindElbowPoint(values = clust_num$distance, index = clust_num$k, sign = 1, xlab = "Cluster number", ylab = "Min. centroid distance")
dev.off()

## Import eset
eset <- table2eset(filename = snakemake@input[["mfuzz_mat"]])

m1 <- mestimate(eset)
print("M1:")
print(m1)
print(cnum)

## Cluster
set.seed(123)
clust <- mfuzz(eset = eset, c = cnum, m = m1)

center <- get_mfuzz_center(data = eset, c = clust, membership_cutoff = 0.5)
rownames(center) <- paste("Cluster", rownames(center), sep = ' ')

cor_mat <- cor(t(center))
pdf(snakemake@output[["corrplot1"]], width = 8, height = 8)
corrplot(corr = cor_mat, type = "upper", diag = T, order = "hclust", hclust.method = "ward.D",
         col = colorRampPalette(colors = rev(brewer.pal(n = 11, name = "Spectral")))(n = 100), number.cex = .7, addCoef.col = "black", tl.col="black")
dev.off()

pdf(snakemake@output[["mfuzzplot1"]], width = 26, height = 2.5)
orig_clusters <- lapply(seq_len(cnum), function(k) names(which(clust$cluster == k)))
plot_clusters(eset, orig_clusters, main_prefix = "Original Cluster")
dev.off()


## ------------------------------
## First merging pass
## ------------------------------
merge_high_corr_clusters <- function(center, eset, prefix = "MergedCluster_", cor_threshold = 0.8) {
  cor_mat <- cor(t(center))
  high_cor_pairs <- which(cor_mat > cor_threshold & upper.tri(cor_mat), arr.ind = TRUE)

  merge_groups <- list()
  for (i in seq_len(nrow(high_cor_pairs))) {
    a <- rownames(center)[high_cor_pairs[i, 1]]
    b <- rownames(center)[high_cor_pairs[i, 2]]
    found_group <- FALSE
    for (g in seq_along(merge_groups)) {
      if (a %in% merge_groups[[g]] || b %in% merge_groups[[g]]) {
        merge_groups[[g]] <- unique(c(merge_groups[[g]], a, b))
        found_group <- TRUE
        break
      }
    }
    if (!found_group) {
      merge_groups[[length(merge_groups) + 1]] <- c(a, b)
    }
  }

  cluster_map <- setNames(rownames(center), rownames(center))
  new_label_id <- 1
  for (grp in merge_groups) {
    label <- paste0(prefix, new_label_id)
    cluster_map[grp] <- label
    new_label_id <- new_label_id + 1
  }
  remaining <- setdiff(rownames(center), unlist(merge_groups))
  for (r in remaining) {
    cluster_map[r] <- paste0(prefix, new_label_id)
    new_label_id <- new_label_id + 1
  }

  merged_centers <- sapply(unique(cluster_map), function(label) {
    vars <- names(cluster_map)[cluster_map == label]
    colMeans(center[vars, , drop = FALSE])
  })
  merged_centers <- t(merged_centers)
  return(list(cluster_map = cluster_map, merged_centers = merged_centers))
}

# ----- Pass 1 -----
merge1 <- merge_high_corr_clusters(center, eset, prefix = "MergedCluster_", cor_threshold = 0.8)
center_merged1 <- merge1$merged_centers
cluster_map1 <- merge1$cluster_map

# ----- Pass 2 -----
merge2 <- merge_high_corr_clusters(center_merged1, eset, prefix = "MergedCluster2_", cor_threshold = 0.8)
center_merged2 <- merge2$merged_centers
cluster_map2 <- merge2$cluster_map

# Combine mapping (original → first merge → second merge)
combined_map <- setNames(cluster_map2[cluster_map1], names(cluster_map1))

## ------------------------------
## Update membership and output final merged clusters
## ------------------------------
membership_df <- clust$membership %>% as.data.frame()
colnames(membership_df) <- rownames(center)

var_cluster_list <- apply(membership_df, 1, function(x) {
  clusters <- names(x)[x >= 0.5]
  if (length(clusters) == 0) return(NA)
  paste(clusters, collapse = ";")
})

var_merged_list <- mapply(function(orig_str) {
  if (is.na(orig_str)) return(NA)
  orig_clusters <- unlist(strsplit(orig_str, ";"))
  merged <- unique(combined_map[orig_clusters])
  paste(merged, collapse = ";")
}, var_cluster_list, USE.NAMES = FALSE)

var_cluster_df <- data.frame(
  variable = rownames(membership_df),
  initial_clusters = var_cluster_list,
  merged_clusters  = var_merged_list,
  stringsAsFactors = FALSE
)

export(var_cluster_df, snakemake@output[["var_cluster_df"]])

## ------------------------------
## Plot and correlation of second-pass merged clusters
## ------------------------------
pdf(snakemake@output[["mfuzzplot2"]], width = 2*length(unique(var_cluster_df$merged_clusters)), height = 2.5)
merged_clusters <- split(var_cluster_df$variable, var_cluster_df$merged_clusters)
plot_clusters(eset, merged_clusters, main_prefix = "Merged Cluster (2-pass)")
dev.off()

merged_cor2 <- cor(t(center_merged2))
pdf(snakemake@output[["corrplot2"]], width = 5, height = 5)
corrplot(
  corr = merged_cor2, type = "upper", diag = TRUE,
  order = "hclust", hclust.method = "ward.D",
  col = colorRampPalette(rev(brewer.pal(11, "Spectral")))(100),
  number.cex = 0.7, addCoef.col = "black", tl.col="black"
)
dev.off()

