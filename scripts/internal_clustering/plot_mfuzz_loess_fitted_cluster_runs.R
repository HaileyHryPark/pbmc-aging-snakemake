library(tidyverse)
library(ggpubr)
library(rio)

# ----------------------------
# Load data
# ----------------------------
cluster_counts <- import(snakemake@input[["cluster_counts"]]) %>% mutate(n_clusters = as.integer(n_clusters))
all_runs <- readRDS(snakemake@input[["all_runs_rds"]])
final_centers <- import(snakemake@input[["final_centers"]], header = T)

# ----------------------------
# Plot 1: cluster number stability
# ----------------------------
p1 <- ggplot(cluster_counts, aes(x = n_clusters)) +
  geom_bar(fill = "black") +
  theme_classic() +
  labs(
    x = "Number of merged clusters",
    y = "Frequency",
  )

ggsave(
  snakemake@output[["plot1"]],
  p1,
  width = 1.3*length(unique(cluster_counts)),
  height = 3
)

# ----------------------------
# Plot 2: trajectory stability
# ----------------------------
final_long <- final_centers %>%
  pivot_longer(-cluster, names_to = "time", values_to = "value") %>%
  mutate(type = "Final", time = as.integer(time))

p2 <- ggplot() +
  geom_line(
    data = final_long,
    aes(time, value, group = cluster, color = cluster),
    linewidth = 1
  ) +
  theme_linedraw(base_size = 13) +
  theme(panel.grid.major=element_blank(), legend.position = "none") +
  scale_x_continuous(breaks=seq(20, 90, 10)) +
  labs(
    x = "Age",
    y = "Scaled expression"
  )

p21 <- ggplot() + geom_line(data = final_long, aes(time, value, group = cluster, color = cluster), linewidth = 1) + 
	theme_linedraw(base_size = 13) + theme(panel.grid.major=element_blank(), legend.position = "bottom") + 
	scale_x_continuous(breaks=seq(20, 90, 10)) +
	labs(x = "Age", y = "Scaled expression")

ggsave(
  snakemake@output[["plot2"]],
  ggarrange(p2, p21, ncol = 1, nrow = 2),
  width = 3.3,
  height = 6
)

run_long <- map_dfr(all_runs, function(r) {
  r$merged_centers %>%
    rownames_to_column("cluster") %>%
    pivot_longer(-cluster, names_to = "time", values_to = "value") %>%
    mutate(type = paste0("Seed_", r$seed), time = as.integer(time))
})

p3 <- ggplot() +
  geom_line(
    data = run_long,
    aes(time, value, group = cluster, color = cluster),
    linewidth = 1
  ) +
  facet_wrap(~ type, scales = "free_y", ncol = 5) +
  theme_linedraw(base_size = 13) +
  theme(panel.grid.major=element_blank(), legend.position = "none", strip.text = element_text(size = 16)) +
  scale_x_continuous(breaks=seq(20, 90, 10)) +
  labs(
    x = "Age",
    y = "Scaled expression"
  )

ggsave(
  snakemake@output[["plot3"]],
  p3,
  width = 16.5,
  height = 6
)

