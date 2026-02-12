library(rio)
library(DMRcate)
library(limma)
library(tidyverse)
library(ggpubr)
library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

### Functions
extract_dmr_cpg_matrix <- function(dmr_ranges, dmr_index, beta_matrix, cpg_gr) {
  # Get the DMR
  dmr <- dmr_ranges[dmr_index]
  
  # Subset CpGs that overlap this DMR
  overlaps <- subsetByOverlaps(cpg_gr, dmr)
  overlaps <- sort(overlaps)
  cpg_ids <- overlaps$probe_id
  
  # Extract CpG betas
  betas <- beta_matrix[cpg_ids, ]
  
  # Return as CpG × Sample matrix
  return(betas)
}

make_long_df <- function(cpg_beta_mat, metadata){
  df <- as.data.frame(cpg_beta_mat)
  df$CpG <- rownames(df)
  
  df_long <- df %>%
    pivot_longer(-CpG, names_to = "Sample", values_to = "Beta") %>%
    left_join(metadata, by = c("Sample" = "Basename"))
  
  return(df_long)
}

plot_cpg_lines <- function(df_long) {
  
  cols <- c("<40" = "black", "40-60" = "orange", ">60" = "forestgreen")
  # summarize mean and SE per CpG per group
  df_summary <- df_long %>%
    group_by(Gender, AgeGroup, CpG, group) %>%
    summarise(
      mean_beta = mean(Beta, na.rm = TRUE),
      se_beta   = sd(Beta, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )
  
  ggplot(df_summary, aes(x = CpG, y = mean_beta, color = AgeGroup, fill = AgeGroup, group = AgeGroup)) +
    facet_wrap(~Gender, ncol = 1, nrow = 2) +
    # geom_ribbon(aes(ymin = mean_beta - se_beta,
    #                 ymax = mean_beta + se_beta),
    #             alpha = 0.2) +   # transparent SE band
    geom_line(linewidth = 0.5) +                     # mean line
    scale_color_manual(values = cols) +
    scale_fill_manual(values = cols) +
    theme_linedraw(base_size = 15) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1), panel.grid.major = element_blank(), axis.text.x.bottom = element_blank()) +
    labs(x = "CpGs", y = "Mean Beta")
}

plot_mean_box <- function(df_long){
  cols <- c("<40" = "black", "40-60" = "orange", ">60" = "forestgreen")
  df_mean <- df_long %>%
    group_by(Sample, AgeGroup, Gender) %>%
    summarize(mean_beta = mean(Beta), .groups = "drop")
  
  ggplot(df_mean, aes(x = AgeGroup, y = mean_beta, fill = AgeGroup)) +
    facet_wrap(~Gender, ncol = 1, nrow = 2) +
    geom_violin(color="black", width = 0.6) +
    geom_boxplot(fill="white", width = 0.2) +
    stat_compare_means(label = "p", comparisons = list(c("<40","40-60"),c("40-60",">60"),c("<40",">60")),
                       step.increase = 0.2, tip.length = 0, bracket.size = 0.5, vjust = -0.4) +
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.15))) +
    scale_fill_manual(values = cols) +
    theme_linedraw(base_size = 15) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1), panel.grid.major = element_blank(), legend.position = "none", axis.title.x = element_blank()) +
    labs(x = "", y = "Mean Beta")
}



### Main
# Read data
metadata <- import(snakemake@input[["metadata"]]) %>% filter(!is.na(Age), Gender != "")
metadata$Gender <- factor(metadata$Gender, levels = c("female", "male"))
metadata <- metadata %>% mutate(group = paste(Gender, AgeGroup, sep = " "), color = case_match(group, "female <40" ~ "orange", "female 40-60" ~ "magenta", "female >60" ~ "red", "male <40" ~ "black", "male 40-60" ~ "forestgreen", "male >60" ~ "blue"))

beta_all <- import(snakemake@input[["beta"]]) %>%
  column_to_rownames("probe_prefix")
beta_all <- beta_all[,metadata$Basename] %>% as.matrix()

res <- readRDS(snakemake@input[["res"]])

dmr_age <- res$dmr$AG_40to60_vs_lt40
dmr_int1 <- res$dmr$Interaction_40to60
dmr_int2 <- res$dmr$Interaction_gt60
dmr_int3 <- res$dmr$Interaction_40to60_gt60


## Pie charts for chromosomes
dmr_int1_df <- res$dmr$Interaction_40to60 %>% as.data.frame() %>% 
  mutate(chr = ifelse(seqnames == "chrX", "chrX", 
                      ifelse(seqnames == "chrY", "chrY", "autosome")))
bychr <- as.data.frame(table(dmr_int1_df$chr)) %>% arrange(Freq)
bychr$col <- c("grey", "#E1556650")
bychr$fraction <- bychr$Freq / sum(bychr$Freq)
bychr$ymax <- bychr$ymax <- cumsum(bychr$fraction)
bychr$ymin <- c(0, head(bychr$ymax, n=-1))

pdf(snakemake@output[["pie"]], width = 3, height = 3)

ggplot(bychr, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=col)) +
  geom_rect() +
  scale_fill_identity(guide = "legend") +
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  theme_void() + 
  theme(legend.position = "right")

dev.off()

dmr_int2_df <- res$dmr$Interaction_gt60 %>% as.data.frame() %>% 
  mutate(chr = ifelse(seqnames == "chrX", "chrX", 
                      ifelse(seqnames == "chrY", "chrY", "autosome")))
bychr <- as.data.frame(table(dmr_int2_df$chr)) %>% arrange(Freq)
bychr$col <- factor(c("#4981BF50","grey", "#E1556650"), levels = c("#E1556650","#4981BF50","grey"))
bychr$fraction <- bychr$Freq / sum(bychr$Freq)
bychr$ymax <- bychr$ymax <- cumsum(bychr$fraction)
bychr$ymin <- c(0, head(bychr$ymax, n=-1))

ggplot(bychr, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=col)) +
  geom_rect() +
  scale_fill_identity(guide = "legend") +
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  theme_void() + 
  theme(legend.position = "right")
dev.off()

## DMRcate plots
colors <- metadata$color
names(colors) <- metadata$group

pdf(snakemake@output[["plot1"]], width = 5, height = 4)
DMR.plot(ranges = dmr_int1, dmr = 22, CpGs = beta_all, what = "Beta", arraytype = "450K", phen.col = colors, genome = "hg19", heatmap = F)
DMR.plot(ranges = dmr_int1, dmr = 47, CpGs = beta_all, what = "Beta", arraytype = "450K", phen.col = colors, genome = "hg19", heatmap = F)
DMR.plot(ranges = dmr_int2, dmr = 193, CpGs = beta_all, what = "Beta", arraytype = "450K", phen.col = colors, genome = "hg19", heatmap = F)
DMR.plot(ranges = dmr_int2, dmr = 208, CpGs = beta_all, what = "Beta", arraytype = "450K", phen.col = colors, genome = "hg19", heatmap = F)
DMR.plot(ranges = dmr_int2, dmr = 244, CpGs = beta_all, what = "Beta", arraytype = "450K", phen.col = colors, genome = "hg19", heatmap = F)
DMR.plot(ranges = dmr_int2, dmr = 295, CpGs = beta_all, what = "Beta", arraytype = "450K", phen.col = colors, genome = "hg19", heatmap = F)
DMR.plot(ranges = dmr_int2, dmr = 293, CpGs = beta_all, what = "Beta", arraytype = "450K", phen.col = colors, genome = "hg19", heatmap = F)
# DMR.plot(ranges = dmr_age, dmr = 2, CpGs = beta_all, what = "Beta", arraytype = "450K", phen.col = colors, genome = "hg19", heatmap = F)
dev.off()

## Line plots and vln plots
anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# Ensure CpG probes only
anno <- anno[rownames(beta_all), ]

# Convert to GRanges
cpg_gr <- GRanges(
  seqnames = anno$chr,
  ranges = IRanges(start = anno$pos, width = 1),
  strand = "*",
  probe_id = rownames(anno)
)

cpg_gr

pdf(snakemake@output[["plot2"]], width = 6.5, height = 5)

lapply(list(22, 47, 193, 208, 244, 295, 293), function(dmr_index){

if(dmr_index %in% c(22, 47))
	cpg_mat <- extract_dmr_cpg_matrix(dmr_int1, dmr_index, beta_all, cpg_gr)
else
	cpg_mat <- extract_dmr_cpg_matrix(dmr_int2, dmr_index, beta_all, cpg_gr)

df_long <- make_long_df(cpg_mat, metadata)
df_long$CpG <- factor(df_long$CpG, levels = rownames(cpg_mat))
df_long <- df_long %>% mutate(AgeGroup = factor(AgeGroup, levels = c("<40", "40-60", ">60")),
                              Gender = factor(Gender, levels = c("female", "male"), labels = c("Female", "Male")))
print(head(df_long))

cols <- c("<40" = "black", "40-60" = "orange", ">60" = "forestgreen")

p1 <- plot_cpg_lines(df_long)
p2 <- plot_mean_box(df_long)

ggarrange(p1, p2, ncol = 2, nrow = 1, widths = c(4,2.5))

})

dev.off()

## Individual plots
dmr_index <- 22
cpg_mat <- extract_dmr_cpg_matrix(dmr_int1, dmr_index, beta_all, cpg_gr)

df_long <- make_long_df(cpg_mat, metadata)
df_long$CpG <- factor(df_long$CpG, levels = rownames(cpg_mat))
df_long <- df_long %>% mutate(AgeGroup = factor(AgeGroup, levels = c("<40", "40-60", ">60")),
                              Gender = factor(Gender, levels = c("female", "male"), labels = c("Female", "Male")))

cols <- c("<40" = "black", "40-60" = "orange", ">60" = "forestgreen")

p1 <- plot_cpg_lines(df_long)
p2 <- plot_mean_box(df_long)

p <- ggarrange(p1, p2, ncol = 2, nrow = 1, widths = c(6,2))
ggsave(snakemake@output[["p1"]], p, width = 8, height = 5)

dmr_index <- 293
cpg_mat <- extract_dmr_cpg_matrix(dmr_int2, dmr_index, beta_all, cpg_gr)

df_long <- make_long_df(cpg_mat, metadata)
df_long$CpG <- factor(df_long$CpG, levels = rownames(cpg_mat))
df_long <- df_long %>% mutate(AgeGroup = factor(AgeGroup, levels = c("<40", "40-60", ">60")),
                              Gender = factor(Gender, levels = c("female", "male"), labels = c("Female", "Male")))

cols <- c("<40" = "black", "40-60" = "orange", ">60" = "forestgreen")

p1 <- plot_cpg_lines(df_long)
p2 <- plot_mean_box(df_long)

p <- ggarrange(p1, p2, ncol = 2, nrow = 1, widths = c(3.5,2))
ggsave(snakemake@output[["p2"]], p, width = 5.5, height = 5)

dmr_index <- 208
cpg_mat <- extract_dmr_cpg_matrix(dmr_int2, dmr_index, beta_all, cpg_gr)

df_long <- make_long_df(cpg_mat, metadata)
df_long$CpG <- factor(df_long$CpG, levels = rownames(cpg_mat))
df_long <- df_long %>% mutate(AgeGroup = factor(AgeGroup, levels = c("<40", "40-60", ">60")),
                              Gender = factor(Gender, levels = c("female", "male"), labels = c("Female", "Male")))

cols <- c("<40" = "black", "40-60" = "orange", ">60" = "forestgreen")

p1 <- plot_cpg_lines(df_long)
p2 <- plot_mean_box(df_long)

p <- ggarrange(p1, p2, ncol = 2, nrow = 1, widths = c(3.5,2))
ggsave(snakemake@output[["p3"]], p, width = 5.5, height = 5)

dmr_index <- 244
cpg_mat <- extract_dmr_cpg_matrix(dmr_int2, dmr_index, beta_all, cpg_gr)

df_long <- make_long_df(cpg_mat, metadata)
df_long$CpG <- factor(df_long$CpG, levels = rownames(cpg_mat))
df_long <- df_long %>% mutate(AgeGroup = factor(AgeGroup, levels = c("<40", "40-60", ">60")),
                              Gender = factor(Gender, levels = c("female", "male"), labels = c("Female", "Male")))

cols <- c("<40" = "black", "40-60" = "orange", ">60" = "forestgreen")

p1 <- plot_cpg_lines(df_long)
p2 <- plot_mean_box(df_long)

p <- ggarrange(p1, p2, ncol = 2, nrow = 1, widths = c(3,2))
ggsave(snakemake@output[["p4"]], p, width = 5, height = 5)

dmr_index <- 295
cpg_mat <- extract_dmr_cpg_matrix(dmr_int2, dmr_index, beta_all, cpg_gr)

df_long <- make_long_df(cpg_mat, metadata)
df_long$CpG <- factor(df_long$CpG, levels = rownames(cpg_mat))
df_long <- df_long %>% mutate(AgeGroup = factor(AgeGroup, levels = c("<40", "40-60", ">60")),
                              Gender = factor(Gender, levels = c("female", "male"), labels = c("Female", "Male")))

cols <- c("<40" = "black", "40-60" = "orange", ">60" = "forestgreen")

p1 <- plot_cpg_lines(df_long)
p2 <- plot_mean_box(df_long)

p <- ggarrange(p1, p2, ncol = 2, nrow = 1, widths = c(3.5,2))
ggsave(snakemake@output[["p5"]], p, width = 5.5, height = 5)

