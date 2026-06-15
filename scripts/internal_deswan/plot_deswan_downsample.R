library(rio)
library(dplyr)
library(ggplot2)
library(svglite)

# -----------------------------
# INPUT
# -----------------------------
both_df <- import(snakemake@input[["res_b"]]) %>% mutate(group = "both")
female_df <- import(snakemake@input[["res_f"]]) %>% mutate(group = "female")
male_df <- import(snakemake@input[["res_m"]]) %>% mutate(group = "male")

df <- bind_rows(both_df, female_df, male_df)

# ensure correct types
df$age_threshold <- as.numeric(df$age_threshold)
df$n_sig <- as.numeric(df$n_sig)
df$group <- factor(df$group, levels = c("both", "female", "male"))

# -----------------------------
# SUMMARY STATISTICS
# -----------------------------
summary_df <- df %>%
  group_by(group, age_threshold) %>%
  summarise(
    mean_nsig = mean(n_sig, na.rm = TRUE),
    sd_nsig   = sd(n_sig, na.rm = TRUE),
    q05       = quantile(n_sig, 0.05, na.rm = TRUE),
    q95       = quantile(n_sig, 0.95, na.rm = TRUE),
    .groups = "drop"
  )


p <- ggplot(summary_df, aes(x = age_threshold, y = mean_nsig, color = group, fill = group)) +
  geom_ribbon(aes(ymin = q05, ymax = q95), alpha = 0.15, color = NA) +
  facet_wrap(~group, ncol = 3) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  labs(
    x = "Age (years)",
    y = "Number of DEFs (q < 0.05)"
  ) +
  scale_color_manual(values = c("both" = "grey", "female" = "#E15566", "male" = "#4981BF")) +
  scale_fill_manual(values = c("both" = "grey", "female" = "#E15566", "male" = "#4981BF")) +
  theme_linedraw(base_size = 17)+
  theme(panel.grid = element_blank(), strip.text = element_text(size=17))

# -----------------------------
# SAVE OUTPUTS
# -----------------------------
ggsave(snakemake@output[["plot"]], p,   width = 11, height = 3.5)
