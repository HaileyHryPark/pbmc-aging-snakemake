print(.libPaths())
.libPaths("resources/r_package")
print(.libPaths())


library(dplyr)
library(tidyverse)
library(rio)
library(broom)
library(qvalue)

# -----------------------------
# INPUT
# -----------------------------
data <- import(snakemake@input[["data"]])

genes <- setdiff(colnames(data),
                  c("rowname", "age", "sex", "dataset", "ethnicity"))

# -----------------------------
# CONFIG
# -----------------------------
window_centers <- seq(20, 90, by = 2)
window_size <- 20

# fixed per-window sample size AFTER downsampling
n_per_window <- 100
alpha <- 0.05

set.seed(123)

# -----------------------------
# WINDOW SAMPLER (BALANCED DOWNSAMPLE)
# -----------------------------
sample_window_balanced <- function(df, wc, window_size, n_target = 100) {

  window_df <- df %>%
    filter(age >= wc - window_size/2,
           age <= wc + window_size/2)

  if (nrow(window_df) < 30) return(NULL)

  # split by strata
  strata <- window_df %>%
    group_by(sex, dataset) %>%
    group_split()

  strata <- Filter(function(x) nrow(x) > 0, strata)

  if (length(strata) == 0) return(NULL)

  # allocate samples per stratum
  per_group <- floor(n_target / length(strata))

  sampled <- lapply(strata, function(x) {
    n_take <- min(nrow(x), per_group)

    x %>%
      slice_sample(n = n_take, replace = n_take < per_group)
  })

  sampled <- bind_rows(sampled)

  # top-up if needed
  if (nrow(sampled) < n_target) {
    extra <- window_df %>%
      slice_sample(n = n_target - nrow(sampled), replace = TRUE)

    sampled <- bind_rows(sampled, extra)
  }

  # final enforce exact size
  sampled <- sampled %>%
    slice_sample(n = n_target, replace = FALSE)

  return(sampled)
}

# -----------------------------
# MAIN ANALYSIS
# -----------------------------
run_sliding_de <- function(df, label) {

  results <- list()

  for (wc in window_centers) {

    cat("\n==============================\n")
    cat("WINDOW CENTER:", wc, "\n")
    cat("==============================\n")

    window_df <- sample_window_balanced(df, wc, window_size, n_per_window)

    if (is.null(window_df)) next

    cat("Total cells:", nrow(window_df), "\n")

    cat("\nSex distribution:\n")
    print(table(window_df$sex))

    cat("\nDataset distribution:\n")
    print(table(window_df$dataset))

    cat("\nAge summary:\n")
    print(summary(window_df$age))

    cat("Age SD:", sd(window_df$age), "\n")

    # sanity check
    if (sd(window_df$age) < 2) {
      cat("SKIP: too narrow age distribution\n")
      next
    }

    pvals <- numeric(0)
    tested_genes <- 0
    failed_genes <- 0

    for (g in genes) {

      df_gene <- window_df %>%
        select(all_of(g), dataset) %>%
        rename(expr = all_of(g))

      if (sd(df_gene$expr, na.rm = TRUE) == 0) next

      fit <- tryCatch({
        glm(expr ~ dataset, data = df_gene)
      }, error = function(e) NULL)

      if (is.null(fit)) {
        failed_genes <- failed_genes + 1
        next
      }

      p <- tryCatch({
        anova(fit, test = "Chisq")$`Pr(>Chi)`[2]
      }, error = function(e) NA)

      if (!is.na(p)) {
        pvals <- c(pvals, p)
        tested_genes <- tested_genes + 1
      }
    }

    pvals <- pvals[is.finite(pvals)]
    
    if (length(pvals) < 50) {
      cat("SKIP: too few p-values\n")
      next
    }
    
    qobj <- tryCatch({
	qvalue::qvalue(pvals)
    }, error = function(e) NULL)

    if (is.null(qobj)) {
      cat("SKIP: qvalue failed\n")
      next
    }

    qvals <- qobj$qvalues

    sig_count <- sum(qvals < 0.05, na.rm = TRUE)

    cat("\nGENE SUMMARY:\n")
    cat("tested genes:", tested_genes, "\n")
    cat("failted genes:", failted_genes, "\n")
    cat("significant genes:", sig_count, "\n")

    results[[length(results) + 1]] <- data.frame(
      window = wc,
      n_sig = sig_count,
      group = label,
      n_cells = nrow(window_df)
    )
  }

  bind_rows(results)
}

# -----------------------------
# RUN FOR GROUPS
# -----------------------------
data_f <- data %>% filter(sex == "female")
data_m <- data %>% filter(sex == "male")

res_f <- run_sliding_de(data_f, "female")
res_m <- run_sliding_de(data_m, "male")
res_b <- run_sliding_de(data, "both")

curve_df <- bind_rows(res_f, res_m, res_b)

export(curve_df, snakemake@output[["res"]])
