print(.libPaths())
.libPaths("resources/r_package")
print(.libPaths())


library(rio)
library(dplyr)
library(tidyverse)
library(DEswan)
library(car)

# =====================================================
# CONFIG
# =====================================================

set.seed(123)

WINDOW_CENTERS <- c(20, 60)
WINDOW_SIZE <- 20
N_PER_WINDOW <- 200
ALPHA <- 0.05

# =====================================================
# LOAD DATA
# =====================================================

data <- import(snakemake@input[["data"]])

cat("\n=========================\n")
cat("RAW DATA SUMMARY\n")
cat("=========================\n")

cat("Rows:", nrow(data), "\n")
cat("Age summary:\n")
print(summary(data$age))

cat("Sex distribution:\n")
print(table(data$sex))

# =====================================================
# WINDOW SAMPLER (DEBUG VERSION)
# =====================================================

sample_window <- function(df, wc, window_size, n_per_window) {

  cat("\n====================================\n")
  cat("SAMPLING WINDOW:", wc, "\n")
  cat("====================================\n")

  window_df <- df %>%
    filter(
      age >= (wc - window_size),
      age != wc,
      age <= (wc + window_size)
    )

  cat("Available samples in window:", nrow(window_df), "\n")

  if (nrow(window_df) < 20) {
    cat("Too few samples → returning NULL\n")
    return(NULL)
  }

  cat("Age distribution in window:\n")
  print(summary(window_df$age))

  cat("Left / Right split BEFORE sampling:\n")
  cat("LEFT:", sum(window_df$age < wc), "\n")
  cat("RIGHT:", sum(window_df$age > wc), "\n")

  strata <- window_df %>%
    group_by(sex, dataset) %>%
    group_split()

  strata <- Filter(function(x) nrow(x) > 0, strata)

  cat("Number of strata:", length(strata), "\n")

  per_group <- floor(n_per_window / length(strata))

  sampled <- lapply(strata, function(x) {

    idx <- sample(
      seq_len(nrow(x)),
      size = min(per_group, nrow(x)),
      replace = TRUE
    )

    x[idx, ]
  })

  sampled <- bind_rows(sampled)

  # top-up
  if (nrow(sampled) < n_per_window) {

    extra_idx <- sample(
      seq_len(nrow(window_df)),
      size = n_per_window - nrow(sampled),
      replace = TRUE
    )

    sampled <- bind_rows(sampled, window_df[extra_idx, ])
  }

  # enforce fixed size
  sampled <- sampled[
    sample(seq_len(nrow(sampled)), n_per_window, replace = FALSE),
  ]

  cat("\nFINAL SAMPLED WINDOW:", wc, "\n")
  cat("Rows:", nrow(sampled), "\n")

  cat("Sex distribution:\n")
  print(table(sampled$sex))

  cat("Dataset distribution:\n")
  print(table(sampled$dataset))

  cat("Age summary (sampled):\n")
  print(summary(sampled$age))

  cat("Left / Right AFTER sampling:\n")
  cat("LEFT:", sum(sampled$age < wc), "\n")
  cat("RIGHT:", sum(sampled$age > wc), "\n")

  return(sampled)
}

# =====================================================
# BUILD COMBINED DATASET
# =====================================================

batch_data <- list()

for (wc in WINDOW_CENTERS) {

  sampled_df <- sample_window(
    data,
    wc,
    WINDOW_SIZE,
    N_PER_WINDOW
  )

  if (!is.null(sampled_df)) {

    sampled_df$window_center <- wc
    batch_data[[as.character(wc)]] <- sampled_df
  }
}

combined_df <- bind_rows(batch_data)

cat("\n====================================\n")
cat("COMBINED DATASET\n")
cat("====================================\n")

cat("Rows:", nrow(combined_df), "\n")
cat("Window centers:\n")
print(table(combined_df$window_center))

cat("Age distribution:\n")
print(summary(combined_df$age))

cat("Sex distribution:\n")
print(table(combined_df$sex))

cat("Dataset distribution:\n")
print(table(combined_df$dataset))

# =====================================================
# EXPRESSIONS + COVARIATES
# =====================================================

genes <- setdiff(
  colnames(combined_df),
  c("rowname", "age", "sex", "dataset", "ethnicity", "window_center")
)

expr <- combined_df[, genes]

cat("\nExpression matrix:\n")
cat("Dimensions:", dim(expr), "\n")

cat("Column types:\n")
print(table(sapply(expr, class)))

qt <- combined_df$age

covariates <- data.frame(
  dataset = factor(combined_df$dataset)
)

# =====================================================
# QT.TMP DIAGNOSTIC (CRITICAL STEP)
# =====================================================

cat("\n====================================\n")
cat("QT.TMP DIAGNOSTIC (DEswan logic check)\n")
cat("====================================\n")

for (wc in WINDOW_CENTERS) {

  qt.tmp <- rep(NA, nrow(combined_df))

  qt.tmp[
    combined_df$age < wc &
    combined_df$age >= (wc - WINDOW_SIZE)
  ] <- 0

  qt.tmp[
    combined_df$age > wc &
    combined_df$age <= (wc + WINDOW_SIZE)
  ] <- 1

  cat("\nWINDOW:", wc, "\n")
  print(table(qt.tmp, useNA = "ifany"))

  cat("Effective samples in window:\n")
  print(sum(!is.na(qt.tmp)))

  cat("Age range in this window:\n")
  print(range(combined_df$age[!is.na(qt.tmp)]))
}

# =====================================================
# DEswan DEBUG WRAPPER
# =====================================================

DEswan_debug <- function(data.df, qt, window.center, buckets.size, covariates = NULL) {

  cat("\n====================================\n")
  cat("ENTERING DEswan\n")
  cat("====================================\n")

  window.center <- sort(unique(window.center))

  pvalues.tot <- NULL
  coefficients.tot <- NULL

  for (k in seq_along(window.center)) {

    wc <- window.center[k]

    cat("\n####################################\n")
    cat("WINDOW:", wc, "\n")
    cat("####################################\n")

    pvalues <- NULL
    coefficients <- NULL

    n_success <- 0
    n_fail <- 0

    for (i in seq_len(ncol(data.df))) {

      qt.tmp <- rep(NA, length(qt))

      qt.tmp[
        qt < wc &
        qt >= (wc - buckets.size)
      ] <- 0

      qt.tmp[
        qt > wc &
        qt <= (wc + buckets.size)
      ] <- 1

      qt.tmp <- factor(qt.tmp)

      if (length(levels(qt.tmp)) < 2) {
        next
      }

      formula <- if (is.null(covariates)) {
        "data.df[, i] ~ qt.tmp"
      } else {
        paste(
          "data.df[, i] ~ qt.tmp +",
          paste("covariates$", colnames(covariates), collapse = " + ")
        )
      }

      if(i == 1){
        cat("\nCOVARIATE vs QT.TMP CHECK\n")
        print(table(qt.tmp, covariates$dataset))
      }

      if(i <= 5){
      
        gene_vals <- data.df[, i]
      
        cat("\nGENE:", colnames(data.df)[i], "\n")
      
        cat("Effective samples after NA removal:\n")
        print(sum(!is.na(gene_vals)))
      
        cat("qt.tmp distribution:\n")
        print(table(qt.tmp, useNA="always"))
      }

      fit <- try(
        glm(as.formula(formula), family = gaussian),
        silent = TRUE
      )

      if (inherits(fit, "try-error")) {
        n_fail <- n_fail + 1
      } else {
        n_success <- n_success + 1
      }
    }

    cat("\nRESULT SUMMARY\n")
    cat("Success:", n_success, "\n")
    cat("Fail:", n_fail, "\n")

  }

  return(list(status = "debug_run_complete"))
}

# =====================================================
# RUN DEswan
# =====================================================

res <- DEswan_debug(
  data.df = expr,
  qt = qt,
  window.center = WINDOW_CENTERS,
  buckets.size = WINDOW_SIZE,
  covariates = covariates
)

cat("\nDONE\n")
print(str(res))

export(data.frame(), snakemake@output[["res"]])
