print(.libPaths())
.libPaths("resources/r_package")
print(.libPaths())

library(rio)
library(dplyr)
library(tidyverse)
library(DEswan)
library(car)
library("DEswan")

# =====================================================
# CONFIG
# =====================================================

set.seed(123)

WINDOW_CENTERS <- seq(40, 44, by = 2)   # your test windows
WINDOW_SIZE <- 20
N_PER_WINDOW <- 200
N_ITER <- 2
ALPHA <- 0.05

# =====================================================
# LOAD DATA
# =====================================================

data <- import(snakemake@input[["data"]])

cat("\n============================\n")
cat("FULL DATASET\n")
cat("============================\n")

cat("Rows:", nrow(data), "\n")
cat("Age summary:\n")
print(summary(data$age))

# =====================================================
# DEswan WITH INTERNAL DOWNSAMPLING
# =====================================================

DEswan_downsampled <- function(
    data.df,
    qt,
    window.center,
    buckets.size,
    covariates = NULL,
    n_per_window = 200
){

  cat("\n====================================\n")
  cat("ENTERING DEswan (DOWNSAMPLED VERSION)\n")
  cat("====================================\n")

  window.center <- sort(unique(window.center))

  pvalues.tot <- NULL
  coefficients.tot <- NULL

  for (k in seq_along(window.center)) {

    wc <- window.center[k]

    cat("\n####################################\n")
    cat("WINDOW CENTER:", wc, "\n")
    cat("####################################\n")

    # =================================================
    # STEP 1: build FULL DEswan window mask
    # =================================================

    idx_left <- which(
      qt < wc &
      qt >= (wc - buckets.size)
    )

    idx_right <- which(
      qt > wc &
      qt <= (wc + buckets.size)
    )

    idx <- c(idx_left, idx_right)

    cat("\nRaw window size:", length(idx), "\n")

    if (length(idx) < n_per_window) {
      cat("Skipping (too small)\n")
      next
    }

    # =================================================
    # STEP 2: DOWNSAMPLE WITHIN WINDOW ONLY
    # =================================================

    idx <- sample(
      idx,
      size = min(n_per_window, length(idx)),
      replace = FALSE
    )

    cat("Downsampled window size:", length(idx), "\n")

    data_sub <- data.df[idx, , drop = FALSE]
    qt_sub <- qt[idx]

    cov_sub <- if (!is.null(covariates)) {
      covariates[idx, , drop = FALSE]
    } else {
      NULL
    }

    # =================================================
    # DOWNSAMPLE DIAGNOSTICS
    # =================================================
    
    left_idx <- which(
      qt_sub < wc &
      qt_sub >= (wc - buckets.size)
    )
    
    right_idx <- which(
      qt_sub > wc &
      qt_sub <= (wc + buckets.size)
    )
    
    cat("\nDOWNSAMPLE DIAGNOSTICS\n")
    
    cat("LEFT samples:", length(left_idx), "\n")
    cat("RIGHT samples:", length(right_idx), "\n")
    
    cat("\nLEFT age summary:\n")
    print(summary(qt_sub[left_idx]))
    
    cat("\nRIGHT age summary:\n")
    print(summary(qt_sub[right_idx]))
    
    if(!is.null(cov_sub)) {
    
      cat("\nLEFT dataset distribution:\n")
      print(table(cov_sub[left_idx, 1]))
    
      cat("\nRIGHT dataset distribution:\n")
      print(table(cov_sub[right_idx, 1]))
    
    }

    # =================================================
    # STEP 3: RUN GLM PER GENE (DEswan CORE)
    # =================================================

    pvalues <- NULL
    coefficients <- NULL

    n_success <- 0
    n_fail <- 0

    for (i in seq_len(ncol(data_sub))) {

      qt.tmp <- rep(NA, length(qt_sub))

      qt.tmp[qt_sub < wc & qt_sub >= (wc - buckets.size)] <- 0
      qt.tmp[qt_sub > wc & qt_sub <= (wc + buckets.size)] <- 1

      qt.tmp <- factor(qt.tmp)

      if (length(levels(qt.tmp)) < 2) {
        next
      }

      formula <- if (is.null(cov_sub)) {
        "data_sub[, i] ~ qt.tmp"
      } else {
        paste(
          "data_sub[, i] ~ qt.tmp +",
          paste("cov_sub$", colnames(cov_sub), collapse = " + ")
        )
      }

      fit <- try(
        glm(as.formula(formula), family = gaussian),
        silent = TRUE
      )

      if (inherits(fit, "try-error")) {
        n_fail <- n_fail + 1
        next
      }

      n_success <- n_success + 1

      glm.res <- try(
        car::Anova(fit, type = "2"),
        silent = TRUE
      )

      if (!inherits(glm.res, "try-error")) {

        pvalues <- rbind(
          pvalues,
          data.frame(
            variable = colnames(data_sub)[i],
            window.center = wc,
            factor = rownames(glm.res),
            pvalue = glm.res$`Pr(>Chisq)`,
            stringsAsFactors = FALSE
          )
        )

        coefficients <- rbind(
          coefficients,
          data.frame(
            variable = colnames(data_sub)[i],
            window.center = wc,
            factor = names(coefficients(fit)),
            coefficient = coefficients(fit),
            stringsAsFactors = FALSE
          )
        )
      }
    }

    cat("\nRESULT SUMMARY\n")
    cat("Success:", n_success, "\n")
    cat("Fail:", n_fail, "\n")

    pvalues.tot <- rbind(pvalues.tot, pvalues)
    coefficients.tot <- rbind(coefficients.tot, coefficients)
  }

   pvalues.tot$factor[which(pvalues.tot$factor=="qt.tmp")]<-"qt"
   pvalues.tot$factor=gsub("^covariates\\$","",pvalues.tot$factor)
   coefficients.tot$factor[which(coefficients.tot$factor=="qt.tmp1")]<-"qt"
   coefficients.tot$factor=gsub("^covariates\\$","",coefficients.tot$factor)

  return(list(p = pvalues.tot, coeff = coefficients.tot))
}

# =====================================================
# BOOTSTRAP WRAPPER
# =====================================================

run_bootstrap_deswan <- function(data, n_iter = 2) {

  genes <- setdiff(
    colnames(data),
    c("rowname", "age", "sex", "dataset", "ethnicity")
  )

  results <- list()

  for (iter in seq_len(n_iter)) {

    cat("\n====================================\n")
    cat("ITERATION:", iter, "\n")
    cat("====================================\n")

    res <- DEswan_downsampled(
      data.df = data[, genes],
      qt = data$age,
      window.center = WINDOW_CENTERS,
      buckets.size = WINDOW_SIZE,
      covariates = data["dataset"],
      n_per_window = N_PER_WINDOW
    )

    if (!is.null(res$p)) {

      p_wide_res <- reshape.DEswan(res,parameter = 1,factor = "qt") %>% as.data.frame()
      print("p_wide_res")
      print(head(p_wide_res))
      q_res <- q.DEswan(p_wide_res, method="BH") %>% as.data.frame()
      print("q_res")
      print(head(q_res))
      sig <- q_res %>%
        pivot_longer(cols = -c(variable), names_to = "age_threshold", values_to = "qvalue") %>% 
        mutate(age_threshold = as.numeric(gsub("X", "", age_threshold))) %>%
        group_by(age_threshold) %>%
        summarise(n_sig = sum(qvalue < ALPHA, na.rm = TRUE))
      print("sig")
      print(head(sig))

      sig$iteration <- iter

      results[[iter]] <- sig
    }
  }

  bind_rows(results)
}

# =====================================================
# RUN
# =====================================================

curve_df <- run_bootstrap_deswan(
  data,
  n_iter = N_ITER
)

# =====================================================
# SAVE
# =====================================================

export(curve_df, snakemake@output[["res"]])
