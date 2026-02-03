library(rio)
library(dplyr)
library(tidyverse)
library(purrr)

options(digits = 15)

k_folds <- 10

## Functions
precompute_zscore <- function(df, feature_cols){
  df %>% group_by(sex) %>% mutate(across(all_of(feature_cols), ~scale(.)[,1])) %>% ungroup()
}

optimize_loess_span_kfold <- function(x, y, span_range = seq(0.5, 1, 0.1), k = k_folds) {
  # Remove missing
  df <- data.frame(x, y) %>% filter(!is.na(x), !is.na(y))
  
  # Assign folds
  set.seed(123)  # reproducibility
  df <- df %>% mutate(fold = (row_number() - 1) %% k + 1)
  
  span_rmse <- map_dfr(span_range, function(s) {
    #print(s)
    rmse_vals <- numeric(k)
    
    for (fold_id in 1:k) {
      train <- df %>% filter(fold != fold_id)
      test  <- df %>% filter(fold == fold_id)
      
      fit <- try(loess(y ~ x, data = train, span = s), silent = TRUE)
      
      if (inherits(fit, "try-error")) {
        rmse_vals[fold_id] <- NA
      } else {
        preds <- try(predict(fit, newdata = test), silent = TRUE)
        if (inherits(preds, "try-error")) {
          rmse_vals[fold_id] <- NA
        } else {
          rmse_vals[fold_id] <- sqrt(mean((test$y - preds)^2, na.rm = TRUE))
        }
      }
    }
    
    data.frame(span = s, rmse = mean(rmse_vals, na.rm = TRUE))
  })
  
  return(span_rmse)
}


limma <- import(snakemake@input[["limma"]])

features <- colnames(limma %>% select(-c(rowname,age,sex,dataset,ethnicity)))
zscore <-  precompute_zscore(limma, features)
export(zscore, snakemake@output[["zscaled"]])

zscore_long <- zscore %>% 
    pivot_longer(cols = all_of(features), names_to = "feature", values_to = "zscore") %>% 
    filter(!is.na(zscore), !is.na(age))

span_res <- lapply(sort(unique(zscore_long$feature)), function(f){
	message("Processing feature: ", f)

		df_fg <- zscore_long %>% filter(feature == f)
		res <- optimize_loess_span_kfold(
      			x = df_fg$age,
      			y = df_fg$zscore,
      			k = k_folds
    		)

		res <- res %>% mutate(feature = f)
		print(head(res))

		opt_span <- res %>% filter(rmse == min(rmse, na.rm = TRUE)) %>% arrange(span) %>% slice(1) %>% pull(span) 
		print(opt_span)

		fit <- loess(zscore ~ age, data = df_fg, span = opt_span)
		new_age <- seq(19, 93)
		preds <- predict(fit, newdata = data.frame(age = new_age))
		fit_df <- data.frame(feature = f, age = new_age, fitted = preds)

		return(list(res, fit_df))
})

export(bind_rows(lapply(span_res, `[[`, 1)), snakemake@output[["span_res"]])

fitted <- bind_rows(lapply(span_res, `[[`, 2))
export(fitted, snakemake@output[["fit_res"]])

fitted_wide <- fitted %>% select(feature, age, fitted) %>% 
  pivot_wider(names_from = age, values_from = fitted) %>% column_to_rownames("feature")

time <- colnames(fitted_wide)
mfuzz_data <- rbind(time, fitted_wide)
row.names(mfuzz_data)[1] <- "time"

write.table(mfuzz_data, file = snakemake@output[["mfuzz_mat"]], sep = "\t", quote = FALSE, col.names = NA)
