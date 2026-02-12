library(dplyr)
library(tidyverse)
library(rio)

gender <- snakemake@params[["gender"]]
print(gender)
f <- ifelse(gender == "both", 2, ifelse(gender == "female", 4, 3))
print(f)
internal_csv <- import(snakemake@input[["pred_i"]]) %>% filter(fold == f) %>% 
	mutate(sample_id = donor_id, disease = "normal", cohort = "internal") %>% 
	select(sample_id, donor_id, cohort, dataset, disease, actual_age, predicted_age, sex, fold)
print(head(internal_csv))
external_csv <- import(snakemake@input[["pred_e"]]) %>% filter(fold == f) %>% 
	mutate(cohort = "external") %>% 
	mutate(actual_age = as.numeric(actual_age)) %>% filter(!is.na(actual_age)) %>%
	select(sample_id, donor_id, cohort, dataset, disease, actual_age, predicted_age, sex, fold)
print(head(external_csv))

if(gender != "both"){
	external_csv <- external_csv %>% filter(sex == gender)
	print("here")
}
# -----------------------------
# Define correction function
# -----------------------------
apply_loess <- function(internal_df, external_df, span = 0.5, degree = 2) {
  
  # Ensure required columns exist
  required_cols <- c("fold", "actual_age", "predicted_age")
  stopifnot(all(required_cols %in% names(internal_df)))
  stopifnot(all(required_cols %in% names(external_df)))
  
  # Add age_diff column
  internal_df <- internal_df %>%
    mutate(age_diff = predicted_age - actual_age, type = "internal")
  
  external_df <- external_df %>%
    mutate(age_diff = predicted_age - actual_age, type = "external")
  
    # Fit loess on internal data
    loess_fit <- loess(age_diff ~ actual_age,
                       data = internal_df,
                       span = span,
                       degree = degree,
                       control = loess.control(surface = "direct"))
    
    # Predict correction term for both
    int_pred <- predict(loess_fit, newdata = internal_df)
    ext_pred <- predict(loess_fit, newdata = external_df)
    
    # Add corrected columns
    internal_df <- internal_df %>%
      mutate(
        loess_fit = int_pred,
        c_age_diff = age_diff - loess_fit,
        c_predicted_age = predicted_age - loess_fit
      )
    
    external_df <- external_df %>%
      mutate(
        loess_fit = ext_pred,
        c_age_diff = age_diff - loess_fit,
        c_predicted_age = predicted_age - loess_fit
      )
    
  return(bind_rows(internal_df, external_df))  
}

# -----------------------------
# Apply correction
# -----------------------------
result <- apply_loess(internal_csv, external_csv, span = 0.5, degree = 2)
result <- result %>% 
	mutate(age_group = ifelse(actual_age > 60, ">60", ifelse(actual_age < 40, "<40", "40-60")),
		age_accel1 = ifelse(c_age_diff > 0, "Accelerated", "Decelerated"),
		c_age_diff_z = (c_age_diff - mean(c_age_diff, na.rm = TRUE)) / sd(c_age_diff, na.rm = TRUE),
         	age_accel2 = case_when(
           		c_age_diff_z > 1 ~ "Accelerated",
           		c_age_diff_z < -1 ~ "Decelerated",
           		TRUE ~ "Intermediate"
         	))
print(head(result))

# -----------------------------
# Save results
# -----------------------------
export(result, snakemake@output[["corrected"]])


