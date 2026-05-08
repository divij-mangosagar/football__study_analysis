library(tidyverse)      # Data manipulation and visualization
library(lme4)           # Linear mixed models (if needed later)
library(broom)          # Tidy model outputs
library(marginaleffects) # Simple slopes and interactions


# FUNCTION 1: Manual ANOVA (One-Way)
# Purpose: Compare means of a continuous variable across 3+ groups
# Inputs:
#   data - data frame containing your data
#   outcome - column name of dependent variable (e.g., "cbcl_change")
#   group - column name of grouping variable (e.g., "sports_group")
# Output: List containing ANOVA table, F-statistic, p-value, effect size, group means

manual_aov <- function(data, outcome, group) {
  
  # Extract vectors
  y <- data[[outcome]]
  g <- as.factor(data[[group]])
  group_levels <- levels(g)
  k <- length(group_levels)
  N <- length(y)
  
  # Remove missing values
  complete_cases <- !is.na(y) & !is.na(data[[group]])
  y <- y[complete_cases]
  g <- g[complete_cases]
  N <- length(y)
  
  # Group statistics
  group_means <- tapply(y, g, mean, na.rm = TRUE)
  group_sizes <- tapply(y, g, length)
  grand_mean <- mean(y, na.rm = TRUE)
  
  # Sum of Squares Between
  SSB <- sum(group_sizes * (group_means - grand_mean)^2)
  
  # Sum of Squares Within
  SSW <- 0
  for (i in 1:k) {
    group_data <- y[g == group_levels[i]]
    group_mean <- group_means[i]
    SSW <- SSW + sum((group_data - group_mean)^2, na.rm = TRUE)
  }
  
  # Total Sum of Squares
  SST <- SSB + SSW
  
  # Degrees of freedom
  between_df <- k - 1
  within_df <- N - k
  total_df <- N - 1
  
  # Mean Squares
  MSB <- SSB / between_df
  MSW <- SSW / within_df
  
  # F-statistic and p-value
  F_stat <- MSB / MSW
  p_value <- 1 - pf(F_stat, between_df, within_df)
  
  # Effect size (eta-squared)
  eta_sq <- SSB / SST
  
  # Return results
  results <- list(
    summary_table = data.frame(
      Source = c("Between", "Within", "Total"),
      SS = c(SSB, SSW, SST),
      df = c(between_df, within_df, total_df),
      MS = c(MSB, MSW, NA),
      F = c(F_stat, NA, NA),
      p = c(p_value, NA, NA),
      eta_sq = c(eta_sq, NA, NA)
    ),
    F_statistic = F_stat,
    p_value = p_value,
    eta_squared = eta_sq,
    group_means = group_means,
    grand_mean = grand_mean,
    group_levels = group_levels,
    sample_sizes = group_sizes
  )
  
  class(results) <- "manual_aov"
  return(results)
}

# Print method for manual_aov
print.manual_aov <- function(x) {
  cat("\n")
  cat("Manual One-Way ANOVA Results\n")
  cat("=============================\n\n")
  print(x$summary_table, row.names = FALSE, digits = 4)
  cat("\n")
  cat("F-statistic:", round(x$F_statistic, 4), "on", 
      length(x$group_levels) - 1, "and", sum(x$sample_sizes) - length(x$group_levels), 
      "df, p-value:", format.pval(x$p_value, digits = 4), "\n")
  cat("Eta-squared:", round(x$eta_squared, 4), "\n")
  cat("\nGroup Means:\n")
  for (i in 1:length(x$group_levels)) {
    cat("  ", x$group_levels[i], ":", round(x$group_means[i], 4), 
        "(n =", x$sample_sizes[i], ")\n")
  }
}

# FUNCTION 2: Manual t-test (Independent Samples, Welch Correction)
# Purpose: Compare means of a continuous variable between two groups
# Inputs:
#   data - data frame containing your data
#   outcome - column name of dependent variable (e.g., "cbcl_change")
#   group - column name of grouping variable (e.g., "sports_group")
#   group1 - value in group column for first group (e.g., "A")
#   group2 - value in group column for second group (e.g., "B")
# Output: List containing t-statistic, df, p-value, CI, effect size

manual_ttest <- function(data, outcome, group, group1, group2) {
  
  # Extract vectors
  y <- data[[outcome]]
  g <- data[[group]]
  
  # Filter to the two groups and remove missing
  y1 <- y[g == group1 & !is.na(y)]
  y2 <- y[g == group2 & !is.na(y)]
  
  # Sample statistics
  n1 <- length(y1)
  n2 <- length(y2)
  mean1 <- mean(y1)
  mean2 <- mean(y2)
  var1 <- var(y1)
  var2 <- var(y2)
  
  # Standard error of the difference
  se <- sqrt(var1/n1 + var2/n2)
  
  # t-statistic
  t_stat <- (mean1 - mean2) / se
  
  # Degrees of freedom (Welch-Satterthwaite)
  numerator <- (var1/n1 + var2/n2)^2
  denominator <- (var1/n1)^2/(n1 - 1) + (var2/n2)^2/(n2 - 1)
  df <- numerator / denominator
  
  # p-value (two-tailed)
  p_value <- 2 * (1 - pt(abs(t_stat), df))
  
  # 95% Confidence interval
  t_critical <- qt(0.975, df)
  ci_lower <- (mean1 - mean2) - t_critical * se
  ci_upper <- (mean1 - mean2) + t_critical * se
  
  # Cohen's d (pooled SD for equal variances assumption)
  pooled_sd <- sqrt(((n1 - 1) * var1 + (n2 - 1) * var2) / (n1 + n2 - 2))
  cohens_d <- (mean1 - mean2) / pooled_sd
  
  # Return results
  results <- list(
    t_statistic = t_stat,
    df = df,
    p_value = p_value,
    conf_interval = c(ci_lower, ci_upper),
    estimate = mean1 - mean2,
    group_means = c(mean1, mean2),
    group_names = c(group1, group2),
    cohens_d = cohens_d,
    sample_sizes = c(n1, n2)
  )
  
  class(results) <- "manual_ttest"
  return(results)
}

# Print method for manual_ttest
print.manual_ttest <- function(x) {
  cat("\n")
  cat("Welch Two Sample t-test\n")
  cat("-----------------------\n")
  cat("t =", round(x$t_statistic, 4), ", df =", round(x$df, 2), 
      ", p =", format.pval(x$p_value, digits = 4), "\n")
  cat("95% CI: [", round(x$conf_interval[1], 4), ", ", 
      round(x$conf_interval[2], 4), "]\n")
  cat("Estimated difference in means:", round(x$estimate, 4), "\n")
  cat("Group means:\n")
  cat("  ", x$group_names[1], ":", round(x$group_means[1], 4), 
      "(n =", x$sample_sizes[1], ")\n")
  cat("  ", x$group_names[2], ":", round(x$group_means[2], 4), 
      "(n =", x$sample_sizes[2], ")\n")
  cat("Cohen's d:", round(x$cohens_d, 4), "\n")
}


# FUNCTION 4: Workflow Wrapper (Read CSV, Compute Change Scores, Run ANOVA on All Measures)
# Purpose: Complete analysis pipeline from raw CSV to results
# Inputs:
#   file_path - path to your CSV file
#   participant_id_col - column name for participant ID
#   time_col - column name indicating timepoint (e.g., "time" with values "T1", "T2")
#   measure_cols - vector of outcome measure column names (e.g., c("cbcl", "stroop", "nback"))
#   group_col - column name for grouping variable (e.g., "sports_group")
# Output: List containing change scores, ANOVA results for each measure, and winning measure

run_analysis <- function(file_path, 
                         participant_id_col, 
                         time_col, 
                         measure_cols, 
                         group_col) {
  
  # Step 1: Read the CSV file
  cat("Reading data from:", file_path, "\n")
  raw_data <- read.csv(file_path, stringsAsFactors = FALSE)
  cat("Loaded", nrow(raw_data), "rows and", ncol(raw_data), "columns\n\n")
  
  # Step 2: Pivot to wide format (one row per participant)
  # Assumes the data is in long format with a time column
  wide_data <- raw_data %>%
    select(all_of(c(participant_id_col, time_col, group_col, measure_cols))) %>%
    pivot_wider(
      id_cols = all_of(c(participant_id_col, group_col)),
      names_from = all_of(time_col),
      values_from = all_of(measure_cols),
      names_sep = "_"
    )
  
  # Step 3: Compute change scores for each measure
  change_scores <- wide_data
  for (measure in measure_cols) {
    t1_col <- paste(measure, "T1", sep = "_")
    t2_col <- paste(measure, "T2", sep = "_")
    change_col <- paste(measure, "change", sep = "_")
    
    if (t1_col %in% colnames(wide_data) && t2_col %in% colnames(wide_data)) {
      change_scores[[change_col]] <- wide_data[[t2_col]] - wide_data[[t1_col]]
    } else {
      warning(paste("Could not find", t1_col, "or", t2_col, "for measure", measure))
    }
  }
  
  # Step 4: Run ANOVA on each change score
  anova_results_list <- list()
  for (measure in measure_cols) {
    change_col <- paste(measure, "change", sep = "_")
    if (change_col %in% colnames(change_scores)) {
      anova_result <- manual_aov(change_scores, change_col, group_col)
      anova_results_list[[measure]] <- data.frame(
        measure = measure,
        F_statistic = anova_result$F_statistic,
        p_value = anova_result$p_value,
        eta_squared = anova_result$eta_squared
      )
    }
  }
  
  # Combine ANOVA results into a single data frame
  anova_results <- do.call(rbind, anova_results_list)
  
  # Step 5: Identify the winning measure (largest F-statistic)
  winning_measure <- anova_results$measure[which.max(anova_results$F_statistic)]
  
  # Step 6: Return results
  results <- list(
    raw_data = raw_data,
    wide_data = wide_data,
    change_scores = change_scores,
    anova_results = anova_results,
    winning_measure = winning_measure,
    measures_tested = measure_cols
  )
  
  class(results) <- "analysis_results"
  return(results)
}

# Print method for analysis_results
print.analysis_results <- function(x) {
  cat("\n")
  cat("Analysis Results\n")
  cat("================\n\n")
  cat("Measures tested:", paste(x$measures_tested, collapse = ", "), "\n")
  cat("\nANOVA Summary:\n")
  print(x$anova_results, digits = 4)
  cat("\n")
  cat("Winning measure (largest F-statistic):", x$winning_measure, "\n")
}


results <- run_analysis(
   file_path = "FootballStudy_data_labels_4.20.26.csv",
   participant_id_col = "record.id",
   time_col = "Repeat.Instrument",  # Or whichever column indicates T1 vs T2
   measure_cols = c("cbcl", "stroop", "nback", "reading"),  # Your actual column names
   group_col = "sports_group"  # Your grouping column
)

print(results)

ttest_result <- manual_ttest(
   results$change_scores,
   outcome = paste(results$winning_measure, "change", sep = "_"),
   group = "sports_group",
   group1 = "A",
   group2 = "B"
)
 print(ttest_result)

 Run regression on the winning measure
 lm_result <- manual_lm(
   results$change_scores,
   outcome = paste(results$winning_measure, "change", sep = "_"),
   predictor = "baseline_t"  # Your brain predictor column
)
print(lm_result)