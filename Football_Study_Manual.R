library(tidyverse)
library(lme4)
library(broom)
library(marginaleffects)


# Manual ANOVA function
manual_aov <- function(data, outcome, group) {
  # data: data frame
  # outcome: column name of dependent variable (quoted string)
  # group: column name of grouping variable (quoted string)
  
  # Step 1: Extract vectors
  y <- data[[outcome]]
  g <- as.factor(data[[group]])  # Ensure group is factor
  group_levels <- levels(g)
  k <- length(group_levels)       # Number of groups
  N <- length(y)                  # Total sample size
  
  # Step 2: Calculate grand mean (mean of all observations)
  grand_mean <- mean(y, na.rm = TRUE)
  
  # Step 3: Calculate group means and group sizes
  group_means <- tapply(y, g, mean, na.rm = TRUE)
  group_sizes <- tapply(y, g, length)
  
  # Step 4: Calculate Sum of Squares Between (SSB)
  # SSB = sum over groups of [group_size * (group_mean - grand_mean)^2]
  SSB <- sum(group_sizes * (group_means - grand_mean)^2)
  
  # Step 5: Calculate Sum of Squares Within (SSW)
  # SSW = sum over groups of sum over observations of (y - group_mean)^2
  SSW <- 0
  for (i in 1:k) {
    group_data <- y[g == group_levels[i]]
    group_mean <- group_means[i]
    SSW <- SSW + sum((group_data - group_mean)^2, na.rm = TRUE)
  }
  
  # Step 6: Calculate Total Sum of Squares (SST)
  SST <- SSB + SSW
  # Alternative: sum((y - grand_mean)^2)
  
  # Step 7: Calculate degrees of freedom
  between_df <- k - 1      # df between groups
  within_df <- N - k       # df within groups
  total_df <- N - 1        # df total
  
  # Step 8: Calculate Mean Squares
  MSB <- SSB / between_df
  MSW <- SSW / within_df
  
  # Step 9: Calculate F-statistic
  F_stat <- MSB / MSW
  
  # Step 10: Calculate p-value
  p_value <- 1 - pf(F_stat, between_df, within_df)
  
  # Step 11: Calculate effect size (eta-squared)
  eta_sq <- SSB / SST
  
  # Step 12: Prepare and return results
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
    grand_mean = grand_mean
  )
  
  return(results)
}

# Manual t-test function (independent samples, Welch correction)
manual_ttest <- function(data, outcome, group, group1, group2) {
  # data: data frame
  # outcome: column name of dependent variable (quoted string)
  # group: column name of grouping variable (quoted string)
  # group1: value in group column for first group
  # group2: value in group column for second group
  
  # Step 1: Extract data for each group
  y <- data[[outcome]]
  g <- data[[group]]
  
  y1 <- y[g == group1]
  y2 <- y[g == group2]
  
  # Step 2: Remove missing values
  y1 <- y1[!is.na(y1)]
  y2 <- y2[!is.na(y2)]
  
  # Step 3: Calculate sample statistics
  n1 <- length(y1)
  n2 <- length(y2)
  mean1 <- mean(y1)
  mean2 <- mean(y2)
  var1 <- var(y1)
  var2 <- var(y2)
  
  # Step 4: Calculate standard error of the difference
  se <- sqrt(var1/n1 + var2/n2)
  
  # Step 5: Calculate t-statistic
  t_stat <- (mean1 - mean2) / se
  
  # Step 6: Calculate degrees of freedom (Welch correction)
  numerator <- (var1/n1 + var2/n2)^2
  denominator <- (var1/n1)^2/(n1 - 1) + (var2/n2)^2/(n2 - 1)
  df <- numerator / denominator
  
  # Step 7: Calculate p-value (two-tailed)
  p_value <- 2 * (1 - pt(abs(t_stat), df))
  
  # Step 8: Calculate confidence interval (95%)
  alpha <- 0.05
  t_critical <- qt(1 - alpha/2, df)
  ci_lower <- (mean1 - mean2) - t_critical * se
  ci_upper <- (mean1 - mean2) + t_critical * se
  
  # Step 9: Calculate Cohen's d effect size
  # Pooled standard deviation (for equal variances assumption)
  pooled_sd <- sqrt(((n1 - 1) * var1 + (n2 - 1) * var2) / (n1 + n2 - 2))
  cohens_d <- (mean1 - mean2) / pooled_sd
  
  # Step 10: Prepare results
  results <- list(
    t_statistic = t_stat,
    df = df,
    p_value = p_value,
    conf_interval = c(ci_lower, ci_upper),
    estimate = mean1 - mean2,
    group_means = c(mean1, mean2),
    names = c(group1, group2),
    cohens_d = cohens_d,
    n = c(n1, n2)
  )
  
  class(results) <- "manual_ttest"
  return(results)
}

# Print method for manual_ttest
print.manual_ttest <- function(x) {
  cat("\n")
  cat("Welch Two Sample t-test\n")
  cat("-----------------------\n")
  cat("t =", round(x$t_statistic, 4), ", df =", round(x$df, 2), ", p =", round(x$p_value, 4), "\n")
  cat("95% CI: [", round(x$conf_interval[1], 4), ",", round(x$conf_interval[2], 4), "]\n")
  cat("Estimated difference in means:", round(x$estimate, 4), "\n")
  cat("Group means:\n")
  cat("  ", x$names[1], ":", round(x$group_means[1], 4), "(n =", x$n[1], ")\n")
  cat("  ", x$names[2], ":", round(x$group_means[2], 4), "(n =", x$n[2], ")\n")
  cat("Cohen's d:", round(x$cohens_d, 4), "\n")
}

