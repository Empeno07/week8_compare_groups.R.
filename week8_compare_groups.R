rm(list = ls())

compare_groups <- function(data, value, group, alpha = 0.05) {

  # Get the numeric column by name
  y <- data[[value]]
  # Turn the group column into a factor (categories)
  g <- factor(data[[group]])
  # Count the groups to choose between t-test and ANOVA
  k <- nlevels(g)
  # Stop if there are fewer than 2 groups
  if (k < 2) stop("Need at least 2 groups to compare.")

  # Fit aov; its residuals are used for the normality check
  fit <- aov(y ~ g)
  # Shapiro-Wilk on the residuals, not the raw data
  p_norm <- shapiro.test(residuals(fit))$p.value
  # Bartlett test for equal spread (equal variances)
  p_var <- bartlett.test(y ~ g)$p.value
  # Both assumptions hold only if both p-values are above alpha
  ok <- p_norm > alpha && p_var > alpha

  # No post-hoc table unless Tukey runs
  post <- NULL
  # Becomes TRUE when a non-parametric test is used
  fallback <- FALSE

  if (k == 2) {
    if (p_norm > alpha) {
      # Welch t-test when normality holds
      res <- t.test(y ~ g)
      test_name <- "Welch t-test"
      stat_lab <- "t"
      stat <- unname(res$statistic)
      p_val <- res$p.value
    } else {
      # Wilcoxon rank-sum when normality fails
      res <- wilcox.test(y ~ g)
      test_name <- "Wilcoxon rank-sum test"
      stat_lab <- "W"
      stat <- unname(res$statistic)
      p_val <- res$p.value
      fallback <- TRUE
    }
  } else if (ok) {
    # ANOVA when both assumptions hold
    res <- summary(fit)
    test_name <- "One-way ANOVA"
    stat_lab <- "F"
    stat <- res[[1]][["F value"]][1]
    p_val <- res[[1]][["Pr(>F)"]][1]
    # Tukey post-hoc on the aov object (not lm)
    post <- TukeyHSD(fit)
  } else {
    # Kruskal-Wallis when an assumption fails
    res <- kruskal.test(y ~ g)
    test_name <- "Kruskal-Wallis test"
    stat_lab <- "chi-squared"
    stat <- unname(res$statistic)
    p_val <- res$p.value
    fallback <- TRUE
  }

  # Main sentence: test, statistic, p-value, careful decision
  msg <- sprintf("%s: %s = %.3f, p = %.4f. %s",
                 test_name, stat_lab, stat, p_val,
                 if (p_val < alpha) "Reject H0: the groups differ."
                 else "Fail to reject H0: not enough evidence of a difference.")

  if (!is.null(post)) {
    # Tukey matrix for the group factor
    tk <- post$g
    # Pairs with adjusted p below alpha
    sig <- rownames(tk)[tk[, "p adj"] < alpha]
    # Add the significant pairs to the message
    msg <- paste(msg, if (length(sig) > 0) 
    paste0("Tukey post-hoc: significant pair(s): ", paste(sig, collapse = ", "), ".")
 else "Tukey post-hoc: no pair is significant.")

  }

  # Warn when a fallback test was used
  if (fallback) msg <- paste(msg, "Assumptions violated; used a non-parametric test.")
  # Show the assumption p-values so the choice of test is clear
  msg <- paste0(msg, sprintf(" [Shapiro p = %.4f; Bartlett p = %.5f]", p_norm, p_var))

  # Return the conclusion as text
  return(msg)
  
}

# Dataset 1: sleep (2 groups) Normality holds (p about 0.097), so Welch t gives p = 0.0794: not enough evidence.
# Design note: same 10 patients on both drugs, so a paired test (p = 0.0028) is better; here the groups are treated as independent as the task says.
print(compare_groups(sleep, "extra", "group"))

# Dataset 2: PlantGrowth (3 groups) Both assumptions hold (p about 0.438 and 0.237). ANOVA: F = 4.846, p = 0.0159.
# Tukey finds only trt2 vs trt1 (p = 0.012).
print(compare_groups(PlantGrowth, "weight", "group"))

# Dataset 3: InsectSprays (6 groups, assumptions violated) Shapiro p about 0.022 and Bartlett p about 0.00009, so the function falls back
# to Kruskal-Wallis (p near 1.5e-10) and warns about it.
print(compare_groups(InsectSprays, "count", "spray"))