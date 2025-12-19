summarize_survival_signal <- function(data, verbose = TRUE, seed = 123, n_bins = 3) {
  require(survival)
  require(dplyr)

  stopifnot(all(c("time", "status") %in% colnames(data)))

  ## Fit full Cox model
  formula <- as.formula(Surv(time, status) ~ .)
  fit_full <- coxph(formula, data = data)
  lp <- predict(fit_full, type = "lp")

  ## Concordance index (C-index)
  cindex <- summary(fit_full)$concordance
  cstat <- cindex[["C"]]
  cstat.se <- cindex[["se(C)"]]

  ## Null model C-index
  shuffled <- data
  shuffled$time <- sample(shuffled$time)
  fit_null <- coxph(formula, data = shuffled)
  cindex_null <- summary(fit_null)$concordance[["C"]]

  ## Likelihood Ratio Test
  fit_null0 <- coxph(Surv(time, status) ~ 1, data = data)
  lr_test <- anova(fit_null0, fit_full, test = "LRT")

  ## Survival by risk groups
  groups <- cut(lp, breaks = quantile(lp, probs = seq(0, 1, length.out = n_bins + 1)),
                include.lowest = TRUE, labels = paste0("G", 1:n_bins))
  surv_km <- survfit(Surv(time, status) ~ groups, data = data)

  ## Summary output
  diagnostics <- list(
    concordance = list(
      full_model = cstat,
      se = cstat.se,
      null_model = cindex_null
    ),
    likelihood_ratio_test = list(
      chisq = lr_test$Chisq[2],
      df = lr_test$Df[2],
      pvalue = lr_test$`P(>|Chi|)`[2]
    ),
    risk_groups = list(
      risk_grouping = groups,
      survival_fit = surv_km
    )
  )

  if (verbose) {
    cat("Concordance (C-index):\n")
    cat(sprintf("  Full model: %.3f ± %.4f\n", cstat, cstat.se))
    cat(sprintf("  Null model (shuffled time): %.3f\n\n", cindex_null))

    cat("Likelihood Ratio Test:\n")
    print(lr_test)

    cat("\nKaplan-Meier curves by risk group shown below:\n")
    plot(surv_km, col = 1:n_bins, lwd = 2, xlab = "Time", ylab = "Survival")
    legend("topright", legend = levels(groups), col = 1:n_bins, lwd = 2)
  }

  invisible(diagnostics)
}

summarize_survival_signal <- function(data, seed = 123) {
  stopifnot(all(c("time", "status") %in% names(data)))
  require(survival)

  # Fit Cox model using all features
  formula <- as.formula(Surv(time, status) ~ .)
  cox_full <- coxph(formula, data = data)
  cox_null <- coxph(Surv(time, status) ~ 1, data = data)

  # C-index
  c_full <- summary(cox_full)$concordance[["C"]]
  c_se <- summary(cox_full)$concordance[["se(C)"]]

  # C-index under null (shuffled time)
  data_shuffled <- data
  data_shuffled$time <- sample(data$time)
  cox_shuffled <- coxph(formula, data = data_shuffled)
  c_null <- summary(cox_shuffled)$concordance[["C"]]

  # LRT between null and full model
  lr_test <- anova(cox_null, cox_full, test = "LRT")
  lr_pval <- lr_test$`Pr(>|Chi|)`[2]

  # Output: interpretable summary
  return(data.frame(
    `C-index` = round(c_full, 3),
    `C-index SE` = round(c_se, 3),
    `C-index (null)` = round(c_null, 3),
    `LRT p-value` = signif(lr_pval, 3)
  ))
}

compare_data_signals <- function(model_data_list) {
  results <- lapply(model_data_list, summarize_survival_signal)
  df <- do.call(rbind, results)
  rownames(df) <- names(model_data_list)
  as.data.frame(df)
}
