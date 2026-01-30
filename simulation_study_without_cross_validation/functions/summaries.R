prepare_summary <- function(res, true_value) {
  setnames(res, "b", "B")

  ## Define the coverage and width of the confidence interval
  res$true_val <- true_value
  res$cov <- 1 * (res$Lower < res$true_val & res$Upper > res$true_val)
  res$cov_b <- 1 * (res$cheap_lower < res$true_val & res$cheap_upper > res$true_val)
  res$width_ic <- res$Upper - res$Lower
  res$width_b <- res$cheap_upper - res$cheap_lower
  res$rel_width <- res$width_b / res$width_ic

  ## Move stuff around to make it easier to calculate both the coverage for the influence curve CI and the bootstrap CI
  res_cov <- melt(res,
    id.vars = c("type_confidence_interval", "type", "B", "eta", "sample_size"),
    measure.vars = c("cov", "cov_b")
  )
  res_cov[variable == "cov", type_confidence_interval := "influence_curve"]

  ## Calculate coverage
  res_cov <- res_cov[, .("coverage" = mean(value)),
    by = c("type_confidence_interval", "type", "B", "eta", "sample_size"),
  ]

  ## Influence curve based coverage should be aggregated across eta and B
  res_cov[type_confidence_interval == "influence_curve", "coverage" := mean(coverage), by = "sample_size"]

  ## Relative width to influence curve-based CI
  res_width <- res[, .("rel_width" = mean(rel_width)),
    by = c("Target", "Protocol", "Time_horizon", "Estimator", "type_confidence_interval", "B", "eta", "sample_size")
  ]

  ## Return list with the results
  list(
    res_cov = res_cov,
    res_width = res_width
  )
}

plot_bootstrap <- function(summary_results,
                           cbbPalette = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#D55E00", "#0072B2", "#CC79A7", "#F0E442"), # color scheme for ggplot in summary_plots
                           ln_width = 1.3, # line width for ggplot in summary_plots
                           b_lower_plot = 5, # dont show the first 1:(b_lower_plot-1) bootstrap samples in the plot with m fixed
                           fix_eta_plot = 0.632, # which eta to fix for the plot where we vary B in summary plots
                           fix_sample_size = 8000) { ## fix the sample size
  res_cov <- summary_results$res_cov

  ## Options for plots
  update_geom_defaults("line", aes(linewidth = ln_width))

  ## m fixed plot; vary B; coverage
  cov_bootstrap <- res_cov[(eta == fix_eta_plot | is.na(eta)) & B >= b_lower_plot & sample_size == fix_sample_size]
  cov_bootstrap$type_confidence_interval <- factor(cov_bootstrap$type_confidence_interval)
  levels(cov_bootstrap$type_confidence_interval) <- c(
    "Cheap Bootstrap",
    "Cheap Subsampling",
    "Standard error based \non influence function",
    "Standard error based \non non-parametric \nbootstrap samples"
  )

  ## remove the "Standard error based \non non-parametric \nbootstrap samples" from the plot
  cov_bootstrap <- cov_bootstrap[type_confidence_interval != "Standard error based \non non-parametric \nbootstrap samples"]

  ggplot(cov_bootstrap, aes(x = B, y = coverage, color = type_confidence_interval)) +
    geom_line() +
    geom_hline(yintercept = 0.95, linetype = "dashed", linewidth = ln_width) +
    ylab("Coverage") +
    xlab("Number of bootstrap samples (B)") +
    scale_color_manual(name = "Confidence interval \nmethod", values = cbbPalette) +
    scale_y_continuous(labels = scales::percent, limits = c(0.935, 0.965)) +
    theme_bw()
}
