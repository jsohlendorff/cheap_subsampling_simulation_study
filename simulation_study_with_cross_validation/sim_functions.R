### sim_functions.R --- 
#----------------------------------------------------------------------
## Author: Johan Sebastian Ohlendorff
## Created: Jan  9 2026 (11:44) 
## Version: 
## Last-Updated: Jan 13 2026 (14:52) 
##           By: Johan Sebastian Ohlendorff
##     Update #: 43
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
simple_sim <- function(n, intervention = NULL){
    sex <- rbinom(n, 1, 0.5)
    L_0 <- rnorm(n)
    if (!is.null(intervention)) {
        A_0 <- as.numeric(as.character(intervention))
    } else {
        A_0 <- rbinom(n, 1, plogis(0.5 * L_0))
    }
    
    Y_1 <- rbinom(n, 1, plogis(-0.5 * A_0 + 0.2 * L_0))
    data.table(id = 1:n,
               sex = sex,
               L_0 = L_0,
               A_0 = factor(A_0, levels = c(0, 1)),
               Y_1 = Y_1,
               Censored_1 = "uncensored",
               Dead_1 = 0,
               L_1 = 1,
               A_1 = 1,
               Y_2 = 1,
               Censored_2 = "uncensored",
               Dead_2 = 0)
}

my_super_learner <- function(character_formula, ...){
    superlearn(folds = 10,
               learners = list("learn_glm",
                               "my_ranger" = list(learner_fun = "learn_ranger", max.depth = 60, min.node.size = 300)),
               id_variable = "id",
               outcome_target_level = "1",
               character_formula = character_formula,
               outcome_variable = as.character(as.formula(character_formula)[2]),
               ...)
}

run_rtmle_setup <- function(data, learner, eta, verbose = FALSE, b){
    x <- rtmle_init(intervals = 2,
                    name_id = "id",
                    name_outcome = "Y",
                    name_competing = "Dead",
                    name_censoring = "Censored",
                    censored_label = "censored")
    x <- add_baseline_data(x, data = data[, .(id, sex)])
    x <- add_wide_data(x, outcome_data = data[, .(id, Y_1, Dead_1, Censored_1, Y_2, Dead_2, Censored_2)],
                       timevar_data = list(A = data[, .(id, A_0, A_1)], L = data[, .(id, L_0, L_1)]))
    x <- protocol(x,name = "Always_A",intervention = data.frame("A" = factor("1",levels = c("0","1"))),verbose = verbose)
    x <- protocol(x,name = "Never_A",intervention = data.frame("A" = factor("0",levels = c("0","1"))),verbose = verbose)
    x <- prepare_data(x)
    x <- target(x,name = "Outcome_risk",estimator = "tmle",protocols = c("Always_A","Never_A"))
    x <- model_formula(x)
    x <- run_rtmle(x,
                   learner = learner,
                   time_horizon = 1,
                   verbose = verbose)
    res_list <- list()
    for (ssp in eta){
        y <- copy(x)
        y <- cheap_bootstrap(y,
                             B = b,
                             M = floor(ssp*NROW(y$prepared_data)),
                             which_subsets = "Main_analysis",
                             replace = FALSE)
        y$estimate$Cheap_bootstrap$Main_analysis$subsample_size_percentage <- ssp
        y$estimate$Cheap_bootstrap$Main_analysis$method <- "subsampling"
        setnames(y$estimate$Cheap_bootstrap$Main_analysis, old = c("Bootstrap_estimate", "Bootstrap_standard_error", "Bootstrap_lower", "Bootstrap_upper"),
        new = c("Estimate", "Standard_error", "Lower", "Upper"))
        res_list[[as.character(ssp)]] <- y$estimate$Cheap_bootstrap$Main_analysis
    }
    y <- copy(x)
    y <- cheap_bootstrap(y,
                         B = b,
                         M = NROW(y$prepared_data),
                         which_subsets = "Main_analysis",
                         replace = TRUE)
    y$estimate$Cheap_bootstrap$Main_analysis$subsample_size_percentage <- NA
    y$estimate$Cheap_bootstrap$Main_analysis$method <- "non-parametric"
    setnames(y$estimate$Cheap_bootstrap$Main_analysis, old = c("Bootstrap_estimate", "Bootstrap_standard_error", "Bootstrap_lower", "Bootstrap_upper"),
             new = c("Estimate", "Standard_error", "Lower", "Upper"))
    res_list[["non-parametric"]] <- y$estimate$Cheap_bootstrap$Main_analysis
    res_main <- x$estimate$Main_analysis
    res_main$subsample_size_percentage <- NA
    res_main$method <- "influence_curve"
    res_main$B <- NA
    res_main[, P_value := NULL]
    res_main[, Main_Estimate := Estimate]
    ## set order columns to match
    setcolorder(res_main, names(res_list[[1]]))
    res_list[["influence_curve"]] <- res_main
    rbindlist(res_list)
}

plot_coverage <- function(coverage){
    ## duplicate values of B = NA and subsample_size_percentage = NA to make horizontal lines for coverage
    coverage_ic <- coverage[method == "influence_curve"]
    coverage_ic <- coverage[, .(B, subsample_size_percentage)][
      , coverage_ic[rep(1, .N)], by = .(B, subsample_size_percentage)
    ][, c(names(coverage_ic)), with = FALSE]

    temp <- list()
    for (eta in unique(coverage$subsample_size_percentage)) {
        temp2 <- coverage[method == "non-parametric"]
        temp2$subsample_size_percentage <- eta
        temp[[as.character(eta)]]<- temp2
    }
    coverage_non <- rbindlist(temp)
    coverage_sub <- coverage[method == "subsampling"]

    coverage <- rbind(coverage_ic, coverage_non, coverage_sub)
    coverage <- na.omit(coverage)
    ggplot(coverage, aes(x = B, y = coverage, color = method)) +
        facet_grid(~ subsample_size_percentage, labeller = labeller(subsample_size_percentage = function(x) paste0("eta = ", x))) +
        geom_line() +
        geom_hline(yintercept = 0.95, linetype = "dashed", color = "black") +
        labs(
            x = "Number of Bootstrap Samples (B)",
            y = "Coverage Probability",
            color = "Method"
        ) +
        theme_minimal()
}

plot_bootstrap <- function(res,
                           cbbPalette = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#D55E00", "#0072B2", "#CC79A7", "#F0E442"), # color scheme for ggplot in summary_plots
                           ln_width = 1.3, # line width for ggplot in summary_plots
                           b_lower_plot = 5, # dont show the first 1:(b_lower_plot-1) bootstrap samples in the plot with m fixed
                           fix_eta_plot = 0.632) { 
  ## Options for plots
  update_geom_defaults("line", aes(linewidth = ln_width))

  ## m fixed plot; vary B; coverage
  cov_bootstrap <- res[(subsample_size_percentage == fix_eta_plot | is.na(subsample_size_percentage)) & (B >= b_lower_plot | is.na(B))]
  cov_bootstrap$method <- factor(cov_bootstrap$method)
  levels(cov_bootstrap$method) <- c(
    "Standard error based \non asymptotics",
    "Cheap Bootstrap",
    paste0("Cheap Subsampling (η = ", scales::percent(fix_eta_plot, accuracy=0.1), ")")
  )
    ## reorder the levels so that cheap_bootstrpa first, then subsampling, then standard error based
    cov_bootstrap$method <- factor(cov_bootstrap$method,
                                   levels = c(
                                       "Cheap Bootstrap",
                                     paste0("Cheap Subsampling (η = ", scales::percent(fix_eta_plot, accuracy=0.1), ")"),
                                     "Standard error based \non asymptotics"
                                     ))
  ## duplicate values of B = NA to make horizontal lines for coverage
  cov_ic <- cov_bootstrap[method == "Standard error based \non asymptotics"]$coverage
  cov_ic <- data.table(B = 1:max(cov_bootstrap$B, na.rm = TRUE),
                       coverage = cov_ic,
                       method = "Standard error based \non asymptotics",
                       subsample_size_percentage = NA,
                       width = NA)
  cov_bootstrap <- rbind(cov_bootstrap, cov_ic)
                         
  ggplot(cov_bootstrap, aes(x = B, y = coverage, color = method)) +
    geom_line() +
    geom_hline(yintercept = 0.95, linetype = "dashed", linewidth = ln_width) +
    ylab("Coverage") +
    xlab("Number of bootstrap samples (B)") +
    scale_color_manual(name = "Confidence interval \nmethod", values = cbbPalette) +
    scale_y_continuous(labels = scales::percent, limits = c(0.86, 0.965)) +
    theme_bw()
}

######################################################################
### sim_functions.R ends here
