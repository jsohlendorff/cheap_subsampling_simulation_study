## Monte Carlo SE (non-parametric) Bootstrap confidence interval
## From \hat{\Psi}_1^*, \dots, \hat{\Psi}_B^* use \hat{SE} = sd(\hat{\Psi}_1^*, \dots, \hat{\Psi}_B^*)
monte_carlo_variance_bootstrap <- function(est, boot_est, alpha) {
  se <- sd(boot_est)
  list(
    cheap_lower = est - qnorm(1 - alpha / 2) * se,
    cheap_upper = est + qnorm(1 - alpha / 2) * se,
    estimate = est
  )
}

rexpit <- function(x) rbinom(n = length(x), size = 1, prob = plogis(x))

## Function to simulate data for the simulation study
simulate_data <- function(n, time_horizon, intervene = NULL) {
  L <- sample(0:1, n, replace = TRUE)
  W_0 <- rnorm(n)
  if (!is.null(intervene)) {
    A_0 <- rep(intervene, n)
  } else {
    A_0 <- rexpit(-0.2 * W_0 + 0.4)
  }
  has_event <- is_censored <- rep(FALSE, n)
  at_risk <- rep(TRUE, n)
  for (t in 1:time_horizon) {
    ## Simulate C_t (censoring)
    ## at_risk <- !is_censored & !has_event
    if (!is.null(intervene)) {
      assign(paste0("C_", t), 1)
    } else {
      C <- rep(1, n)
      C[is_censored] <- 0
      C[at_risk] <- rexpit(3.5 + get(paste0("W_", t - 1))[at_risk])
      is_censored <- C == 0
      assign(paste0("C_", t), C)
      at_risk <- !is_censored & at_risk
    }

    ## Simulate Y_t (outcome)
    Y <- rep(0, n)
    Y[is_censored] <- NA
    Y[has_event] <- 1
    Y[at_risk] <- rexpit(-1.4 + 0.1 * get(paste0("W_", t - 1))[at_risk] - 1.5 * get(paste0("A_", t - 1))[at_risk])
    has_event <- Y == 1
    assign(paste0("Y_", t), Y)
    at_risk <- !has_event & at_risk

    ## Simulate W_t (time-varying covariate)
    W <- rep(NA, n)
    W[at_risk] <- rnorm(0.5 * get(paste0("W_", t - 1))[at_risk] + 0.2 * get(paste0("A_", t - 1))[at_risk])
    assign(paste0("W_", t), W)

    ## Simulate A_t (time-varying treatment)
    A <- rep(NA, n)
    if (!is.null(intervene)) {
      A[at_risk] <- rep(intervene, sum(at_risk))
    } else {
      A[at_risk] <- rexpit(-0.4 * get(paste0("W_", t - 1))[at_risk] + 0.8 * get(paste0("A_", t - 1))[at_risk])
    }
    A <- as.integer(A)
    assign(paste0("A_", t), A)
  }
  list(
    outcome = as.data.table(dplyr::bind_cols(pnr = 1:n, lapply(1:time_horizon, function(t) {
      dt <- data.table::data.table(
        C = get(paste0("C_", t)),
        Y = get(paste0("Y_", t))
      )
      setnames(dt, paste0(names(dt), "_", t))
      dt
    }))),
    regimen = as.data.table(dplyr::bind_cols(pnr = 1:n, A_0 = A_0, lapply(seq_len(time_horizon - 1), function(t) {
      dt <- data.table::data.table(
        A = get(paste0("A_", t))
      )
      setnames(dt, paste0(names(dt), "_", t))
      dt
    }))),
    timevarying_covariates = as.data.table(dplyr::bind_cols(pnr = 1:n, W_0 = W_0, lapply(seq_len(time_horizon - 1), function(t) {
      dt <- data.table::data.table(
        W = get(paste0("W_", t))
      )
      setnames(dt, paste0(names(dt), "_", t))
      dt
    }))),
    baseline_data = data.table(pnr = 1:n, L = L)
  )
}

## True value of the parameter for intervention
true_value_rtmle <- function(sample_size, time_horizon, intervention = 1) {
  simulated_data <- simulate_data(sample_size, time_horizon, intervene = intervention)
  mean(simulated_data$outcome[[paste0("Y_", time_horizon)]])
}

## Main function to run the simulation study
run_simulation <- function(b_max = 10,
                           eta = 0.5,
                           time_horizon = 2,
                           type = "subsampling", ## "subsampling" or "non_parametric"
                           sample_size) {
  simulated_data <- simulate_data(sample_size, time_horizon, intervene = NULL)

  ## RTMLE initialization
  x <- rtmle_init(
    intervals = time_horizon,
    name_id = "pnr",
    name_outcome = "Y",
    name_competing = NULL,
    name_censoring = "C",
    censored_levels = c(1, 0),
    censored_label = 0
  )
  
  # Add data
  x <- add_baseline_data(x, data = simulated_data$baseline_data)
  x$data$outcome_data <- simulated_data$outcome
  x$data$timevar_data <- simulated_data$timevarying_covariates
  x$data$timevar_data <- x$data$timevar_data[simulated_data$regimen, on = "pnr"]
  x <- protocol(x, name = "Always_A", intervention = data.frame("A" = factor("1", levels = c("0", "1"))))
  x <- prepare_data(x)
  
  # Specify the target parameter
  x <- target(
    x,
    name = "Outcome_risk",
    strategy = "additive",
    estimator = "tmle",
    protocols = "Always_A"
  )
  x <- model_formula(x)
  rtmle_arguments <- list(learner = "learn_glm", time_horizon = time_horizon)

  res_cheap_bootstrap <- cheap_bootstrap_rtmle(
    x = x,
    cheap_bootstrap_arguments = list(
      b = b_max,
      size = round(eta * sample_size),
      parallelize = FALSE,
      cores = 1,
      type = type
    ),
    rtmle_arguments = rtmle_arguments,
    id_name = "pnr"
  )
  ## Add the monte carlo se bootstrap confidence interval if applicable
  if (type == "non_parametric") {
    monte_carlo_confidence_interval <- list()
    for (bs in seq_len(b_max)) {
      monte_carlo_confidence_interval[[bs]] <- res_cheap_bootstrap$boot_estimates[b %in% 1:bs,
        monte_carlo_variance_bootstrap(estimate[1], b_est, 0.05),
        by = c("Target", "Protocol", "Time_horizon", "Estimator")
      ]
      monte_carlo_confidence_interval[[bs]]$b <- bs
      monte_carlo_confidence_interval[[bs]]$type_confidence_interval <- "monte_carlo_se_bootstrap"
    }
    res_cheap_bootstrap$res$type_confidence_interval <- "cheap_bootstrap"
    res_cheap_bootstrap$res <- rbind(rbindlist(monte_carlo_confidence_interval), res_cheap_bootstrap$res)
  } else {
    res_cheap_bootstrap$res$type_confidence_interval <- "cheap_subsampling"
  }
  merge(res_cheap_bootstrap$res, res_cheap_bootstrap$original_value, by = c("Target", "Protocol", "Time_horizon", "Estimator"))
}
