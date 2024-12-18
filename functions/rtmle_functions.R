cheap_bootstrap_rtmle <- function(x,
                                  cheap_bootstrap_arguments,
                                  rtmle_arguments,
                                  id_name) {
  ## Recursively apply rbind
  recursive_rbind <- function(x) if (inherits(x, "list")) do.call("rbind", lapply(x, recursive_rbind)) else x

  fun <- function(data) {
    ## Assume that the ids across x$followup and x$prepared_data are ordered the same way
    matched_ids <- match(data[[id_name]], x$prepared_data[[id_name]])

    ## Subset data according to the resampled data]
    x$prepared_data <- x$prepared_data[matched_ids]
    x$followup <- x$followup[matched_ids]

    ## Run rtmle and suppress messages
    suppressMessages(res <- recursive_rbind(do.call(
      "run_rtmle",
      append(list(x = x), rtmle_arguments)
    )$estimate))

    ## For correct usage of CheapSubsampling package, no numerics are allowed in paramater_column_names
    ## Thus Time_horizon should be changed to, say, a factor
    res$Time_horizon <- as.factor(res$Time_horizon)
    res
  }

  do.call(
    "cheap_bootstrap",
    append(
      list(
        fun = fun,
        estimate_column_name = "Estimate",
        parameter_column_names = c("Target", "Protocol", "Time_horizon", "Estimator"),
        data = x$prepared_data
      ),
      cheap_bootstrap_arguments
    )
  )
}
