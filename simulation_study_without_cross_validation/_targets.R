## Simulation study for the Cheap Subsampling confidence interval
## We vary the hyperparameters of the Cheap Subsampling confidence interval
## - Number of bootstrap samples (b=1,...,b_max)
## - Subsample percentage (etas) (m = floor(eta * n))
## - Sample size (sample_sizes)

## The functions are in the functions folder, with the following content:
## - simulation_functions.R: main simulations, cheap subsampling/bootstrap cis, functions for plots, true values
## - rtmle_functions.R: functions for LTMLE, a CheapSubsampling wrapper for rtmle

## Install necessary packages (e.g., on server)
# install.packages(c("data.table", "ggplot2", "tarchetypes", "crew", "crew.cluster","devtools", "dplyr", "tidyr", "tibble", "gt", "ggpubr"))
# devtools::install_github("jsohlendorff/CheapSubsampling")
# devtools::install_github("tagteam/rtmle")

library(targets)
library(tarchetypes)
library(crew)
library(crew.cluster)

## Try to load on server
try(setwd("/projects/biostat01/people/snf991/phd/cheap_subsampling_article/simulation_study"),
  silent = TRUE
)

## Try to load locally
try(setwd("~/phd/cheap_subsampling_simulation_study_repo/simulation_study_without_cross_validation//"), silent = TRUE)

tar_source("functions")

## Set up the crew controller
if (dir.exists("/projects/biostat01/people/snf991/phd/cheap_subsampling_simulation_study_repo/simulation_study_without_cross_validation")) {
  controller <- crew_controller_slurm(
    workers = 100,
    seconds_idle = 15,
    options_cluster = crew_options_slurm(partition = "long") # start on markov
  )
} else {
  controller <- NULL
}

## targets set packages
tar_option_set(
  packages = c(
    "data.table",
    "ggplot2",
    "CheapSubsampling",
    "rtmle",
    "tidyr",
    "dplyr",
    "tibble",
    "gt",
    "ggpubr",
    "latex2exp"
  ),
  controller = controller
)
## Intervetion
intervention <- 1

## Number of time points
time_horizon <- 2

## Sample sizes
sample_sizes <- c(250, 500, 1000, 2000, 8000)

## Bootstrap settings
etas <- c(0.5, 0.632, 0.8, 0.9) # subsample percentage
b_max <- 20 # number of bootstrap samples

## Repititions (= reps*batches)
batches <- 1 # i.e., number of slurm jobs for each value of eta
reps <- 1 # how many simulation repititions in each batch; to use parallel with ncpus argument to slurm, use rep_workers instead to make sure that the repititions are paralleized

## Vary the parameters of the simulation study
vals_subsampling <- expand.grid(
  eta = etas,
  sample_size = sample_sizes,
  type = "subsampling",
  stringsAsFactors = FALSE
)
vals_non_parametric <- expand.grid(
  eta = NA,
  sample_size = sample_sizes,
  type = "non_parametric",
  stringsAsFactors = FALSE
)
vals <- rbind(vals_subsampling, vals_non_parametric)

## main targets string applicable for each type of simulation study
list(
  ## True value
  tar_target(true_value, true_value_rtmle(
    sample_size = 10000,#10000000,
    time_horizon = time_horizon,
    intervention = intervention
  )),

  ## Smulation study
  tar_map_rep(
    sims,
    command = run_simulation(
      b_max = b_max,
      eta = eta,
      time_horizon = time_horizon,
      type = type,
      sample_size = sample_size
    ),
    values = vals,
    batches = batches,
    reps = reps,
    cue = tar_cue(mode = "never")
  ),

  ## Prepare for summary
  tar_target(summary_results, prepare_summary(sims, true_value)),

  tar_target(bootstrap_plot, plot_bootstrap(
    summary_results = summary_results,
    cbbPalette = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#D55E00", "#0072B2", "#CC79A7", "#F0E442"), # color scheme for ggplot in summary_plots
    ln_width = 0.9, # line width for ggplot in summary_plots
    b_lower_plot = 5, # dont show the first 1:(b_lower_plot-1) bootstrap samples in the plot with m
    fix_eta_plot = 0.8, # which eta to fix for the plot where we vary B in summary plots
    fix_sample_size = 8000
  ))
)
