## _targets.R --- 
#----------------------------------------------------------------------
## Author: Johan Sebastian Ohlendorff
## Created: Jan  9 2026 (11:44) 
## Version: 
## Last-Updated: Jan 13 2026 (14:57) 
##           By: Johan Sebastian Ohlendorff
##     Update #: 97
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
library(targets)
library(tarchetypes)
library(crew)
library(crew.cluster)

## Try to load locally
try(setwd("~/phd/cheap_subsampling_article/simulation_study_tmle_cross_validation/"), silent = TRUE)

tar_source("sim_functions.R")

## Server set / crew setup
if (dir.exists("/projects/biostat01/people/snf991/phd/cheap_subsampling_article/simulation_study_tmle_cross_validation/")) {
    controller <- crew_controller_slurm(
        workers = 100,
        seconds_idle = 15,
        options_cluster = crew_options_slurm(partition = "long") # start on markov
    )
} else {
    controller <- crew_controller_local(workers = 16,
                                        options_local = crew_options_local(log_directory = "log"))
}

tar_option_set(
    packages = c(
        "rtmle",
        "data.table",
        "ranger",
        "reshape2",
        "ggplot2"
    ),
    controller = controller,
    error = "null" # produce a result even if the target errored out.
)

## Sample sizes
sample_size <- 1000

etas <- c(0.5, 0.632, 0.8, 0.9) # subsample percentage
b_max <- 100 # number of bootstrap samples

batches <- 200 # i.e., number of slurm jobs for each value of eta
reps <- 10 # how many simulation repititions in each batch; to use parallel with ncpus argument to slurm, use rep_workers instead to make sure that the repititions are paralleized

list(
    ## Find true value
    tar_target(true_value, {
        d_interventional <- simple_sim(1000000, intervention = 1)
        mean(d_interventional$Y_1)
    }),
    ## Run simulation study
    tar_rep(
        name = simulation_study_simple,
        command = {
            d <- simple_sim(sample_size)
            run_rtmle_setup(data = d,
                            learner = "my_super_learner",
                            eta = etas,
                            verbose = FALSE,
                            b = b_max)
        },
        batches = batches,
        reps = reps
    ),
    ## Find coverage
    tar_target(
        coverage,
        {
            simulation_study_simple[Protocol == "Always_A", .(coverage = mean((true_value >= Lower) & (true_value <= Upper)), width = mean(Upper - Lower)
                                                              ), by = .(subsample_size_percentage, method, B)]
        }
    ),
    ## Make plots
    tar_target(simple_plot_632, plot_bootstrap(coverage, fix_eta_plot = 0.632)),
    tar_target(simple_plot_5, plot_bootstrap(coverage, fix_eta_plot = 0.5)),
    tar_target(simple_plot_8, plot_bootstrap(coverage, fix_eta_plot = 0.8)),
    tar_target(simple_plot_9, plot_bootstrap(coverage, fix_eta_plot = 0.9))
)

######################################################################
### _targets.R ends here
