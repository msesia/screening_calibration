# Load required libraries
suppressMessages(library(tidyverse))
library(survival)
#library(confsurv)
source("../conf_surv/utils_diagnostics.R")
source("../conf_surv/utils_semi_synthetic_data.R")
source("../conf_surv/survival_wrappers.R")
source("../conf_surv/censoring_wrappers.R")

source("../conf_surv/utils_weights_scores.R")
source("../conf_surv/selection_conformal.R")
source("../conf_surv/selection_calibration.R")
source("../conf_surv/selection_utils.R")
source("../conf_surv/utils_evaluation.R")

## Source utility functions for data generation and analysis
source("../utils/utils_misc.R")

set.seed(1)

######################
## Input parameters ##
######################

## Flag to determine if input should be parsed from command line
parse_input <- TRUE

if(parse_input) {
    ## Reading command line arguments
    args <- commandArgs(trailingOnly = TRUE)

    ## Checking if the correct number of arguments is provided
    if (length(args) < 8) {
        stop("Insufficient arguments provided. Expected 9 arguments.")
    }

    ## Assigning command line arguments to variables
    setup <- args[1]
    real_data <- as.integer(args[2])
    gen_model_type <- args[3]
    surv_model_type <- args[4]
    num_samples_train <- as.integer(args[5])
    num_samples_cal <- as.integer(args[6])
    screening_time <- as.numeric(args[7])
    batch <- as.integer(args[8])

} else {
    setup <- "v1"
    real_data <- 1
    gen_model_type <- "grf"
    surv_model_type <- "grf"
    num_samples_train <- 2000
    num_samples_cal <- 1000
    screening_time <- 2
    batch <- 1
}

######################
## Fixed parameters ##
######################

## Censoring model type
cens_model_type <- "grf"

## Do not use weights (default: TRUE, use FALSE only for debugging)
use_weights <- TRUE

## Test sample size
num_samples_test <- 1000

## Number of repetitions
batch_size <- 20

## Number of boostrap samples for LTT method (slow if >0)
B_boot <- 0

## Use built-in model CI (slow, due to bootstrapping)
use_model_ci <- FALSE

## Fixed time points and screening probabilities
screening_time <- screening_time
screening_prob <- 0.9
screening_crit <- "low risk"

####################
## Prepare output ##
####################

## Store important parameters including model types
header <- tibble(real_data = real_data,
                 gen_model_type = gen_model_type,
                 surv_model_type = surv_model_type,
                 cens_model_type = cens_model_type,
                 n_train = num_samples_train,
                 n_cal = num_samples_cal,
                 n_test = num_samples_test,
                 batch = batch)

## Generate a unique and interpretable file name based on the input parameters
output_file <- paste0("results/", setup, "/",
                      "real_", real_data,
                      "_gen_", gen_model_type,
                      "_surv_", surv_model_type,
                      "_train", num_samples_train,
                      "_cal", num_samples_cal,
                      "_time", screening_time,
                      "_batch", batch, ".txt")

## Print the output file name to verify
cat("Output file name:", output_file, "\n")

###############################################################
## Load the raw data and initialize semi-synthetic generator ##
###############################################################

load_raw_data <- function() {
    combined <- readRDS("../../data/fhrd_data.rds")
    combined <- as_tibble(do.call("rbind", combined))
    combined <- combined %>%
        mutate(time=month, status=event) %>%
        select(-month, -event) %>%
        select(time, status, everything())
    colnames(combined) <- c("time", "status", paste("X", seq(1,ncol(combined)-2), sep=""))
    return(combined)
}
data.raw <- load_raw_data()

## Instantiate generator (trains models once)
if(real_data) {
    data.generator <- RealDataGenerator$new(data = data.raw)
} else {
    data.generator <- SemiSyntheticDataGenerator$new(data = data.raw, surv_model_type = gen_model_type, cens_model_type = gen_model_type)
}


if(FALSE) {
    ## For debugging: test the data generator
    ## Generate a single semi-synthetic dataset
    data.synthetic <- data.generator$sample(shuffle=TRUE)
    data_list <- list(
        "Original" = data.raw,
        "Synthetic" = data.synthetic
    )
    compare_data_signals(data_list)

}

split_data_n <- function(data, n_train, n_cal, n_test, shuffle = TRUE) {
  n_total <- nrow(data)
  stopifnot(n_train + n_cal + n_test <= n_total)

  if (shuffle) {
    idx <- sample(seq_len(n_total))
  } else {
    idx <- seq_len(n_total)
  }

  idx_train <- idx[1:n_train]
  idx_cal <- idx[(n_train + 1):(n_train + n_cal)]
  idx_test <- idx[(n_train + n_cal + 1):(n_train + n_cal + n_test)]

  return(list(
    data.train = data[idx_train, , drop = FALSE],
    data.cal   = data[idx_cal,   , drop = FALSE],
    data.test  = data[idx_test,  , drop = FALSE]
  ))
}

###################################################
## Instantiate the survival and censoring models ##
###################################################

surv_model <- init_surv_model(surv_model_type)
surv_model_tune <- init_surv_model(surv_model_type)
surv_model_large <- init_surv_model(surv_model_type)

# List of covariates to use for censoring model
num_features <- ncol(data.raw)-2
use.covariates <- paste("X", 1:min(num_features, num_features), sep="")

# Instantiate censoring model based on the specified type
cens_base_model <- init_censoring_model(cens_model_type, use_covariates=use.covariates)
cens_model <- CensoringModel$new(model = cens_base_model)

#######################################
# Define function to analyze the data #
#######################################

analyze_data <- function(data.train, data.cal, data.test, surv_model, cens_model, data.test.oracle) {
    ## Fit the Kaplan-Meier survival model
    surv_object <- Surv(time = data.cal$time, event = data.cal$status)
    km_fit <- survival::survfit(surv_object ~ 1)

    ## Fit the survival model on the training data
    surv_model$fit(Surv(time, status) ~ ., data = data.train)

    ## Fit the survival model on all training and calibration data
    data.supervised <- rbind(data.train, data.cal)
    surv_model_large$fit(Surv(time, status) ~ ., data = data.supervised)
   
    ## Predict with large survival model
    model_pred <- as.numeric(surv_model_large$predict(data.test, screening_time)$predictions)
    sel_model <- select_patients_threshold(model_pred, screening_prob, screening_crit)
    sel_model_avg <- select_patients_average(model_pred, screening_prob, screening_crit)
   
    ## Predict with oracle model (if available)
    if(!real_data) {
        oracle_pred <- as.numeric(data.generator$surv_model$predict(data.test, screening_time)$predictions)

    } else {
        oracle_pred <- NULL
    }

    ## Fit the censoring model on a subset of the training data
    cens_model$fit(data = data.train)

    ## Selections with oracle model (if available)
    if(!is.null(oracle_pred)) {
        sel_oracle <- select_patients_threshold(oracle_pred, screening_prob, screening_crit)
        sel_oracle_group <- select_patients_average(oracle_pred, screening_prob, screening_crit)
    } else {
        sel_oracle <- NULL
        sel_oracle_group <- NULL
    }

    ## Selection with black-box model
    sel_model <- select_patients_threshold(model_pred, screening_prob, screening_crit)
    sel_model_group <- select_patients_average(model_pred, screening_prob, screening_crit)

    if(screening_crit == "low risk") {
        alpha <- 1 - screening_prob
        score_type <- "survival"
    } else {
        alpha <- screening_prob
        score_type <- "one_minus_survival"
    }

    ## Compute conformity scores for calibration data
    scores.cal <- compute_scores_from_model(data.cal, surv_model, times = screening_time, score_type = score_type)

    ## Compute conformity scores for test data
    scores.test <- compute_scores_from_model(data.test, surv_model, times = screening_time, score_type = score_type)

    ## Compute censoring weights at event times
    weights.cal.event <- compute_ipcw_weights(data.cal, cens_model, data.cal$time, data.cal$status, ipcw_method="et")
    ## Compute censoring weights at fixed times
    weights.cal.fixed <- compute_ipcw_weights(data.cal, cens_model, data.cal$time, data.cal$status, ipcw_method="ft", t0=screening_time)

    ## Compute censoring weights for training data (for tuning conformal method)
    weights.train <- compute_ipcw_weights(data.train, cens_model, data.train$time, data.train$status, ipcw_method="et")

    ## Selections with conformal p-values
    res_conformal <- screening_conformal(data.test, data.cal, weights.cal.event, surv_model,
                                         screening_time, screening_prob, screening_crit,
                                         data.tune=data.train, weights.tune=weights.train, refit_model=FALSE,
                                         p.sel.accept = 0.9)
    sel_conformal <- res_conformal$selections
    sel_conformal_stable <- res_conformal$selections.stable

    
    ## Selections with calibration methods
    selections_calibration <- list()

    delta_grid <- c(0.05, 0.10, 0.20, 0.50)  ## loop over these
    ipcw_grid  <- c("et", "ft")

    for (delta in delta_grid) {
        for (ipcw_method in ipcw_grid) {
            ## 1) Calibration table
            table.cal <- compute_calibration_table(
                scores.cal, data.cal$time, data.cal$status, screening_time,
                screening_crit = screening_crit,
                ipcw_method    = ipcw_method,
                weights_event  = weights.cal.event,
                weights_fixed  = weights.cal.fixed,
                num_lambda     = 100,
                delta          = delta,
                include_uniform   = TRUE,
                include_bootstrap = FALSE,
                B_boot            = 1000,
                min_n_asymptotics = 10
            )
            key_base <- sprintf("%s|%.2f", ipcw_method, delta)
            
            ## 2a) Greedy on delta-method pointwise estimate
            lam_pt <- select_lambda_greedy(
                lambda = table.cal$lambda,
                risk   = table.cal$risk_ub_point_delta,
                alpha  = alpha,
                size   = table.cal$n_selected
            )$lambda
            selections_calibration[[sprintf("%s|pt_delta", key_base)]] <- select_patients_threshold(scores.test, lam_pt, screening_crit)

            ## 2b) Greedy on uniform band
            lam_uni <- select_lambda_greedy(
                lambda = table.cal$lambda,
                risk   = table.cal$risk_ub_uniform,
                alpha  = alpha,
                size   = table.cal$n_selected
            )$lambda
            selections_calibration[[sprintf("%s|uniform", key_base)]] <- select_patients_threshold(scores.test, lam_uni, screening_crit)

            ## 2c) LTT from anchor on delta-method band
            lam_ltt <- select_lambda_LTT_from_anchors(
                lambda  = table.cal$lambda,
                risk    = table.cal$risk_ub_point_delta,
                alpha   = alpha,                         ## = 1 - screening_prob
                anchors = c(1-alpha, 1-alpha/2),
                size    = table.cal$n_selected
            )$lambda
            selections_calibration[[sprintf("%s|ltt_delta", key_base)]] <- select_patients_threshold(scores.test, lam_ltt, screening_crit)
        }
    }

    ## Combine all selections
    selections <- list("oracle (point)" = sel_oracle,
                       "oracle (group)" = sel_oracle_group,
                       "model (point)" = sel_model,
                       "model (group)" = sel_model_group,
                       "conformal" = sel_conformal,
                       "conformal (stable)" = sel_conformal_stable)
    selections <- c(selections, selections_calibration)


    ## Evaluate and format results
    evaluated.no.oracle <- map2_dfr(selections, names(selections), function(selected, method_name) {
        res.raw <- evaluate_selections_without_oracle(data.test, selected, screening_time, screening_prob, screening_crit)
        res.raw %>%
            transmute(
                Method = method_name,
                Time = Screening.time,
                Criterion = screening_crit,
                Probability = screening_prob,
                Screened = Screened,
                Survival_lower = Proportion.survived.lower,
                Survival_upper = Proportion.survived.upper
            )
    })
    evaluated.oracle <- map2_dfr(selections, names(selections), function(selected, method_name) {
        res.raw <- evaluate_selections(data.test.oracle, selected, screening_time, screening_prob, screening_crit)
        res.raw %>%
            transmute(
                Method = method_name,
                Time = Screening.time,
                Criterion = screening_crit,
                Probability = screening_prob,
                Screened = Screened,
                Survival = Proportion.survived
            )
    }) %>% as_tibble()
    if(real_data) {
        evaluated.oracle$Survival <- NA
    }
    evaluated <- inner_join(evaluated.oracle, evaluated.no.oracle, by = join_by(Method, Time, Criterion, Probability, Screened))
    print(evaluated)
    return(evaluated)
}

#######################################
## Define function to run experiment ##
#######################################

run_experiment <- function(random.state) {
    set.seed(random.state)

    ## Generate training, calibration, and test data
    data.synthetic.oracle <- data.generator$sample(shuffle=TRUE, return_oracle=TRUE)$oracle
    splits <- split_data_n(data.synthetic.oracle, n_train = num_samples_train, n_cal = num_samples_cal, n_test = num_samples_test)
    data.train.oracle <- splits$data.train
    data.cal.oracle <- splits$data.cal
    data.test.oracle <- splits$data.test

    ## Remove true event and censoring times from the data (right-censoring)
    data.train <- data.train.oracle |> select(-event_time, -censoring_time)
    data.cal <- data.cal.oracle |> select(-event_time, -censoring_time)
    data.test <- data.test.oracle |> select(-event_time, -censoring_time)

    ## Run analysis
    results <- analyze_data(data.train, data.cal, data.test, surv_model, cens_model, data.test.oracle)

    return(results)
}


## Function to run multiple experiments and gather results
## Args:
##   batch_size: Number of repetitions for each experimental setting
## Returns:
##   A tibble containing the combined results of all experiments
run_multiple_experiments <- function(batch_size) {
    results_df <- data.frame()  # Initialize an empty data frame to store cumulative results

    # Print a progress bar header
    cat("Running experiments\n")
    pb <- txtProgressBar(min = 0, max = batch_size, style = 3)  # Initialize progress bar

    # Loop over each repetition
    for (i in 1:batch_size) {
        random.state <- batch*1000 + i
        res <- run_experiment(random.state)  # Run experiment and get the result
        
        ## Combine the results with experiment metadata
        result_df <- tibble(Seed = random.state) |> cbind(header) |> cbind(res)

        # Add the result to the cumulative data frame
        results_df <- rbind(results_df, result_df) %>% as_tibble()

        # Write the cumulative results to the CSV file
        write.csv(results_df, output_file, row.names = FALSE)

        setTxtProgressBar(pb, i)  # Update progress bar
    }

    close(pb)  # Close the progress bar

    return(results_df)  # Return the cumulative results data frame
}

#####################
## Run experiments ##
#####################

## Run the experiments with specified parameters
results <- run_multiple_experiments(batch_size)

#browser()
##results |> filter(Time==5, Probability==0.8, Criterion=="low risk") |> select(Seed, Criterion, Method, Screened, Survival)

##results |> filter(Time==5, Probability==0.8, Criterion=="low risk") |> select(Seed, Criterion, Method, Screened, Survival) |>
##group_by(Criterion, Method) |> summarise(Screened=mean(Screened), Survival=mean(Survival))
