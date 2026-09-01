# Load required libraries
suppressMessages(library(tidyverse))
library(survival)
#library(confsurv)
source("../conf_surv/utils_diagnostics.R")
source("../conf_surv/utils_semi_synthetic_data.R")
source("../conf_surv/utils_tilted_data.R")
source("../conf_surv/survival_wrappers.R")
source("../conf_surv/censoring_wrappers.R")

source("../conf_surv/utils_weights_scores.R")
source("../conf_surv/pseudo_outcomes.R")
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
        stop("Insufficient arguments provided. Expected 8 or 9 arguments.")
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
    cens_mix <- if (length(args) >= 9) as.numeric(args[9]) else 0
} else {
    setup <- "v0c"
    real_data <- 0
    gen_model_type <- "grf"
    surv_model_type <- "grf"
    num_samples_train <- 5000
    num_samples_cal <- 1000
    screening_time <- 2
    batch <- 1
    cens_mix <- 0
}

######################
## Fixed parameters ##
######################

## Censoring model type
cens_model_type <- "grf"

## Test sample size
num_samples_test <- 1000

## Number of repetitions
batch_size <- 20

## Number of boostrap samples for LTT method (slow if >0)
B_boot <- 0

# Smallest selected sample size for asymptotics
min_n_sel <- 15

## Use built-in model CI (slow, due to bootstrapping)
use_model_ci <- FALSE

## Tune the conformal time shift gamma separately for each pseudo-outcome?
##
## gamma only sets the time shift of the nonconformity score, and that score
## is the SAME function of x for both flavours -- the two arms rank patients
## identically and differ only in how H(z) is estimated. So one ET-based
## tuning is shared by default.
##
## TRUE tunes each arm on its own p-values, which is arguably more faithful
## to Algorithm 1, but the AIPCW tuner costs ~10 s per replicate against
## ~0.2 s for the shared one, i.e. ~3 extra minutes per 20-replicate batch.
tune_per_flavor <- FALSE

## Fixed time shifts for the conformal method, in ADDITION to the tuned arm.
## Under v0c this sweeps gamma on a hard-coded grid, so its effect can be read
## off directly rather than inferred from the tuner. Only the "dr" arm is
## swept, since gamma acts on the score and not on the pseudo-outcome. Empty
## for the other setups, which keep the columns they already have.
gamma_grid_fixed <- if (setup == "v0c") c(0, 1, 10, 100, 1000) else numeric(0)

## Fixed time points and screening probabilities
screening_time <- screening_time
screening_prob <- 0.9
screening_crit <- "low risk"

###############################################################
## Load the raw data and initialize semi-synthetic generator ##
###############################################################

## load_raw_data_old <- function() {
##     e <- new.env()
##     load("../../data/oncology_data.RData", envir = e)
##     ## Combine all data
##     combined <- rbind(
##         cbind(e$data.train, .source = "train"),
##         cbind(e$data.cal,   .source = "cal"),
##         cbind(e$data.test,  .source = "test")
##     ) %>%
##         as_tibble() %>%
##         select(-.source)
##     return(combined)
## }

load_raw_data <- function() {
    combined <- readRDS("../../data/sample_data_full.rds")
    combined <- as_tibble(do.call("rbind", combined))
    combined <- combined %>%
        mutate(time=month, status=event) %>%
        select(-month, -event) %>%
        select(time, status, everything())
    colnames(combined) <- c("time", "status", paste("X", seq(1,ncol(combined)-2), sep=""))
    return(combined)
}
data.raw <- load_raw_data()

## ---- setup "v0t": a deliberately misspecified censoring model ---------
q.z <- as.numeric(quantile(data.raw$X32, 0.5))
tilt_surv_analysis <- list(f = ~ -1 + 4 * as.numeric(X32 > q.z))
tilt_cens_analysis <- list(f = ~ -1 + 6 * as.numeric(X32 > q.z))
tilt_cens_alt      <- list(f = ~ -1 + 0 * as.numeric(X32 > q.z))

## Instantiate generators (each trains its models once).
##
## Outside v0t the analysis data is ordinary semi-synthetic data and the
## censoring model is fitted on the training data itself, as before.
gen.cens.alt <- NULL
if (real_data) {
    data.generator <- RealDataGenerator$new(data = data.raw)
} else if (setup == "v0t") {
    ## analysis distribution: supplies train / calibration / test
    data.generator <- TiltedSemiSyntheticDataGenerator$new(
        data = data.raw, surv_model_type = gen_model_type, cens_model_type = gen_model_type,
        tilt_surv = tilt_surv_analysis, tilt_cens = tilt_cens_analysis)
    ## alternative distribution: same survival tilt, different censoring.
    ## Only ever used to corrupt the censoring model's training sample.
    gen.cens.alt <- TiltedSemiSyntheticDataGenerator$new(
        data = data.raw, surv_model_type = gen_model_type, cens_model_type = gen_model_type,
        tilt_surv = tilt_surv_analysis, tilt_cens = tilt_cens_alt)
    cat(sprintf("Setup v0t: censoring model fitted on %.0f%% analysis law + %.0f%% alternative law%s\n",
                100*(1 - cens_mix), 100*cens_mix,
                if (cens_mix == 0) "  (-> Ghat correctly specified)" else "  (-> Ghat misspecified)"))
    compare_generators(data.generator, gen.cens.alt, screening_time,
                       labels = c("ANALYSIS law", "ALTERNATIVE law"))
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
use.covariates <- paste("X", seq_len(num_features), sep="")

# Instantiate censoring model based on the specified type
cens_base_model <- init_censoring_model(cens_model_type, use_covariates=use.covariates)
cens_model <- CensoringModel$new(model = cens_base_model)

####################
## Prepare output ##
####################

## Store important parameters including model types
header <- tibble(real_data = real_data,
                 gen_model_type = gen_model_type,
                 surv_model_type = surv_model_type,
                 cens_model_type = cens_model_type,
                 n_train = num_samples_train,
                 cens_mix = cens_mix,
                 n_cal = num_samples_cal,
                 n_test = num_samples_test,
                 batch = batch)

## Tag every v0t file with its mixture, including 0, so the control point
## has a file of its own. Other setups keep the names they have.
## Must match TC_TAG in submit_experiment_1.sh.
traincens_tag <- if (setup == "v0t") paste0("_mix", cens_mix) else ""

## Generate a unique and interpretable file name based on the input parameters
output_file <- paste0("results/", setup, "/",
                      "real_", real_data,
                      "_gen_", gen_model_type,
                      "_surv_", surv_model_type,
                      "_train", num_samples_train,
                      traincens_tag,
                      "_cal", num_samples_cal,
                      "_time", screening_time,
                      "_batch", batch, ".txt")

## Make sure the output directory exists. The submit script also does this,
## but the R script must not depend on that: running a NEW setup (or running
## interactively) otherwise fails at write.csv() with "cannot open the
## connection" only after the first replicate has already been computed.
dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
if (!dir.exists(dirname(output_file))) {
    stop("could not create output directory: ", dirname(output_file))
}

## Print the output file name to verify
cat("Output file name:", output_file, "\n")


#######################################
# Define function to analyze the data #
#######################################

analyze_data <- function(data.train, data.cal, data.test, surv_model, cens_model, data.test.oracle,
                         data.train.cens = data.train) {
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

    ## Fit the censoring model. Under v0t with a non-zero tilt this is an
    ## auxiliary sample whose censoring mechanism differs from the analysis
    ## data -- which is exactly what makes Ghat misspecified. Otherwise it
    ## is the training data itself.
    cens_model$fit(data = data.train.cens)

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

    ## Compute censoring weights for training data (for tuning conformal method)
    weights.train <- compute_ipcw_weights(data.train, cens_model, data.train$time, data.train$status, ipcw_method="et")

    ## Pseudo-outcomes: IPCW (et) and AIPCW (dr). These share the same IPCW
    ## weights, so they differ only in the augmentation terms.
    po <- list(
        et = build_pseudo_outcomes(data.cal, cens_model, data.cal$time, data.cal$status,
                                   screening_time, flavor="et",
                                   weights_event=weights.cal.event),
        dr = build_pseudo_outcomes(data.cal, cens_model, data.cal$time, data.cal$status,
                                   screening_time, flavor="dr", surv_model=surv_model,
                                   weights_event=weights.cal.event)
    )

    ## Tune the time shift on the training data (see tune_per_flavor above).
    gamma.opt <- if (tune_per_flavor) {
        sapply(c("et","dr"), function(fl)
            .conformal_autotune(data.train, weights.train, surv_model, cens_model,
                                screening_time, screening_prob, flavor = fl)$gamma_opt)
    } else {
        g <- .conformal_autotune(data.train, weights.train, surv_model, cens_model,
                                 screening_time, screening_prob, flavor = "et")$gamma_opt
        c(et = g, dr = g)
    }

    ## Selections with conformal p-values
    res_conformal <- screening_conformal(data.test, data.cal, po$et, surv_model,
                                         screening_time, screening_prob, screening_crit,
                                         gamma = gamma.opt[["et"]], p.sel.accept = 0.9)
    sel_conformal <- res_conformal$selections
    sel_conformal_stable <- res_conformal$selections.stable

    res_conformal_dr <- screening_conformal(data.test, data.cal, po$dr, surv_model,
                                            screening_time, screening_prob, screening_crit,
                                            gamma = gamma.opt[["dr"]], p.sel.accept = 0.9)
    sel_conformal_dr <- res_conformal_dr$selections
    sel_conformal_dr_stable <- res_conformal_dr$selections.stable

    ## Conformal at each FIXED gamma on the grid, alongside the tuned arm
    ## above. Keys are "conformal-dr-g<gamma>"; gamma = 0 reproduces what was
    ## previously reported as "conformal-dr-notune".
    selections_gamma <- list()
    for (g in gamma_grid_fixed) {
        rg <- screening_conformal(data.test, data.cal, po$dr, surv_model,
                                  screening_time, screening_prob, screening_crit,
                                  gamma = g, p.sel.accept = 0.9)
        key <- sprintf("conformal-dr-g%g", g)
        selections_gamma[[key]] <- rg$selections
        selections_gamma[[sprintf("%s (stable)", key)]] <- rg$selections.stable
    }

    
    ## Selections with calibration methods
    selections_calibration <- list()

    delta_grid <- c(0.05, 0.10, 0.20, 0.50)  ## loop over these
    flavor_grid <- c("et", "dr")             ## IPCW and AIPCW pseudo-outcomes

    for (delta in delta_grid) {
        for (flavor in flavor_grid) {
            ## 1) Calibration table
            table.cal <- compute_calibration_table(
                scores.cal, po[[flavor]],
                screening_crit = screening_crit,
                num_lambda     = 100,
                delta          = delta,
                include_uniform   = TRUE,
                include_bootstrap = FALSE,
                B_boot            = B_boot,
                min_n_sel = min_n_sel
            )
            ## 1b) second table at the per-path level delta/2 (Algorithm S4, step 3)
            table.cal.div2 <- compute_calibration_table(
                scores.cal, po[[flavor]],
                screening_crit = screening_crit,
                num_lambda     = 100,
                delta          = delta/2,
                include_uniform   = FALSE,
                include_bootstrap = FALSE,
                min_n_sel = min_n_sel
            )
            key_base <- sprintf("%s|%.2f", flavor, delta)
            
            ## 2a) Greedy on delta-method pointwise estimate
            lam_pt <- select_lambda_greedy(
                lambda = table.cal$lambda,
                risk   = table.cal$risk_ub_point_delta,
                alpha  = alpha
            )$lambda
            selections_calibration[[sprintf("%s|pt_delta", key_base)]] <- select_patients_threshold(scores.test, lam_pt, screening_crit)

            ## 2b) Greedy on uniform band
            lam_uni <- select_lambda_greedy(
                lambda = table.cal$lambda,
                risk   = table.cal$risk_ub_uniform,
                alpha  = alpha
            )$lambda
            selections_calibration[[sprintf("%s|uniform", key_base)]] <- select_patients_threshold(scores.test, lam_uni, screening_crit)

            ## 2c) LTT on delta-method band, using two anchors
            lam_ltt <- select_lambda_LTT_from_anchors(
                lambda  = table.cal.div2$lambda,
                risk    = table.cal.div2$risk_ub_point_delta,
                alpha   = alpha,
                anchors = c(1-alpha, 1-alpha/2)
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
                       "conformal (stable)" = sel_conformal_stable,
                       "conformal-dr" = sel_conformal_dr,
                       "conformal-dr (stable)" = sel_conformal_dr_stable,
                       "conformal-dr (stable)" = sel_conformal_dr_stable)
    selections <- c(selections, selections_gamma, selections_calibration)


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

    ## Sample the censoring model is fitted on. Under v0t it is a mixture:
    ## a fraction cens_mix of the rows are drawn from the ALTERNATIVE
    ## censoring law and the rest from the analysis law. Both components
    ## carry the same survival tilt, so only the censoring differs.
    data.train.cens <- if (is.null(gen.cens.alt)) {
        data.train
    } else {
        n_alt <- round(cens_mix * num_samples_train)
        n_ana <- num_samples_train - n_alt
        parts <- list()
        if (n_ana > 0) {
            a <- data.generator$sample(shuffle = TRUE)
            parts[[length(parts) + 1L]] <- a[seq_len(min(n_ana, nrow(a))), , drop = FALSE]
        }
        if (n_alt > 0) {
            b <- gen.cens.alt$sample(shuffle = TRUE)
            parts[[length(parts) + 1L]] <- b[seq_len(min(n_alt, nrow(b))), , drop = FALSE]
        }
        do.call(rbind, parts)
    }

    ## Run analysis
    results <- analyze_data(data.train, data.cal, data.test, surv_model, cens_model, data.test.oracle,
                            data.train.cens = data.train.cens)

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
                                        
