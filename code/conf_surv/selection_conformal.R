this_file <- parent.frame(2)$ofile  # works when sourced via `source()`
this_dir <- dirname(this_file)

screening_conformal <- function(data.test, data.cal, weights.cal, surv_model,
                                screening_time, screening_prob, screening_crit,
                                lambda.scores=NULL, data.tune=NULL, weights.tune=NULL, refit_model=FALSE,
                                p.sel.accept = 0.5) {

    ## Tune the lambda parameter using tuning data (if provided)
    if(is.null(lambda.scores) && !is.null(data.tune) && !is.null(weights.tune)) {
        lambda_grid <- unique(c(0,sort(quantile(data.cal$time, seq(0.5, 0.95, length.out=50)))))
        lambda.tuning <- .conformal_autotune(data.tune, weights.tune, surv_model$clone(),
                                             screening_time, screening_prob, screening_crit,
                                             refit_model=refit_model, lambda_grid=lambda_grid, criterion="fisher", n_repeats=5)
        lambda.scores <- lambda.tuning$lambda_opt

        ##cat(sprintf("best lambda: %.2f.\n", lambda.scores))
    } else if (is.null(lambda.scores)) {
        lambda.scores <- 0
    }

    if(screening_crit == "low risk") {
        alpha <- 1 - screening_prob
        score_type <- "survival"
    } else {
        alpha <- screening_prob
        score_type <- "one_minus_survival"
    }

    ## Compute conformity scores for calibration data
    idx.cal.events <- which(data.cal$status==1)
    scores.cal <- compute_scores_from_model(data.cal[idx.cal.events,], surv_model,
                                            times = data.cal$time[idx.cal.events] + lambda.scores,
                                            score_type = score_type)
    if(screening_crit == "low risk") {
        scores.cal[data.cal$time[idx.cal.events] > screening_time] <- 0
    } else {
        scores.cal[data.cal$time[idx.cal.events] < screening_time] <- 0
    }

    ## Compute conformity scores for test data
    scores.test <- compute_scores_from_model(data.test, surv_model, times = screening_time + lambda.scores, score_type = score_type)

    ## Pre-compute denominator for conformal p-values
    idx.event.cal <- which(data.cal$status==1)
    den <- 1+nrow(data.cal)

    ## Initialize matrix for p-values
    n.test <- nrow(data.test)
    time_points <- c(screening_time)

    ## Compute score comparisons efficiently using outer()
    score_comparison <- outer(scores.cal, scores.test, FUN = match.fun(">="))

    ## Compute the numerator using matrix operations
    num_values <- 1 + colSums(weights.cal[idx.event.cal] * score_comparison, na.rm = TRUE)
    ## Compute and store p-values for the current horizon
    pvals <- pmin(1, num_values / den)

    ## Apply BH
    pvals.adj <- p.adjust(pvals, method = "BH")
    sel <- which(pvals.adj<=alpha)

    ## Check stability of selections. Remove if too unstable
    p.sel.boot <- mean(sapply(1:1000, function(b) {
        sum(p.adjust(sample(pvals, replace=T), method = "BH")<0.1)>0
    }))
    if(p.sel.boot<p.sel.accept && length(sel)>0) {
        sel.stable <- c()
    } else {
        sel.stable <- sel
    }
    return(list(selections = sel, selections.stable = sel.stable,
                p.values = pvals, p.values.adj = pvals.adj, lambda.scores = lambda.scores,
                p.sel.boot=p.sel.boot))
}

## ------------------------------------------------------------
## FAST conformal autotune: precompute curves, interpolate for λ
## ------------------------------------------------------------
## data.tune, weights.tune : tuning frame + IPC weights (same rows)
## surv_model              : trainable model with $predict(data, times)$predictions -> n x |times|
## screening_time/_prob/_crit : passed-through configuration
## lambda_grid             : vector of λ (if NULL → data-driven default)
## n_repeats               : number of random splits
## criterion               : "bh" = #discoveries after BH; "fisher" = -2*sum log p
## refit_model             : refit model on train split each repeat?
## seed                    : RNG seed
## grid_size               : #time points for precomputed grid (trade speed/accuracy)
##
## Returns: list(lambda_opt, summary, per_repeat)
.conformal_autotune <- function(
  data.tune, weights.tune,
  surv_model,
  screening_time, screening_prob, screening_crit,
  lambda_grid = NULL,
  n_repeats = 5,
  criterion = c("bh","fisher"),
  refit_model = FALSE,
  seed = NULL,
  grid_size = 400
) {
  criterion <- match.arg(criterion)
  if (!is.null(seed)) set.seed(seed)

  # default λ grid (upper tail of observed times)
  if (is.null(lambda_grid)) {
    lambda_grid <- sort(unique(stats::quantile(data.tune$time, seq(0, 1, length.out = 200))))
  }
  G <- length(lambda_grid)
  N <- nrow(data.tune)

  # split proportions
  split_fracs <- if (refit_model) c(train=.7, cal=.15, test=.15) else c(train=0, cal=.5, test=.5)
  stopifnot(length(split_fracs)==3L, abs(sum(split_fracs)-1) < 1e-8)

  split_once <- function(N, fracs) {
    idx <- sample.int(N)
    n_tr <- max(0L, floor(fracs["train"] * N))
    n_ca <- max(1L, floor(fracs["cal"]   * N))
    tr <- if (n_tr>0) idx[seq_len(n_tr)] else integer(0)
    ca <- idx[seq_len(n_ca) + n_tr]
    te <- setdiff(idx, c(tr, ca))
    list(train=tr, cal=ca, test=te)
  }

  # gating rule (matches your earlier code)
  .gate_cal <- function(scores_cal_ev, time_ev, t_screen, crit) {
    out <- scores_cal_ev
    if (crit == "low risk")  out[ time_ev >  t_screen] <- 0
    if (crit == "high risk") out[ time_ev <  t_screen] <- 0
    out
  }

  # per-repeat stats
  r_stats <- vector("list", n_repeats)

  for (rep in seq_len(n_repeats)) {
    sp <- split_once(N, split_fracs)

    # (re)fit model
    if (refit_model && length(sp$train)>0) {
      surv_model$fit(stats::as.formula(Surv(time, status) ~ .), data = data.tune[sp$train, , drop=FALSE])
    }

    # split data/weights
    data_cal   <- data.tune[sp$cal,  , drop=FALSE]
    data_test  <- data.tune[sp$test, , drop=FALSE]
    w_cal_full <- weights.tune[sp$cal]
    time_cal   <- data_cal$time
    status_cal <- data_cal$status
    idx_ev_cal <- which(status_cal == 1L)
    n_cal      <- nrow(data_cal)
    n_test     <- nrow(data_test)

    # If no events in cal or no test rows → skip safely
    if (length(idx_ev_cal) == 0L || n_test == 0L) {
      r_stats[[rep]] <- data.frame(rep=rep, lambda=lambda_grid, stat=0)
      next
    }
    time_ev <- time_cal[idx_ev_cal]
    w_cal_ev <- w_cal_full[idx_ev_cal]

    # ----- time grid (single) -----
    # cover [0, max target time] where targets are time_ev + max(λ) and screening_time + max(λ)
    t_max <- max(0, max(time_ev, na.rm=TRUE) + max(lambda_grid, 0), screening_time + max(lambda_grid, 0))
    t_grid <- seq(0, t_max, length.out = grid_size)

    # ----- precompute survival curves once per split -----
    # cal EVENTS only (we only need event rows for cal)
    surv_cal_ev <- surv_model$predict(data_cal[idx_ev_cal, , drop=FALSE], t_grid)$predictions  # n_ev x |grid|
    # test (all rows)
    surv_test   <- surv_model$predict(data_test, t_grid)$predictions                           # n_test x |grid|

    # ----- build score matrices for all λ via interpolation -----
    # cal events: per-row vector xout = pmax(0, time_ev + λ)
    t_out_cal <- outer(time_ev, lambda_grid, function(t, l) pmax(0, t + l))  # n_ev x G
    scores_cal_mat <- matrix(NA_real_, nrow = length(idx_ev_cal), ncol = G)
    for (i in seq_len(nrow(surv_cal_ev))) {
        scores_cal_mat[i, ] <- stats::approx(t_grid, surv_cal_ev[i, ], xout = t_out_cal[i, ], rule = 2)$y
    }
    # gate by screening_time (as in your code)
    scores_cal_mat <- apply(scores_cal_mat, 2, .gate_cal, time_ev = time_ev,
                            t_screen = screening_time, crit = screening_crit)

    # test: same xout for every row (screening_time + λ)
    t_out_test <- pmax(0, screening_time + lambda_grid)  # length G
    # vectorized per-row interpolation
    scores_test_mat <- matrix(NA_real_, nrow = n_test, ncol = G)
    for (i in seq_len(nrow(surv_test))) {
      scores_test_mat[i, ] <- stats::approx(t_grid, surv_test[i, ], xout = t_out_test, rule = 2)$y
    }

    # ----- compute conformal stats for all λ (no extra predicts) -----
    alpha <- if (screening_crit == "low risk") (1 - screening_prob) else screening_prob
    den   <- 1 + n_cal

    stats_rep <- numeric(G)
    for (g in seq_len(G)) {
      s_cal  <- scores_cal_mat[, g]            # length n_ev
      s_test <- scores_test_mat[, g]           # length n_test

      # outer compare (n_ev x n_test) → TRUE if cal >= test
      comp <- outer(s_cal, s_test, FUN = ">=")
      # numerator per test obs: 1 + Σ w_cal_ev * 1{ s_cal >= s_test_j }
      num  <- 1 + as.numeric(t(comp) %*% w_cal_ev)  # length n_test
      pval <- pmin(1, num / den)

      if (criterion == "bh") {
        stats_rep[g] <- sum(stats::p.adjust(pval, method = "BH") <= alpha)
      } else { # fisher
        pv <- pval[is.finite(pval) & pval > 0]
        stats_rep[g] <- if (length(pv)) -2 * sum(log(pv)) else 0
      }
    }

    r_stats[[rep]] <- data.frame(rep = rep, lambda = lambda_grid, stat = stats_rep, row.names = NULL)
  }

  per_repeat <- do.call(rbind, r_stats)
  # summarize across repeats
  agg <- aggregate(stat ~ lambda, data = per_repeat,
                   FUN = function(z) c(mean = mean(z), sd = stats::sd(z)))
  summary <- data.frame(
    lambda    = agg$lambda,
    mean_stat = agg$stat[,1],
    sd_stat   = agg$stat[,2],
    n         = n_repeats,
    row.names = NULL
  )
  #best_lambda <- summary$lambda[which.max(summary$mean_stat)]
  best_lambda <- select_lambda_1se(summary, direction = "maximize")$lambda_1se

  list(
    lambda_opt = best_lambda,
    summary    = summary,
    per_repeat = per_repeat
  )
}

## Select λ by the 1-SE rule (glmnet-style)
## - direction="maximize" (your case): pick largest λ with mean >= max_mean - SE_at_max
## - direction="minimize"             : pick largest λ with mean <= min_mean + SE_at_min
select_lambda_1se <- function(df,
                              lambda_col = "lambda",
                              mean_col   = "mean_stat",
                              sd_col     = "sd_stat",
                              n_col      = "n",
                              direction  = c("maximize","minimize"),
                              already_se = FALSE) {
    direction <- match.arg(direction)
    x <- df[complete.cases(df[, c(lambda_col, mean_col, sd_col, n_col)]), ]
    if (nrow(x) == 0) stop("No complete rows.")

    ## Standard error at each λ
    se <- if (already_se) x[[sd_col]] else x[[sd_col]] / sqrt(x[[n_col]])
    se[!is.finite(se)] <- NA_real_

    if (direction == "maximize") {
        i_best <- which.max(x[[mean_col]])
        thr <- x[[mean_col]][i_best] - se[i_best]
        idx <- which(x[[mean_col]] >= thr)
    } else {
        i_best <- which.min(x[[mean_col]])
        thr <- x[[mean_col]][i_best] + se[i_best]
        idx <- which(x[[mean_col]] <= thr)
    }

    ## Among candidates within 1-SE, choose the largest λ
    j <- idx[which.max(x[[lambda_col]][idx])]
    res <- list(
        lambda_1se = x[[lambda_col]][j],
        index_1se  = j,
        lambda_best = x[[lambda_col]][i_best],
        index_best  = i_best,
        threshold   = thr,
        mean_best   = x[[mean_col]][i_best],
        se_at_best  = se[i_best]
    )
    return(res)
}
