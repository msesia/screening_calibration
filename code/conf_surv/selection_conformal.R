## ============================================================
## Conformal screening with FDR control (paper Section 4).
##
## Flavor-agnostic: the IPCW p-values (19)/(22) and the augmented ones
## (24) are the SAME computation applied to different pseudo-outcomes.
##
## Because R(x,t) = Shat(t + gamma | x) I(t <= t0) is non-increasing in t,
## {R(X,T) >= z} = {T <= b_z(X)} with b_z of (23), so
##
##   H(z) = P[R(X,T) >= z] = P[T <= b_z(X)]
##
## has the same form as theta(lambda) with t0 replaced by the
## covariate-dependent horizon b_z(X).  Hence
##
##   Hhat(z) = n^{-1} sum_i Psi_i( b_z(X_i) )
##
## for ANY flavor: with "et" this is exactly the weighted rank of (19),
## with "dr" it is the augmented estimator of Section 4.3.  The p-value is
##
##   phat(z) = Pi( (1 + n Hhat(z)) / (1 + n) ),
##   Pi f(z) = min{1, sup_{z' >= z} f(z')}                            (24)
##
## which is a no-op for "et" (already non-increasing in z) and the least
## non-increasing majorant for "dr".
##
## Depends on: pseudo_outcomes.R, utils_weights_scores.R.
## ============================================================


## ---------------------------------------------------------
## Survival curves evaluated at times.
## If `curves` (on `t_grid`) is supplied the values are interpolated rather
## than re-predicted -- the gamma tuner sweeps many time shifts and must
## not re-call the model for each one.
## `times` is either one common vector, or one time per row (row_times).
## ---------------------------------------------------------
## Linear interpolation of every row of `curves` (given on `t_grid`),
## vectorised.  Equivalent to approx(..., rule = 2) applied row by row, but
## without the per-row call: the tuner would otherwise make O(n |Gamma| R)
## calls to approx(), which dominated its running time.
.interp_common <- function(curves, t_grid, x) {          # one common time vector
    x <- pmin(pmax(x, t_grid[1L]), t_grid[length(t_grid)])          # rule = 2
    j <- findInterval(x, t_grid, all.inside = TRUE)
    w <- (x - t_grid[j]) / (t_grid[j + 1L] - t_grid[j])
    W <- matrix(w, nrow(curves), length(x), byrow = TRUE)
    curves[, j, drop = FALSE] * (1 - W) + curves[, j + 1L, drop = FALSE] * W
}

.interp_rowwise <- function(curves, t_grid, x) {         # one time per row
    x <- pmin(pmax(x, t_grid[1L]), t_grid[length(t_grid)])
    j <- findInterval(x, t_grid, all.inside = TRUE)
    w <- (x - t_grid[j]) / (t_grid[j + 1L] - t_grid[j])
    rows <- seq_len(nrow(curves))
    curves[cbind(rows, j)] * (1 - w) + curves[cbind(rows, j + 1L)] * w
}

.eval_curves <- function(data, surv_model, times, curves = NULL, t_grid = NULL, eps = 1e-12) {
    times <- pmax(0, as.numeric(times))
    S <- if (is.null(curves)) {
        matrix(as.numeric(surv_model$predict(data, times)$predictions), nrow = nrow(data))
    } else {
        .interp_common(curves, t_grid, times)
    }
    pmin(pmax(S, eps), 1)
}

.eval_curves_rowwise <- function(data, surv_model, row_times, curves = NULL, t_grid = NULL,
                                 eps = 1e-12) {
    row_times <- pmax(0, as.numeric(row_times))
    S <- if (is.null(curves)) {
        compute_scores_from_model(data, surv_model, times = row_times, score_type = "survival")
    } else {
        .interp_rowwise(curves, t_grid, row_times)
    }
    pmin(pmax(as.numeric(S), eps), 1)
}

## Shifted curves Shat(u_k + gamma | X_i) on the grid, monotonised so that
## the horizon search below is a binary search.
.shifted_curves <- function(data, surv_model, u, gamma = 0, curves = NULL, t_grid = NULL) {
    .row_apply(.eval_curves(data, surv_model, u + gamma, curves, t_grid), cummin)
}

## ---------------------------------------------------------
## b_z(X_i) of (23), as a 0-based grid index: the largest u_k with
## Shat(u_k + gamma | X_i) >= z, or -1 when no such point exists.
## Returns an n x length(z) integer matrix.
## ---------------------------------------------------------
.horizon_index <- function(Sshift, z, K) {
    z <- as.numeric(z)
    idx <- matrix(NA_integer_, nrow = nrow(Sshift), ncol = length(z))
    for (i in seq_len(nrow(Sshift))) {
        ## count of grid values strictly below z; rev() makes it non-decreasing
        idx[i, ] <- K - findInterval(z, rev(Sshift[i, ]), left.open = TRUE)
    }
    idx
}

## ---------------------------------------------------------
## Everything the p-value computation needs from the survival model at a
## given time shift.  WHICH pieces are needed depends on the pseudo-outcome,
## and this is the only place that distinction appears:
##
##   Psi^et is an INDICATOR in the horizon,
##       Psi_i(b_z(X_i)) = w_i I(T_i <= b_z(X_i)) = w_i I(z_cal_i >= z),
##   so each calibration subject is summarised by the single number
##       z_cal_i = Shat(T_i + gamma | X_i),  gated to -Inf when T_i > t0,
##   and no grid is involved at all.
##
##   Psi_dr instead has to be evaluated AT the horizon b_z(X_i), because
##   Qhat and the prefix sums depend on where it falls, so the whole shifted
##   curve is required.
##
## Both routes compute the same Hhat; test 7f checks they agree exactly.
## Using the indicator route where it applies is what keeps the gamma tuner
## affordable: O(n log n + m log n) per candidate instead of O(n(K + m)).
## ---------------------------------------------------------
conformal_inputs <- function(po, data.cal, data.test, surv_model, gamma = 0,
                             curves.cal = NULL, curves.test = NULL, t_grid = NULL,
                             force_horizons = FALSE) {
    z_test <- as.numeric(.eval_curves(data.test, surv_model, po$t0 + gamma,
                                      curves.test, t_grid))
    indicator <- (po$flavor == "et") && !force_horizons

    if (indicator) {
        z_cal <- .eval_curves_rowwise(data.cal, surv_model, po$time + gamma,
                                      curves.cal, t_grid)
        z_cal[po$time > po$t0] <- -Inf          # R(x,t) = 0 beyond the horizon
        list(z_cal = z_cal, z_test = z_test, Sshift = NULL)
    } else {
        list(z_cal = NULL, z_test = z_test,
             Sshift = .shifted_curves(data.cal, surv_model, po$u, gamma,
                                      curves.cal, t_grid))
    }
}

## Restrict the model-side inputs to a calibration / test split.  Both
## z_cal and Sshift are row-wise in the calibration subject, and z_test is
## row-wise in the cohort subject, so for a fixed gamma the inputs can be
## built once on a whole pool and then subset.
.subset_conformal_inputs <- function(inp, idx_cal, idx_test) {
    list(z_cal  = if (is.null(inp$z_cal))  NULL else inp$z_cal[idx_cal],
         z_test = inp$z_test[idx_test],
         Sshift = if (is.null(inp$Sshift)) NULL else inp$Sshift[idx_cal, , drop = FALSE])
}

## ---------------------------------------------------------
## Conformal p-values for a cohort.  Flavor-agnostic: swap po$flavor to
## switch between (19) and (24).
## ---------------------------------------------------------
conformal_pvalues <- function(po, inp) {
    n <- po$n
    z_test <- inp$z_test
    m <- length(z_test)
    if (!isTRUE(po$varying_horizon)) {
        stop("conformal screening needs pseudo-outcomes valid at varying horizons ",
             "(flavor 'et' or 'dr'); got '", po$flavor, "'.")
    }
   
    Hhat <- if (is.null(inp$Sshift)) {
        ## Indicator route: accumulate the weighted right tail by sorting.
        ord  <- order(inp$z_cal, decreasing = TRUE)
        cumw <- cumsum(po$w_event[ord])
        k    <- n - findInterval(z_test, sort(inp$z_cal), left.open = TRUE)  # #{z_cal >= z}
        ifelse(k > 0, cumw[pmax(k, 1L)], 0) / n
    } else {
        ## General route: evaluate Psi at the horizon b_z(X_i).  Done in
        ## column blocks -- one vectorised pass each, with the block width
        ## chosen so the temporaries stay near 1e6 entries.
        ELL <- .horizon_index(inp$Sshift, z_test, po$K)
        H <- numeric(m)
        width <- max(1L, min(m, floor(1e6 / max(n, 1))))
        for (start in seq(1L, m, by = width)) {
            j <- start:min(m, start + width - 1L)
            H[j] <- colMeans(pseudo_psi(po, ELL[, j, drop = FALSE]))
        }
        H
    }

    p <- (1 + n * Hhat) / (1 + n)
    ## Pi of (24): sup over z' >= z, then clip to [0,1].  The lower clip is
    ## not part of Pi but is needed because Hhat can be negative under
    ## augmentation; it preserves monotonicity and only enlarges p-values.
    ord <- order(z_test, decreasing = TRUE)
    g <- rep(NA_real_, m)
    g[ord] <- cummax(p[ord])
    list(p.values = pmin(1, pmax(0, g)), H.hat = Hhat)
}

## ---------------------------------------------------------
## screening_conformal()
##
##   po           : pseudo-outcome object built on data.cal with GRID
##                  ENDPOINT t0.  t0 is the end of the time grid, not the
##                  horizon at which Psi is evaluated: this routine
##                  evaluates Psi at the covariate-dependent horizons
##                  b_z(X_i) of (23), which vary with both the calibration
##                  subject and the candidate threshold z (typically a few
##                  hundred distinct horizons per cohort).  A single object
##                  suffices because R(x,t) = Shat(t+gamma|x) I[t <= t0] is
##                  zero beyond t0, so b_z(X) <= t0 always, and the object
##                  carries Ghat, Shat and the prefix sums across the whole
##                  grid [0, t0].
##                  flavor "et" -> eq. (19); flavor "dr" -> eq. (24).
##   gamma        : time shift of (17).  If NULL and tuning data are
##                  supplied it is tuned by .conformal_autotune().
##   p.sel.accept : nu of Section 4.5; abstain below it
##
## Only low-risk screening is supported.  For high-risk screening the
## score is non-DEcreasing in t and the relevant horizons live in
## [t0, Inf), so both the grid and the equivalence in (23) would have to
## be rebuilt; the function stops rather than returning something wrong.
## ---------------------------------------------------------
screening_conformal <- function(data.test, data.cal, po, surv_model,
                                screening_time, screening_prob,
                                screening_crit = "low risk",
                                gamma = NULL,
                                data.tune = NULL, weights.tune = NULL, cens_model = NULL,
                                p.sel.accept = 0.9, n.boot = 1000L) {

    if (screening_crit != "low risk") {
        stop("screening_conformal() supports screening_crit = 'low risk' only; ",
             "see the note in this file for why the high-risk case needs a different grid.")
    }
    stopifnot(nrow(data.cal) == po$n)
    ## the horizon comes from the pseudo-outcome object, so a mismatched
    ## screening_time would otherwise be silently ignored
    stopifnot(isTRUE(all.equal(as.numeric(screening_time), as.numeric(po$t0))))
    alpha <- 1 - screening_prob

    if (is.null(gamma)) {
        gamma <- if (!is.null(data.tune) && !is.null(weights.tune)) {
            .conformal_autotune(data.tune, weights.tune, surv_model, cens_model,
                                screening_time, screening_prob)$gamma_opt
        } else 0
    }

    inp <- conformal_inputs(po, data.cal, data.test, surv_model, gamma = gamma)
    res <- conformal_pvalues(po, inp)
    pvals <- res$p.values
    browser()

    pvals.adj <- p.adjust(pvals, method = "BH")
    sel <- which(pvals.adj <= alpha)

    ## Abstention (Section 4.5): bootstrap P(|S| > 0 | Dcal).  phat is a
    ## fixed function of Dcal, so resampling cohort rows is the same as
    ## resampling their p-values.
    p.sel.boot <- mean(vapply(seq_len(n.boot), function(b)
        any(p.adjust(sample(pvals, replace = TRUE), method = "BH") <= alpha), logical(1)))
    sel.stable <- if (p.sel.boot < p.sel.accept) integer(0) else sel

    list(selections = sel, selections.stable = sel.stable,
         p.values = pvals, p.values.adj = pvals.adj, H.hat = res$H.hat,
         gamma = gamma, p.sel.boot = p.sel.boot)
}

## ---------------------------------------------------------
## Cross-validation for the time shift gamma (Algorithm S5).
##
## Uses the same conformal_pvalues() core as the main routine.  Since the
## models are fixed, NOTHING here needs rebuilding per split: the
## pseudo-outcomes and the fitted curves are computed once for the whole
## tuning set, the model-side inputs once per gamma, and each repeat is a
## pure subsetting operation.
##
## The `flavor` argument selects which p-values gamma is tuned for; it
## should match the flavor the method will actually run with, since gamma is
## a power heuristic.  The guarantees hold conditionally on gamma however it
## was chosen, provided the tuning data are independent of Dcal.
## ---------------------------------------------------------
.conformal_autotune <- function(data.tune, weights.tune, surv_model, cens_model = NULL,
                                screening_time, screening_prob,
                                gamma_grid = NULL, n_repeats = 5,
                                criterion = c("fisher","bh"), flavor = "et",
                                max_cal = 1000L, max_test = 250L,
                                tune_grid = 100L, grid_size = 400,
                                seed = NULL, one_se = FALSE) {
    criterion <- match.arg(criterion)
    if (!is.null(seed)) set.seed(seed)
    alpha <- 1 - screening_prob
    N <- nrow(data.tune)

    if (is.null(gamma_grid)) {
        gamma_grid <- unique(c(0, as.numeric(stats::quantile(data.tune$time,
                                                             seq(0.5, 0.95, length.out = 50),
                                                             na.rm = TRUE))))
    }
    G <- length(gamma_grid)

    ## The models are fixed, so Psi_i depends only on subject i: build the
    ## pseudo-outcomes ONCE for the whole tuning set, and let each repeat
    ## subset them.  A coarse grid is enough here -- this chooses one scalar.
    po_full <- build_pseudo_outcomes(data.tune, cens_model, data.tune$time, data.tune$status,
                                     screening_time, flavor = flavor,
                                     surv_model = surv_model,
                                     weights_event = weights.tune,
                                     num_grid = tune_grid, max_grid = tune_grid)

    ## ... and one prediction pass for the whole tuning set; every gamma is
    ## then an interpolation.
    t_max <- max(data.tune$time, screening_time, na.rm = TRUE) + max(gamma_grid, 0)
    t_grid <- seq(0, t_max, length.out = grid_size)
    curves <- matrix(as.numeric(surv_model$predict(data.tune, t_grid)$predictions), nrow = N)

    ## Splits are drawn once and reused across gamma, so the candidates are
    ## compared on exactly the same partitions.
    n_ca <- min(max_cal, floor(N/2))
    splits <- lapply(seq_len(n_repeats), function(r) {
        idx <- sample.int(N)
        list(ca = idx[seq_len(n_ca)],
             te = idx[n_ca + seq_len(min(max_test, N - n_ca))])
    })

    stats_mat <- matrix(NA_real_, nrow = n_repeats, ncol = G)

    for (g in seq_len(G)) {
        ## z_cal and z_test are row-wise in the subject, so they too are
        ## computed once per gamma rather than once per (gamma, repeat).
        inp_full <- conformal_inputs(po_full, data.tune, data.tune, surv_model,
                                     gamma = gamma_grid[g],
                                     curves.cal = curves, curves.test = curves,
                                     t_grid = t_grid)

        for (r in seq_len(n_repeats)) {
            ca <- splits[[r]]$ca
            te <- splits[[r]]$te
            if (!length(te) || !any(po_full$status[ca] == 1L)) next

            pv <- conformal_pvalues(subset_pseudo_outcomes(po_full, ca),
                                    .subset_conformal_inputs(inp_full, ca, te))$p.values

            stats_mat[r, g] <- if (criterion == "bh") {
                sum(p.adjust(pv, method = "BH") <= alpha)
            } else {
                pvp <- pv[is.finite(pv) & pv > 0]
                if (length(pvp)) -2 * sum(log(pvp)) else 0
            }
        }
    }

    mean_stat <- colMeans(stats_mat, na.rm = TRUE)
    sd_stat   <- apply(stats_mat, 2, stats::sd, na.rm = TRUE)
    summary <- data.frame(gamma = gamma_grid, mean_stat = mean_stat,
                          sd_stat = sd_stat, n = n_repeats)

    ## Algorithm S5 selects the plain argmax (ties -> smallest gamma).
    ## one_se = TRUE gives the more conservative 1-SE variant instead.
    gamma_opt <- if (one_se) {
        select_lambda_1se(summary, lambda_col = "gamma", direction = "maximize")$lambda_1se
    } else {
        cand <- which(mean_stat == max(mean_stat, na.rm = TRUE))
        gamma_grid[min(cand)]
    }

    list(gamma_opt = gamma_opt, summary = summary)
}

## Select by the 1-SE rule (glmnet-style); kept for the optional variant above.
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
    se <- if (already_se) x[[sd_col]] else x[[sd_col]] / sqrt(x[[n_col]])
    se[!is.finite(se)] <- NA_real_
    if (direction == "maximize") {
        i_best <- which.max(x[[mean_col]]); thr <- x[[mean_col]][i_best] - se[i_best]
        idx <- which(x[[mean_col]] >= thr)
    } else {
        i_best <- which.min(x[[mean_col]]); thr <- x[[mean_col]][i_best] + se[i_best]
        idx <- which(x[[mean_col]] <= thr)
    }
    j <- idx[which.max(x[[lambda_col]][idx])]
    list(lambda_1se = x[[lambda_col]][j], index_1se = j,
         lambda_best = x[[lambda_col]][i_best], index_best = i_best,
         threshold = thr, mean_best = x[[mean_col]][i_best], se_at_best = se[i_best])
}
