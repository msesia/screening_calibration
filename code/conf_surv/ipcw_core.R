## ============================================================
## Calibration for the conditional risk upon selection.
##
## Flavor-agnostic: everything here operates on a pseudo-outcome object
## (see pseudo_outcomes.R) and is identical for IPCW and AIPCW.  The
## estimator is always the ratio
##
##   rhat(lambda) = mean_i A_lambda(X_i) Psi_i / mean_i A_lambda(X_i),
##
## with influence contribution (13)
##
##   psi_i(lambda) = A_lambda(X_i) {Psi_i - rhat(lambda)} / muhat(lambda),
##
## which is what the delta-method bound and the multiplier band consume.
##
## Depends on: pseudo_outcomes.R.
## ============================================================

## ---------- Utilities ----------
.event_by_t0 <- function(time, status, t0) as.integer(status == 1 & time <= t0)
.at_risk_t0  <- function(time, t0)        as.integer(time >= t0)

.select_by_score <- function(score, lambda, screening_crit = c("low risk","high risk")) {
    screening_crit <- match.arg(screening_crit)
    if (screening_crit == "low risk") as.integer(score >= lambda) else as.integer(score <= lambda)
}

## Binomial one-sided CP bounds
binom_lower_bound <- function(k, n, delta) {
    if (n == 0 || k <= 0) return(0)
    stats::qbeta(delta, k, n - k + 1)
}
binom_upper_bound <- function(k, n, delta) {
    if (n == 0) return(1)
    if (k >= n) return(1)
    stats::qbeta(1 - delta, k + 1, n - k)
}

## Optional bound on Z_i (used for finite-sample UBs)
estimate_M <- function(Z) {
    Zf <- Z[is.finite(Z)]
    if (!length(Zf)) return(1)
    pmax(1, max(Zf))
}

## ---------- Core per-lambda estimator and influence function ----------
## The ONLY input that depends on the flavor is Psi.
.core_pseudo_risk <- function(A, Psi) {
    n <- length(A)
    A <- as.integer(A)
    mu_bar <- mean(A)
    if (mu_bar <= 0) return(list(risk_hat = NA_real_, mu_bar = 0, phi = rep(0, n)))
    theta_bar <- mean(A * Psi)
    risk_hat  <- theta_bar / mu_bar
    ## (13); note that non-selected subjects contribute exactly 0
    phi <- A * (Psi - risk_hat) / mu_bar
    list(risk_hat = risk_hat, mu_bar = mu_bar, phi = phi)
}

## ---------- Pointwise bounds ----------
## Delta-method one-sided UB
pointwise_delta_ub <- function(risk_hat_vec, Phi_mat, delta) {
    n <- nrow(Phi_mat)
    se <- sqrt(pmax(apply(Phi_mat, 2, stats::var), 0) / n)
    z  <- stats::qnorm(1 - delta)
    ub <- pmin(1, risk_hat_vec + z * se)
    ub[is.na(ub)] <- 1
    ub
}

## Nonparametric bootstrap, Algorithm S1: resamples the pairs
## {(X_i, Psi_i)}, so it is flavor-agnostic by construction.  The models
## are NOT refitted inside the loop, matching the conditional perspective.
pointwise_bootstrap_ub <- function(Amat, Psi, delta, B_boot = 1000L) {
    n <- nrow(Amat); K <- ncol(Amat)
    boots <- matrix(NA_real_, nrow = B_boot, ncol = K)
    for (b in seq_len(B_boot)) {
        idx <- sample.int(n, n, replace = TRUE)
        Ab  <- Amat[idx, , drop = FALSE]
        num <- colMeans(Ab * Psi[idx])
        den <- colMeans(Ab)
        boots[b, ] <- ifelse(den > 0, num / den, NA_real_)
    }
    apply(boots, 2, function(x) {
        q <- stats::quantile(x, probs = 1 - delta, na.rm = TRUE, names = FALSE)
        ifelse(is.na(q), 1, min(1, q))
    })
}

## Finite-sample UB: empirical Bernstein for the numerator + exact
## binomial lower bound for the denominator (Supplement S1.1.4).
##
## This always uses the IPCW-ET pseudo-outcome, as specified in
## Supplement S1.1.2, because Psi_dr is not bounded in [0, M].
estimate_risk_fs <- function(A, po, delta = 0.05, M = NULL,
                             method = c("empirical_bernstein","hoeffding")) {
    method <- match.arg(method)
    n <- length(A); A <- as.integer(A)
    if (sum(A) == 0) return(1)

    Z <- po$w_event * (po$time <= po$t0) * A     # in [0, M]
    theta_hat <- mean(Z)
    var_Z <- stats::var(Z)
    mu_low <- binom_lower_bound(sum(A), n, delta/2)
    if (mu_low <= 0) return(1)
    if (is.null(M)) M <- estimate_M(Z)

    phi_upp <- if (method == "empirical_bernstein") {
        theta_hat + sqrt(2 * var_Z * log(4/delta) / n) + (7 * M * log(4/delta)) / (3 * max(1, n - 1))
    } else {
        theta_hat + M * sqrt(log(2/delta) / (2 * n))
    }
    ub <- min(1, phi_upp / mu_low)
    if (!is.finite(ub)) ub <- 1
    ub
}


## ---------- Studentization factors ----------
## Upper confidence bound for the standard deviation of the influence
## contributions at each lambda.  Studentizing by sigma_hat itself is
## unstable: sigma_hat(lambda) can be arbitrarily small, or exactly zero
## (one selected subject, or all selected subjects censored before t0
## under IPCW), which would let a handful of noisy thresholds dictate the
## critical value.  Using an upper bound floors the denominator.
##
## The degrees of freedom are n_eff(lambda), the number of SELECTED
## subjects, not n: phi_i = 0 for every non-selected subject, so
## colMeans(Phi^2) is an average of n_eff nonzero terms over n.  Using
## df = n would give a flat ~4% inflation at every threshold and no
## protection where it is needed -- and an under-estimated sd(lambda)
## narrows the band at exactly the threshold that gets selected.
.sd_upper_bound <- function(Phi, n_eff, gamma = 0.05) {
    s2 <- colMeans(Phi^2)
    sqrt(pmax(nrow(Phi) * s2 / stats::qchisq(gamma, df = pmax(n_eff - 1, 1)), 0))
}

## ---------- Simultaneous (uniform) band ----------
## Studentized Gaussian multiplier band (Supplement S1.1.6).  Unlike the
## unstudentized version, the half-width is lambda-specific:
##
##   UCB(lambda) = rhat(lambda) + q_{1-delta} * sd(lambda) / sqrt(n),
##
## with q_{1-delta} the (1-delta) quantile of max_lambda Z(lambda) and
## Z(lambda) the multiplier process divided by sd(lambda).  Only an upper
## bound on r is needed, so the max is one-sided by default.
##
## `ok` selects the thresholds the band is built on; the rest are handled
## by the finite-sample fallback in compute_calibration_table().
## `ok` selects the thresholds the band is built on; the rest are handled
## by the vacuous value 1 in compute_calibration_table().  `n_eff` is the
## number of SELECTED subjects at each threshold: phi_i is exactly zero
## for non-selected subjects, so sigma_hat(lambda) is estimated from
## n_eff(lambda) nonzero terms, not from n.
.multiplier_halfwidth <- function(Phi, ok = NULL, n_eff = NULL, delta = 0.05,
                                  center = TRUE, B = NULL, one_sided = TRUE,
                                  sd_gamma = 0.05) {
    n <- nrow(Phi); K <- ncol(Phi)
    hw <- rep(NA_real_, K); sd_ub <- rep(NA_real_, K)
    if (is.null(ok))    ok <- rep(TRUE, K)
    if (is.null(n_eff)) n_eff <- rep(n, K)
    if (K == 0 || !any(ok)) return(list(hw = hw, crit = NA_real_, sd = sd_ub))

    Phi_ok    <- Phi[, ok, drop = FALSE]
    sd_ub[ok] <- .sd_upper_bound(Phi_ok, n_eff = n_eff[ok], gamma = sd_gamma)
    K_ok      <- sum(ok)
    if (is.null(B)) B <- max(1000L, ceiling(200 * log(K_ok + 10)))

    xi <- matrix(stats::rnorm(n * B), nrow = n, ncol = B)
    if (center) xi <- sweep(xi, 2, colMeans(xi), "-")
    Z <- t(xi) %*% Phi_ok / sqrt(n)                 # B x K_ok
    Z <- sweep(Z, 2, sd_ub[ok], FUN = "/")          # studentize

    Tstat <- if (one_sided) apply(Z, 1, max) else apply(abs(Z), 1, max)
    q     <- as.numeric(stats::quantile(Tstat, probs = 1 - delta, names = FALSE))

    hw[ok] <- q * sd_ub[ok] / sqrt(n)
    list(hw = hw, crit = q, sd = sd_ub)
}

## ---------- Binomial CP (conservative) upper bound ----------
## Treats censored-before-t0 as failures; does not use pseudo-outcomes.
cp_conservative_ub <- function(A, time, status, t0, delta) {
    A <- as.integer(A)
    m <- sum(A)
    if (m == 0) return(1)
    err_star <- as.integer((status == 1 & time <= t0) | (status == 0 & time < t0))
    k <- sum(A * err_star)
    if (k >= m) 1 else stats::qbeta(1 - delta, k + 1, m - k)
}

## Diagnostic only: number of SELECTED subjects with an OBSERVED event
## before t0.  Reported in the calibration table but NOT used to gate
## anything: at a short horizon a genuinely low-risk subgroup has no
## events, and zero events out of many selected is strong evidence that
## r <= alpha, not an absence of information.
.n_eff_events <- function(Amat, po) {
    ev <- .event_by_t0(po$time, po$status, po$t0)
    as.integer(crossprod(Amat, ev))
}

## ---------- Public API ----------
## compute_calibration_table(scores, po, ...)
##
## `po` is a pseudo-outcome object from build_pseudo_outcomes(); switching
## between IPCW and AIPCW is a change of its `flavor` and nothing else.
## The returned columns are unchanged, so the selectors in
## selection_calibration.R work as before.
##
## Note: with flavor = "dr", rhat(lambda) need not lie in [0,1] (paper
## Section 3.1.2).  It is left unclipped; the bounds are clipped at 1.
compute_calibration_table <- function(
                                      scores, po,
                                      screening_crit = c("low risk","high risk"),
                                      lambda_seq     = NULL,
                                      num_lambda = 100,
                                      delta = 0.05,
                                      include_uniform = TRUE,
                                      include_bootstrap = FALSE, B_boot = 1000L,
                                      include_cp_conservative = FALSE,
                                      M = NULL, min_n_sel = 10,
                                      uniform_one_sided = TRUE,
                                      uniform_sd_gamma = 0.05
                                      ) {
    screening_crit <- match.arg(screening_crit)
    n <- length(scores)
    stopifnot(n == po$n)

    if (is.null(lambda_seq)) {
        lambda_seq <- as.numeric(stats::quantile(scores, probs = seq(0, 1, length.out = num_lambda),
                                                 na.rm = TRUE))
        lambda_seq <- unique(lambda_seq)
    }
    K <- length(lambda_seq)

    ## The pseudo-outcomes at the fixed horizon t0 do not depend on lambda,
    ## so they are computed once (Supplement S1.2).
    Psi <- pseudo_psi_at_time(po, po$t0)

    risk_hat <- numeric(K); mu_bar <- numeric(K)
    Phi  <- matrix(0,  nrow = n, ncol = K)
    Amat <- matrix(0L, nrow = n, ncol = K)

    for (k in seq_len(K)) {
        A <- .select_by_score(scores, lambda_seq[k], screening_crit)
        est <- .core_pseudo_risk(A, Psi)
        risk_hat[k] <- est$risk_hat
        mu_bar[k]   <- est$mu_bar
        Phi[, k]    <- est$phi
        Amat[, k]   <- A
    }
    risk_hat[is.na(risk_hat)] <- 1
    n_sel <- colSums(Amat)
    n_ev  <- .n_eff_events(Amat, po)

    ## Always: pointwise delta UB
    risk_ub_point_delta <- if (K > 0) pointwise_delta_ub(risk_hat, Phi, delta) else numeric(0)

    ## Always: finite-sample UB (used as the small-sample fallback below)
    risk_ub_point_fs <- vapply(seq_len(K), function(k)
        estimate_risk_fs(Amat[, k], po, delta = delta, M = M), numeric(1))

    ## Optional: pointwise bootstrap UB
    risk_ub_point_boot <- NULL
    if (include_bootstrap && K > 0) {
        risk_ub_point_boot <- pointwise_bootstrap_ub(Amat, Psi, delta, B_boot)
    }

    ## Optional: conservative CP UB
    risk_ub_cp_conservative <- NULL
    if (include_cp_conservative && K > 0) {
        risk_ub_cp_conservative <- vapply(seq_len(K), function(k)
            cp_conservative_ub(Amat[, k], po$time, po$status, po$t0, delta), numeric(1))
    }

    ## Optional: simultaneous band, restricted to thresholds where the
    ## asymptotics are justified and the influence contributions are non-degenerate.
    risk_ub_uniform <- NULL
    uniform_crit <- NA_real_
    if (include_uniform && K > 0) {
        sigma_hat <- sqrt(colMeans(Phi^2))
        ok_band   <- (n_sel >= min_n_sel) & is.finite(risk_hat) & (sigma_hat > 0)
        mb <- .multiplier_halfwidth(Phi, ok = ok_band, n_eff = n_sel, delta = delta,
                                    one_sided = uniform_one_sided,
                                    sd_gamma = uniform_sd_gamma)
        risk_ub_uniform <- rep(1, K)
        risk_ub_uniform[ok_band] <- pmin(1, risk_hat[ok_band] + mb$hw[ok_band])
        uniform_crit <- mb$crit
    }
    
    calib <- data.frame(
        lambda        = lambda_seq,
        n_selected    = as.integer(n_sel),
        n_events      = n_ev,
        frac_selected = mu_bar,
        risk_hat      = risk_hat,
        risk_ub_point_delta = risk_ub_point_delta,
        risk_ub_point_fs    = risk_ub_point_fs,
        stringsAsFactors = FALSE
    )
    if (!is.null(risk_ub_point_boot))      calib$risk_ub_point_boot      <- risk_ub_point_boot
    if (!is.null(risk_ub_cp_conservative)) calib$risk_ub_cp_conservative <- risk_ub_cp_conservative
    if (!is.null(risk_ub_uniform))         calib$risk_ub_uniform         <- risk_ub_uniform
    if (!is.null(risk_ub_uniform)) attr(calib, "uniform_crit") <- uniform_crit
    
    ## Fall back to the finite-sample bound where the asymptotics are not
    ## justified (too few selected subjects).
    idx <- which(n_sel < min_n_sel)    
    if (length(idx) > 0) {
        for (ub_col in c("risk_ub_point_delta", "risk_ub_point_boot")) {
            if (ub_col %in% names(calib)) calib[[ub_col]][idx] <- calib[["risk_ub_point_fs"]][idx]
        }
    }

    calib
}
