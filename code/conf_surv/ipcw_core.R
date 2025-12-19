## ============================================================
## Calibration for conditional risk upon selection (model-agnostic)
##   - score = survival probability at t0 (higher = safer)
##   - HT denominators (unweighted E[A])
##   - IPCW(1) = "et" (PPV-style), IPCW(2) = "ft" (NPV-style)
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
    pmax(1,max(Zf))
}

## ---------- Core per-λ estimator (HT) and influence function ----------
.core_ipcw_risk <- function(A, time, status, t0,
                            ipcw_method = c("ft","et"),
                            weights_event = NULL, weights_fixed = NULL) {
    ipcw_method <- match.arg(ipcw_method)
    n <- length(A)
    A <- as.integer(A)
    if (sum(A) == 0) return(list(risk_hat = NA_real_, mu_bar = 0, phi = rep(0, n)))

    if (ipcw_method == "et") {
        stopifnot(!is.null(weights_event))
        E <- .event_by_t0(time, status, t0)
        w <- weights_event
        w[is.na(w) | E == 0] <- 0                  # non-events contribute 0
        X <- A * w                                  # = A * w * E
        Y <- A
        theta_bar <- mean(X);  mu_bar <- mean(Y)
        risk_hat  <- if (mu_bar > 0) theta_bar / mu_bar else NA_real_
        Xc <- X - theta_bar; Yc <- Y - mu_bar
        phi <- (Xc / mu_bar) - (theta_bar / (mu_bar^2)) * Yc          # IF for θ/μ

    } else { # ft (NPV-style)
        stopifnot(!is.null(weights_fixed))
        R <- .at_risk_t0(time, t0)
        X <- A * weights_fixed * R
        Y <- A
        theta_bar <- mean(X);  mu_bar <- mean(Y)
        S_hat <- if (mu_bar > 0) theta_bar / mu_bar else NA_real_
        risk_hat <- if (!is.na(S_hat)) 1 - S_hat else NA_real_
        Xc <- X - theta_bar; Yc <- Y - mu_bar
        phi <- - (Xc / mu_bar - (theta_bar / (mu_bar^2)) * Yc)        # IF for 1 - θ/μ
    }

    risk_hat[is.na(risk_hat)] <- 1
    list(risk_hat = risk_hat, mu_bar = mu_bar, phi = phi)
}

## ---------- Pointwise bounds ----------
## Always compute: delta-method one-sided UB
pointwise_delta_ub <- function(risk_hat_vec, Phi_mat, delta) {
    n <- nrow(Phi_mat)
    se <- sqrt(pmax(apply(Phi_mat, 2, stats::var), 0) / n)
    z  <- stats::qnorm(1 - delta)
    ub <- pmin(1, risk_hat_vec + z * se)
    ub[is.na(ub)] <- 1
    return(ub)
}

## Optional: subject-level nonparametric bootstrap (slow; conditions on provided weights)
pointwise_bootstrap_ub <- function(scores, time, status, t0, screening_crit,
                                   ipcw_method, weights_event, weights_fixed,
                                   lambda_seq, delta, B_boot = 1000L) {
    n <- length(scores); K <- length(lambda_seq)
    boots <- matrix(NA_real_, nrow = B_boot, ncol = K)
    for (b in seq_len(B_boot)) {
        idx <- sample.int(n, n, replace = TRUE)
        scores_b <- scores[idx]; time_b <- time[idx]; status_b <- status[idx]
        A_b <- sapply(lambda_seq, function(l) .select_by_score(scores_b, l, screening_crit))
        if (ipcw_method == "et") {
            if (is.null(weights_event)) stop("weights_event required for event_time method.")
            w_evt_b <- weights_event[idx]
            for (k in seq_len(K)) {
                estb <- .core_ipcw_risk(A_b[,k], time_b, status_b, t0, "et", weights_event = w_evt_b)
                boots[b, k] <- estb$risk_hat
            }
        } else {
            if (is.null(weights_fixed)) stop("weights_fixed required for fixed_time method.")
            w_fix_b <- weights_fixed[idx]
            for (k in seq_len(K)) {
                estb <- .core_ipcw_risk(A_b[,k], time_b, status_b, t0, "ft", weights_fixed = w_fix_b)
                boots[b, k] <- estb$risk_hat
            }
        }
    }
    apply(boots, 2, function(x) {
        q <- stats::quantile(x, probs = 1 - delta, na.rm = TRUE, names = FALSE)
        ifelse(is.na(q), 1, min(1, q))
    })
}

## Optional: finite-sample (empirical-Bernstein or Hoeffding) UB
estimate_risk_ipcw_fs <- function(selections, time, status, t0,
                                  ipcw_method = c("et","ft"),
                                  weights_event = NULL, weights_fixed = NULL,
                                  delta = 0.05, M = NULL,
                                  method = c("empirical_bernstein","hoeffding")) {
    ipcw_method <- match.arg(ipcw_method)
    method <- match.arg(method)
    n <- length(selections)
    A <- as.integer(selections)
    if (sum(A) == 0) return(1)

    if (ipcw_method == "et") {
        if (is.null(weights_event)) stop("weights_event required for event_time.")
        E <- .event_by_t0(time, status, t0)
        w <- weights_event
        w[is.na(w) | E == 0] <- 0
        Z <- w * A                           # = w * A * E ∈ [0, M]
    } else {
        if (is.null(weights_fixed)) stop("weights_fixed required for fixed_time.")
        R <- .at_risk_t0(time, t0)
        w <- weights_fixed
        w[!is.finite(w) | R == 0] <- 0
        Z <- w * A * R                       # ∈ [0, M]
    }

    theta_hat <- mean(Z)
    var_Z <- stats::var(Z)
    m <- sum(A)
    mu_low <- binom_lower_bound(m, n, delta/4)
    if (mu_low <= 0) return(1)

    if (is.null(M)) M <- estimate_M(Z)

    if (method == "empirical_bernstein") {
        phi_upp <- theta_hat + sqrt(2 * var_Z * log(4/delta) / n) + (7 * M * log(4/delta)) / (3 * max(1, n - 1))
    } else { # Hoeffding
        phi_upp <- theta_hat + M * sqrt(log(4/delta) / (2 * n))
    }
    ub <- min(1, phi_upp / mu_low)
    if (!is.finite(ub)) ub <- 1
    ub[is.na(ub)] <- 1
    ub
}

## ---------- Simultaneous (uniform) band ----------
.multiplier_halfwidth <- function(Phi, delta = 0.05, center = TRUE, B = NULL) {
    n <- nrow(Phi); K <- ncol(Phi)
    if (K == 0) return(list(hw = NA_real_, crit = NA_real_))
    if (is.null(B)) B <- max(1000L, ceiling(200 * log(K + 10)))
    xi <- matrix(rnorm(n * B), nrow = n, ncol = B)
    if (center) xi <- sweep(xi, 2, colMeans(xi), "-")
    Z <- t(xi) %*% Phi / sqrt(n)              # B x K
    Tstat <- apply(abs(Z), 1, max)            # sup over λ
    q <- as.numeric(stats::quantile(Tstat, probs = 1 - delta, names = FALSE))
    list(hw = q / sqrt(n), crit = q)
}

## ---------- Binomial CP (conservative) upper bound ----------
## Treat censored-before-t0 as failures; fully observed 0/1 within A.
cp_conservative_ub <- function(A, time, status, t0, delta) {
  A <- as.integer(A)
  m <- sum(A)
  if (m == 0) return(1)
  # error* within A: 1 if (event by t0) OR (censored before t0)
  err_star <- as.integer( (status == 1 & time <= t0) | (status == 0 & time < t0) )
  k <- sum(A * err_star)
  # one-sided CP upper bound for p = P(error | A=1)
  if (k >= m) 1 else stats::qbeta(1 - delta, k + 1, m - k)
}

## ---------- Public API ----------
## All enabled bounds are included as columns in the returned data.frame
compute_calibration_table <- function(
                                      scores, time, status, t0,
                                      screening_crit = c("low risk","high risk"),
                                      ipcw_method    = c("ft","et"),
                                      weights_event  = NULL,
                                      weights_fixed  = NULL,
                                      lambda_seq     = NULL,
                                      num_lambda = 100,
                                      delta = 0.05,
                                      include_uniform = TRUE,
                                      include_bootstrap = FALSE, B_boot = 1000L,
                                      include_cp_conservative = FALSE,
                                      M = NULL, min_n_asymptotics = 20
                                      ) {
    screening_crit <- match.arg(screening_crit)
    ipcw_method    <- match.arg(ipcw_method)
    fs_method      <- "empirical_bernstein"

    n <- length(scores)
    stopifnot(length(time)==n, length(status)==n)

    ## λ grid
    if (is.null(lambda_seq)) {
        lambda_seq <- as.numeric(stats::quantile(scores, probs = seq(0, 1, length.out = num_lambda), na.rm = TRUE))
        lambda_seq <- unique(lambda_seq)
    }

    K <- length(lambda_seq)
    risk_hat <- numeric(K)
    mu_bar   <- numeric(K)
    Phi      <- matrix(0, nrow = n, ncol = K)

    ## Per-λ estimates & IF
    for (k in seq_len(K)) {
        A <- .select_by_score(scores, lambda_seq[k], screening_crit)
        est <- .core_ipcw_risk(A, time, status, t0, ipcw_method,
                               weights_event = weights_event, weights_fixed = weights_fixed)
        risk_hat[k] <- est$risk_hat
        mu_bar[k]   <- est$mu_bar
        Phi[,k]     <- est$phi
    }
    risk_hat[is.na(risk_hat)] <- 1
    
    ## Always: pointwise delta UB
    risk_ub_point_delta <- if (K > 0) pointwise_delta_ub(risk_hat, Phi, delta) else numeric(0)

    ## Optional: pointwise bootstrap UB
    risk_ub_point_boot <- NULL
    if (include_bootstrap && K > 0) {
        risk_ub_point_boot <- pointwise_bootstrap_ub(
            scores, time, status, t0, screening_crit,
            ipcw_method, weights_event, weights_fixed,
            lambda_seq, delta, B_boot
        )
    }

    ## Optional: finite-sample UB
    risk_ub_point_fs <- NULL
    if (fs_method != "none" && K > 0) {
        risk_ub_point_fs <- vapply(seq_len(K), function(k) {
            A_k <- .select_by_score(scores, lambda_seq[k], screening_crit)
            estimate_risk_ipcw_fs(
                selections = A_k, time = time, status = status, t0 = t0,
                ipcw_method = if (ipcw_method=="et") "et" else "ft",
                weights_event = weights_event, weights_fixed = weights_fixed,
                delta = delta, M = M, method = if (fs_method=="empirical_bernstein") "empirical_bernstein" else "hoeffding"
            )
        }, numeric(1))
    }

    ## Optional: conservative CP UB (no surrogate needed)
    risk_ub_cp_conservative <- NULL
    if (include_cp_conservative && K > 0) {
        risk_ub_cp_conservative <- vapply(seq_len(K), function(k) {
            A_k <- .select_by_score(scores, lambda_seq[k], screening_crit)
            cp_conservative_ub(A_k, time, status, t0, delta)
        }, numeric(1))
    }

    ## Simultaneous (uniform) UB
    risk_ub_uniform <- NULL
    if (include_uniform && K > 0) {
        mb <- .multiplier_halfwidth(Phi, delta = delta)
        risk_ub_uniform <- pmin(1, risk_hat + mb$hw)
        risk_ub_uniform[is.na(risk_ub_uniform)] <- 1
    }

    ## Assemble the single calibration table (only columns for methods you used)
    calib <- data.frame(
        lambda        = lambda_seq,
        n_selected    = as.integer(round(mu_bar * n)),
        frac_selected = mu_bar,
        risk_hat      = risk_hat,
        risk_ub_point_delta = risk_ub_point_delta,
        stringsAsFactors = FALSE
    )
    if (!is.null(risk_ub_point_boot))       calib$risk_ub_point_boot <- risk_ub_point_boot
    if (!is.null(risk_ub_point_fs))         calib$risk_ub_point_fs   <- risk_ub_point_fs
    if (!is.null(risk_ub_cp_conservative))  calib$risk_ub_cp_conservative  <- risk_ub_cp_conservative
    if (!is.null(risk_ub_uniform))          calib$risk_ub_uniform    <- risk_ub_uniform

    ## Pick which UB column you actually use downstream
    for(ub_col in c("risk_ub_point_delta", "risk_ub_point_boot", "risk_ub_uniform")) {
        if ("risk_ub_point_fs" %in% names(calib)) {
            fallback <- calib[["risk_ub_point_fs"]]
        } else {
            fallback <- rep(1, nrow(calib))
        }
        idx <- which(calib$n_selected < min_n_asymptotics)
        if(length(idx)>0) {
            if(ub_col %in% names(calib)) {
                calib[[ub_col]][idx] <- fallback[idx]
            }
        }
    }
    
    return(calib)
}
