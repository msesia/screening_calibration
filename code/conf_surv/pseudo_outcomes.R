## ---------------------------------------------------------
## AIPCW preprocessing: one O(nK) pass giving Ghat, Shat and the prefix
## sums (S10), so that each later evaluation of Psi_dr costs O(1).
## ---------------------------------------------------------

## Row-wise cumulative helper.  apply(M, 1, f) returns a *vector* when M
## has a single column, so the result is re-shaped explicitly.
.row_apply <- function(M, f) {
    out <- apply(M, 1, f)
    if (is.matrix(out)) t(out) else matrix(out, nrow = nrow(M))
}


## ---------------------------------------------------------
## Time grid  0 = u_0 < u_1 < ... < u_K = t0.
##
## The observed times are always grid points.  This is what makes the
## indicator I(T_i <= t) exact when a covariate-dependent horizon is
## rounded to the grid, which is what the conformal p-values need.
## `max_grid` bounds the memory of the n x (K+1) tables built for the
## "dr" flavor; at n_cal = 1000 it never binds.
## ---------------------------------------------------------
pseudo_time_grid <- function(time, t0, num_grid = 100, max_grid = 2000) {
    t0 <- as.numeric(t0)
    obs <- as.numeric(time)[time > 0 & time < t0]
    if (length(obs) + 2L > max_grid) {
        obs <- as.numeric(stats::quantile(obs, probs = seq(0, 1, length.out = max_grid - 2L),
                                          na.rm = TRUE))
    }
    u <- c(0, t0, obs)
    if (is.finite(num_grid) && num_grid > 0) {
        u <- c(u, seq(0, t0, length.out = num_grid + 1))
    }
    sort(unique(u))
}


## ---------------------------------------------------------
## AIPCW preprocessing: one O(nK) pass giving Ghat, Shat and the prefix
## sums (S10), so that each later evaluation of Psi_dr costs O(1).
##
##   A_l(x) = sum_{k<=l} dLambda^C_k(x) / G_{k-1}(x)
##   B_l(x) = sum_{k<=l} dLambda^C_k(x) / (G_{k-1}(x) S_{k-1}(x))
##
## Term (c) of (9) is then A_m - S(t|x) B_m, which is exact because
## Q_k(x,t) is affine in S(t|x).
##
## Matrices are n x (K+1); column j holds the value at u_{j-1}.  All
## k / m / ell indices in this file are 0-based, converted with +1L at
## lookup time.
## ---------------------------------------------------------
.aipcw_tables <- function(data, surv_model, cens_model, u, w_event, eps = 1e-12) {
    n <- nrow(data)
    K <- length(u) - 1L

    Shat <- matrix(as.numeric(surv_model$predict(data, u)$predictions), nrow = n)
    Ghat <- matrix(as.numeric(cens_model$predict(data, u)$predictions), nrow = n)
    stopifnot(ncol(Shat) == K + 1L, ncol(Ghat) == K + 1L)

    ## clip to (0,1] and enforce monotone non-increasing rows (fitted step
    ## functions can wobble slightly under interpolation)
    Shat <- .row_apply(pmin(pmax(Shat, eps), 1), cummin)
    Ghat <- .row_apply(pmin(pmax(Ghat, eps), 1), cummin)

    ## the winsorization cap on the IPCW weights is applied to the 1/Ghat
    ## factors in all three terms of (9) (Supplement S1.1.2), by flooring Ghat
    cap <- suppressWarnings(max(w_event[is.finite(w_event)], na.rm = TRUE))
    if (is.finite(cap) && cap > 0) Ghat <- pmax(Ghat, 1 / cap)

    ## discrete censoring hazard: mass 1 - G_k/G_{k-1} at u_k
    dLam <- 1 - Ghat[, -1, drop = FALSE] / Ghat[, -(K + 1L), drop = FALSE]
    dLam <- pmin(pmax(dLam, 0), 1 - eps)

    Gp <- Ghat[, -1, drop = FALSE]             # G_k,     k = 1..K
    Sm <- Shat[, -(K + 1L), drop = FALSE]      # S_{k-1}, k = 1..K
    zero <- matrix(0, nrow = n, ncol = 1L)

    A <- cbind(zero, .row_apply(dLam / Gp, cumsum))
    B <- cbind(zero, .row_apply(dLam / (Gp * Sm), cumsum))

    ## A_l telescopes to 1/G_l - 1; exact unless the dLam clipping bound
    ## above is active.
    stopifnot(max(abs(A - sweep(1 / Ghat, 1, 1 / Ghat[, 1], "-"))) < 1e-6)
    
    list(Shat = Shat, Ghat = Ghat, A = A, B = B)
}


## ---------------------------------------------------------
## Constructor.
##
##   t0            : end of the time grid = largest evaluable horizon
##   flavor        : "et", "ft" or "dr"
##   weights_event : optional pre-computed IPCW-ET weights.  Passing the
##                   ones already used elsewhere makes term (a) of the
##                   AIPCW pseudo-outcome identical to the IPCW one.
##   surv_model    : only needed for "dr"
## ---------------------------------------------------------
build_pseudo_outcomes <- function(data, cens_model = NULL, time, status, t0,
                                  flavor = c("et", "ft", "dr"),
                                  surv_model = NULL,
                                  weights_event = NULL, weights_fixed = NULL,
                                  num_grid = 100, max_grid = 2000,
                                  trim = 0.01, eps = 1e-12) {
    flavor <- match.arg(flavor)
    n <- nrow(data)
    time <- as.numeric(time)
    status <- as.integer(status)
    stopifnot(length(time) == n, length(status) == n)

    u <- pseudo_time_grid(time, t0, num_grid = num_grid, max_grid = max_grid)
    K <- length(u) - 1L

    ## IPCW-ET weights: term (a) of every flavor's numerator, and the
    ## finite-sample bound
    if (is.null(weights_event)) {
        if (is.null(cens_model)) stop("supply either cens_model or weights_event.")
        weights_event <- compute_ipcw_weights(data, cens_model, time, status,
                                              ipcw_method = "et", trim = trim)
    }
    w <- as.numeric(weights_event)
    w[!is.finite(w) | status != 1L] <- 0       # non-events contribute 0

    po <- list(u = u, K = K, n = n, time = time, status = status, t0 = t0,
               w_event = w, flavor = flavor,
               kT = findInterval(time, u, left.open = TRUE),   # min{k : T_i <= u_k}
               varying_horizon = (flavor != "ft"))

    if (flavor == "ft") {
        if (is.null(weights_fixed)) {
            weights_fixed <- compute_ipcw_weights(data, cens_model, time, status,
                                                  ipcw_method = "ft", t0 = t0, trim = trim)
        }
        wf <- as.numeric(weights_fixed)
        wf[!is.finite(wf)] <- 0
        po$w_fixed <- wf
    }

    if (flavor == "dr") {
        if (is.null(surv_model)) stop("flavor = 'dr' requires a fitted surv_model.")
        po <- c(po, .aipcw_tables(data, surv_model, cens_model, u, w, eps = eps))
    }

    po
}


## ---------------------------------------------------------
## Restrict a pseudo-outcome object to a subset of its rows.
##
## Because the survival and censoring models are treated as fixed
## (Assumption (A2)), Psi_i depends only on subject i: w_event, kT and the
## augmentation tables Shat, Ghat, A, B are all row-wise.  The only
## sample-level quantities are the time grid u and the winsorization cap
## folded into Ghat, and BOTH ARE DELIBERATELY INHERITED from the parent
## object, so that every subset shares one grid and one cap.
##
## Consequence: a resampling loop over splits of a fixed data set does not
## need to rebuild anything.  Build once, subset per split.
## ---------------------------------------------------------
subset_pseudo_outcomes <- function(po, idx) {
    idx <- as.integer(idx)
    out <- po
    out$n <- length(idx)
    for (f in c("time", "status", "w_event", "w_fixed", "kT")) {
        if (!is.null(po[[f]])) out[[f]] <- po[[f]][idx]
    }
    for (f in c("Shat", "Ghat", "A", "B")) {
        if (!is.null(po[[f]])) out[[f]] <- po[[f]][idx, , drop = FALSE]
    }
    out
}


## ---------------------------------------------------------
## Psi_i at horizon u[ell_i], vectorised over individuals.
## ---------------------------------------------------------
pseudo_psi <- function(po, ell) {
    n <- po$n
    was_matrix <- is.matrix(ell)
    ell <- if (was_matrix) {
        stopifnot(nrow(ell) == n); matrix(as.integer(ell), n)
    } else {
        e <- as.integer(ell)
        if (length(e) == 1L) e <- rep(e, n)
        stopifnot(length(e) == n)
        matrix(e, n, 1L)
    }
    J <- ncol(ell)

    ## Length-n vectors below recycle down the columns of the n x J matrices,
    ## which is what makes a whole block of horizons cost one vectorised pass.
    pos <- ell >= 1L                          # horizon strictly positive
    elp <- pmax(ell, 1L)                      # guarded index for lookups
    tt <- matrix(po$u[elp + 1L], n, J)
    before <- (po$time <= tt)

    out <- switch(po$flavor,

        ## ---- (7): w_i is already 0 unless E_i = 1
        "et" = po$w_event * before,

        ## ---- fixed-time form; only meaningful at t0
        "ft" = {
            if (!all(ell == po$K)) {
                stop("flavor = 'ft' is only defined at the fixed horizon t0.")
            }
            1 - po$w_fixed * (po$time >= tt)
        },

        ## ---- (9)-(10) with the prefix sums of (S10)
        "dr" = {
            rows <- rep.int(seq_len(n), J)
            k <- pmin(elp, po$kT)             # k_i(t): index of T_i ^ t on the grid
            m <- pmax(k - 1L, 0L)             # left-limit index, for Qhat only
            i_t <- cbind(rows, as.vector(elp) + 1L)
            i_k <- cbind(rows, as.vector(k) + 1L)
            i_m <- cbind(rows, as.vector(m) + 1L)

            S_t <- matrix(po$Shat[i_t], n, J)
            S_m <- matrix(po$Shat[i_m], n, J)
            Q <- pmin(pmax((S_m - S_t) / S_m, 0), 1)          # Qhat(u, t | X_i), eq. (10)

            ## Observed times are grid points, so T_i ^ t = u_k exactly and the
            ## atom at u_k belongs in term (c); indexing A and B at i_m would
            ## drop it, including the atom at t0 for every subject still at risk.
            po$w_event * before +                                        # (a)
                (po$status == 0L) * before * Q / matrix(po$Ghat[i_k], n, J) -  # (b)
                (matrix(po$A[i_k], n, J) - S_t * matrix(po$B[i_k], n, J))      # (c)
        },

        stop("unknown pseudo-outcome flavor: ", po$flavor)
    )

    out[!pos] <- 0
    if (was_matrix) out else as.numeric(out)
}


## ---------------------------------------------------------
## Psi_i at an explicit horizon t (scalar or length n).
## t is rounded up to the next grid point: ell = min{k : t <= u_k}.
## ---------------------------------------------------------
pseudo_psi_at_time <- function(po, t) {
    t <- as.numeric(t)
    if (length(t) == 1L) t <- rep(t, po$n)
    pseudo_psi(po, pmin(findInterval(t, po$u, left.open = TRUE), po$K))
}
