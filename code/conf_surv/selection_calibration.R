this_file <- parent.frame(2)$ofile  ## works when sourced via `source()`
this_dir <- dirname(this_file)
source(file.path(this_dir, "ipcw_core.R"))

## =========================================================
## Selectors for calibrated screening thresholds
##   Inputs are plain vectors from your calibration table:
##     - lambda : grid of thresholds (length K)
##     - risk   : risk estimates or *upper bounds* (length K)
##     - alpha  : target risk level (scalar or length K)
##     - size   : OPTIONAL size metric per λ (e.g., n_selected or frac_selected)
##     - mask   : OPTIONAL logical length-K to pre-filter eligible λ (e.g., size constraints)
## =========================================================

## Helper: vectorize alpha to length K
.align_alpha <- function(alpha, K) {
    if (length(alpha) == 1L) rep(alpha, K) else {
                                               stopifnot(length(alpha) == K); alpha
                                           }
}

## (Optional) build straight sequences from starting points for LTT
.build_sequences_from_starts <- function(starts, K, step = 1L) {
    starts <- unique(starts[starts >= 1 & starts <= K])
    lapply(starts, function(s) seq.int(from = s, to = K, by = step))
}

## -------------------------
## 1) Greedy selector
## -------------------------
## Pick the λ with the LARGEST 'size' among indices with risk ≤ alpha (and mask==TRUE).
## If 'size' is NULL, tie-break by lambda according to 'tie'.
## Algorithm S2, lines 7-8: among the thresholds whose upper bound meets the
## target, take the SMALLEST. Yield is non-increasing in lambda, so the
## smallest feasible threshold is the one that selects the most patients.
select_lambda_greedy <- function(lambda, risk, alpha) {
    stopifnot(length(lambda) == length(risk))
    feasible <- which(is.finite(risk) & risk <= .align_alpha(alpha, length(lambda)))
    if (length(feasible) == 0) {
        return(list(index = NA_integer_, lambda = NA_real_, feasible = integer(0)))
    }
    pick <- feasible[which.min(lambda[feasible])]
    list(index = pick, lambda = lambda[pick], feasible = feasible)
}

## Algorithm S4, lines 4-15: walk each path in order, stop it at its first
## failure, and keep the most liberal (smallest) threshold accepted on any
## path. Within a path the last accepted threshold is already the smallest
## accepted on it, so this is exactly min{lambda_hat^(1), lambda_hat^(2)}.
select_lambda_LTT <- function(lambda, risk, alpha, sequences) {
    stopifnot(length(lambda) == length(risk))
    K <- length(lambda)
    alpha_vec <- .align_alpha(alpha, K)

    pick <- NA_integer_
    for (path in sequences) {
        for (j in unique(path[path >= 1 & path <= K])) {
            if (!(is.finite(risk[j]) && risk[j] <= alpha_vec[j])) break   ## stop this path
            if (is.na(pick) || lambda[j] < lambda[pick]) pick <- j        ## keep the most liberal
        }
    }
    list(index  = pick,
         lambda = if (is.na(pick)) NA_real_ else lambda[pick])
}

.closest_anchor_index <- function(lambda, anchor,
                                  snap = c("closest","floor","ceiling")) {
    snap <- match.arg(snap)
    o <- order(lambda)                 ### work on sorted λ for robust floor/ceiling
    lam <- lambda[o]
    K <- length(lam)

    j <- switch(snap,
                "closest" = {
                    j0 <- which.min(abs(lam - anchor))
                    if (length(j0) == 0) 1 else j0[1]
                },
                "floor"   = {
                    idx <- which(lam <= anchor)
                    if (length(idx) == 0) 1 else max(idx)
                },
                "ceiling" = {
                    idx <- which(lam >= anchor)
                    if (length(idx) == 0) K else min(idx)
                }
                )
    list(anchor_sorted = j, order = o)
}

## Main builder: returns indices in the ORIGINAL lambda indexing
build_sequences_from_anchors <- function(lambda, anchors,
                                         step_up = 1L, step_down = 1L,
                                         snap = c("closest","floor","ceiling"),
                                         include_anchor = TRUE,
                                         return = c("both","up","down")) {
    snap <- match.arg(snap)
    return <- match.arg(return)

    sequences <- lapply(anchors, function(anchor) {
        aux <- .closest_anchor_index(lambda, anchor, snap = snap)
        j   <- aux$anchor_sorted
        o   <- aux$order
        K   <- length(lambda)
        down_sorted  <- seq(j, 1L, by = -step_down)
        if (!include_anchor) {
            if (length(down_sorted) > 0) down_sorted <- down_sorted[-1]
        }
        seq_down <- o[down_sorted]  ## map back to original indices
        return(seq_down)
    })
    return(sequences)
}

## -------------------------------------------------
## One-call LTT wrapper from an anchor
## (uses your existing select_lambda_LTT)
## -------------------------------------------------
select_lambda_LTT_from_anchors <- function(lambda, risk, alpha,
                                           anchors,
                                           step_up = 1L, step_down = 1L,
                                           snap = c("closest","floor","ceiling")) {
    snap <- match.arg(snap)
    seqs <- build_sequences_from_anchors(lambda, anchors,
                                          step_up = step_up, step_down = step_down,
                                          snap = snap, include_anchor = TRUE, return = "both")
    select_lambda_LTT(lambda, risk, alpha, sequences = seqs)
}
