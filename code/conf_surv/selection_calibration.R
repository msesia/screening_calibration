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
select_lambda_greedy <- function(lambda, risk, alpha, size = NULL, mask = NULL,
                                 tie = c("largest_size_then_largest_lambda",
                                         "largest_size_then_smallest_lambda",
                                         "largest_lambda",
                                         "smallest_lambda")) {
    stopifnot(length(lambda) == length(risk))
    K <- length(lambda)
    if (!is.null(size)) stopifnot(length(size) == K)
    if (!is.null(mask)) {
        stopifnot(length(mask) == K)
    } else {
        mask <- rep(TRUE, K)
    }
    tie <- match.arg(tie)
    alpha_vec <- .align_alpha(alpha, K)

    eligible <- which(mask & is.finite(risk) & (risk <= alpha_vec))
    if (length(eligible) == 0) {
        return(list(index = NA_integer_, lambda = NA_real_, eligible = integer(0)))
    }

    if (is.null(size)) {
        ## no size provided → pick by lambda per tie rule
        pick <- switch(tie,
                       "largest_lambda"                    = eligible[which.max(lambda[eligible])],
                       "smallest_lambda"                   = eligible[which.min(lambda[eligible])],
                       "largest_size_then_largest_lambda"  = eligible[which.max(lambda[eligible])],
                       "largest_size_then_smallest_lambda" = eligible[which.min(lambda[eligible])]
                       )
    } else {
        df <- data.frame(idx = eligible, size = size[eligible], lambda = lambda[eligible])
        ord <- switch(tie,
                      "largest_size_then_largest_lambda"  = order(-df$size, -df$lambda),
                      "largest_size_then_smallest_lambda" = order(-df$size,  df$lambda),
                      "largest_lambda"                    = order(-df$lambda),
                      "smallest_lambda"                   = order(df$lambda)
                      )
        pick <- df$idx[ord][1]
    }

    list(index = pick, lambda = lambda[pick], eligible = eligible)
}

## -------------------------------
## 2) Learn-then-Test (fixed seqs)
## -------------------------------

## LTT over user-provided sequences (list of integer vectors of row indices).
## Tests each sequence in order; within a sequence, continue while risk[j] ≤ alpha[j];
## stop that sequence at the first failure. Each index is tested at most once overall.
## 'mask' (if provided) skips indices that are FALSE.
select_lambda_LTT <- function(lambda, risk, alpha, sequences,
                              size = NULL, mask = NULL,
                              tie = c("largest_size_then_largest_lambda",
                                      "largest_lambda",
                                      "largest_size_then_smallest_lambda",
                                      "smallest_lambda")) {
    stopifnot(length(lambda) == length(risk))
    K <- length(lambda)
    if (!is.null(size)) stopifnot(length(size) == K)
    if (!is.null(mask)) {
        stopifnot(length(mask) == K)
    } else {
        mask <- rep(TRUE, K)
    }
    tie <- match.arg(tie)
    alpha_vec <- .align_alpha(alpha, K)

    tested   <- rep(FALSE, K)
    rejected <- rep(FALSE, K)

    for (path in sequences) {
        path <- unique(path[path >= 1 & path <= K])
        for (j in path) {
            if (!mask[j]) next                   ## skip ineligible indices
            if (tested[j]) next                  ## already tested elsewhere
            tested[j] <- TRUE
            if (is.finite(risk[j]) && risk[j] <= alpha_vec[j]) {
                rejected[j] <- TRUE                ## keep walking this path
            } else {
                break                              ## stop this path at first failure
            }
        }
    }
    rej_idx <- which(rejected)
    if (length(rej_idx) == 0) {
        return(list(index = NA_integer_,
                    lambda = NA_real_,
                    rejected_indices = integer(0),
                    tested_indices   = which(tested)))
    }

    ## Choose one λ among rejections
    if (is.null(size)) {
        pick <- switch(tie,
                       "largest_lambda"                    = rej_idx[which.max(lambda[rej_idx])],
                       "smallest_lambda"                   = rej_idx[which.min(lambda[rej_idx])],
                       "largest_size_then_largest_lambda"  = rej_idx[which.max(lambda[rej_idx])],
                       "largest_size_then_smallest_lambda" = rej_idx[which.min(lambda[rej_idx])]
                       )
    } else {
        df <- data.frame(idx = rej_idx, size = size[rej_idx], lambda = lambda[rej_idx])
        ord <- switch(tie,
                      "largest_size_then_largest_lambda"  = order(-df$size, -df$lambda),
                      "largest_lambda"                    = order(-df$lambda),
                      "largest_size_then_smallest_lambda" = order(-df$size,  df$lambda),
                      "smallest_lambda"                   = order(df$lambda)
                      )
        pick <- df$idx[ord][1]
    }

    list(index = pick,
         lambda = lambda[pick],
         rejected_indices = rej_idx,
         tested_indices   = which(tested))
}

## -------------------------------------------------
## Build LTT sequences from an anchor lambda value
## -------------------------------------------------

## Internal: find anchor index (on the lambda *scale*)
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
                                           snap = c("closest","floor","ceiling"),
                                           size = NULL, mask = NULL,
                                           tie = c("largest_size_then_largest_lambda",
                                                   "largest_lambda",
                                                   "largest_size_then_smallest_lambda",
                                                   "smallest_lambda")) {
    snap <- match.arg(snap)
    seqs <- build_sequences_from_anchors(lambda, anchors,
                                          step_up = step_up, step_down = step_down,
                                          snap = snap, include_anchor = TRUE, return = "both")
    select_lambda_LTT(lambda, risk, alpha, sequences = seqs, size = size, mask = mask, tie = tie)
}
