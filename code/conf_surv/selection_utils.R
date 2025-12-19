## ---------------------------------------------
## 1) Per-patient threshold selection
##    - low risk  : keep scores >= cutoff
##    - high risk : keep scores <= cutoff
## ---------------------------------------------
select_patients_threshold <- function(scores, cutoff,
                                      screening_crit = c("low risk","high risk")) {
    screening_crit <- match.arg(screening_crit)
    x <- as.numeric(scores)
    keep <- is.finite(x)

    if (screening_crit == "low risk") {
        idx <- which(keep & x >= cutoff)
    } else { ## high risk
        idx <- which(keep & x <= cutoff)
    }
    if (length(idx) == 0) integer(0) else idx
}

## -----------------------------------------------------------
## 2) Largest subset whose average meets the cutoff
##    - low risk  : largest S with mean(scores[S]) >= cutoff
##    - high risk : largest S with mean(scores[S]) <= cutoff
                                        #
## Greedy is optimal:
##   low risk  → sort d_i =  scores[i] - cutoff  (desc), take longest prefix with cumsum >= 0
##   high risk → sort d_i =  cutoff - scores[i]  (desc), take longest prefix with cumsum >= 0
## -----------------------------------------------------------
select_patients_average <- function(scores, cutoff,
                                    screening_crit = c("low risk","high risk")) {
    screening_crit <- match.arg(screening_crit)
    x <- as.numeric(scores)
    keep <- which(is.finite(x))
    if (!length(keep)) return(integer(0))

    d <- if (screening_crit == "low risk") (x[keep] - cutoff) else (cutoff - x[keep])
    ord <- order(d, decreasing = TRUE)
    cs  <- cumsum(d[ord])
    k   <- max(which(cs >= 0), na.rm = TRUE)

    if (!is.finite(k) || k <= 0) integer(0) else keep[ord[seq_len(k)]]
}
