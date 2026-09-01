## ============================================================
## Semi-synthetic generator with tilts on T | X and C | X.
##
## A survival curve is S(t|x) = exp(-H(t|x)), so raising it to a power
## multiplies the hazard. That is the whole mechanism:
##
##     S_tilt(t|x) = S(t|x) ^ exp(f_S(x))
##     G_tilt(t|x) = G(t|x) ^ exp(f_C(x))
##
## f is any R expression in the covariates and is the LOG hazard
## multiplier. Nothing is centred, normalised or rescaled: exp(f(x)) is
## used exactly as written.
##
##     tilt_surv = list(f = ~  1.5 * X12)          # 1.5 log-HR per unit X12
##     tilt_cens = list(f = ~  log(12) - 1.5 * X12) # 12x baseline, opposite slope
##     tilt_cens = list(f = ~  2.5 + 0.8*(X2 == "white"))
##
## A constant term shifts everyone's hazard: +log(2) doubles it. A slope
## makes the hazard depend on x. Both change the marginal rates, which is
## intended -- you are specifying the hazard, not a normalised version.
##
## f = NULL means no tilt and reproduces the untilted generator exactly.
##
## Sampling is exact:
##     T = inf{t : S_tilt(t|x) <= U} = inf{t : S(t|x) <= U^(1/r)},
## so draw U ~ Unif(0,1), raise to the power 1/r, invert the original
## fitted curve.
##
## $surv_model and $cens_model return the TILTED curves, so the oracle
## benchmark in experiment_1.R stays correct.
## ============================================================


## Curve raised to the per-patient power exp(f(x)).
TiltedCurve <- R6::R6Class("TiltedCurve",
  public = list(
    base = NULL, tilt = NULL,
    initialize = function(base_model, tilt) { self$base <- base_model; self$tilt <- tilt },
    predict = function(new_data, time.points = NULL) {
      pr <- self$base$predict(new_data, time.points = time.points)
      list(predictions = pr$predictions ^ self$tilt$rate(new_data),
           time.points  = pr$time.points)
    }
  )
)


## r(x) = exp(f(x)).
TiltSpec <- R6::R6Class("TiltSpec",
  public = list(
    expr = NULL, env = NULL, label = "none",

    initialize = function(f = NULL) {
      if (inherits(f, "formula")) {
        self$expr <- f[[length(f)]]          ## rhs of  ~ expr
        self$env  <- environment(f)
      } else if (!is.null(f)) {
        self$expr <- f                       ## quote(...) or a call
        self$env  <- parent.frame()
      }
      self$label <- if (is.null(self$expr)) "none" else paste(deparse(self$expr), collapse = " ")
      invisible(self)
    },

    ## f(x), with the covariates in scope
    lp = function(d) {
      if (is.null(self$expr)) return(rep(0, nrow(d)))
      v <- tryCatch(eval(self$expr, as.data.frame(d), self$env),
                    error = function(e)
                        stop(sprintf("tilt: could not evaluate  %s  -- %s",
                                     self$label, conditionMessage(e)), call. = FALSE))
      v <- as.numeric(v)
      if (length(v) == 1L) v <- rep(v, nrow(d))
      if (length(v) != nrow(d))
          stop(sprintf("tilt:  %s  returned %d values for %d rows",
                       self$label, length(v), nrow(d)), call. = FALSE)
      if (any(!is.finite(v)))
          stop(sprintf("tilt:  %s  produced non-finite values", self$label), call. = FALSE)
      v
    },

    rate = function(new_data) exp(self$lp(new_data)),

    describe = function() sprintf("f = %s", self$label)
  )
)


TiltedSemiSyntheticDataGenerator <- R6::R6Class("TiltedSemiSyntheticDataGenerator",
  public = list(
    data_raw = NULL, surv_model = NULL, cens_model = NULL,
    base_surv = NULL, base_cens = NULL, tilt_surv = NULL, tilt_cens = NULL,
    grid = NULL, curves_surv = NULL, curves_cens = NULL,

    #' @param tilt_surv,tilt_cens  list(f = ~ <expression>)
    initialize = function(data, surv_model_type, cens_model_type,
                          tilt_surv = list(), tilt_cens = list()) {
      stopifnot("time" %in% colnames(data), "status" %in% colnames(data))
      self$data_raw <- data

      self$base_surv <- init_surv_model(surv_model_type)
      self$base_surv$fit(Surv(time, status) ~ ., data = data)

      self$base_cens <- CensoringModel$new(model = init_censoring_model(cens_model_type))
      self$base_cens$fit(data = data)

      self$tilt_surv <- TiltSpec$new(tilt_surv$f)
      self$tilt_cens <- TiltSpec$new(tilt_cens$f)
      self$surv_model <- TiltedCurve$new(self$base_surv, self$tilt_surv)
      self$cens_model <- TiltedCurve$new(self$base_cens, self$tilt_cens)

      ps <- self$base_surv$predict(data);  pc <- self$base_cens$predict(data)
      self$curves_surv <- ps$predictions;  self$curves_cens <- pc$predictions
      self$grid <- list(surv = ps$time.points, cens = pc$time.points)
      invisible(self)
    },

    ## Inverse-CDF sampling from the tilted curve. Rows of `curves` are
    ## non-increasing, so the number of entries >= u is the time index.
    .draw = function(curves, times, rate) {
      u <- stats::runif(nrow(curves)) ^ (1 / rate)
      idx <- rowSums(curves >= u)
      idx[idx < 1] <- 1L
      times[idx]
    },

    sample = function(shuffle = TRUE, return_oracle = FALSE) {
      n <- nrow(self$data_raw)
      T_new <- self$.draw(self$curves_surv, self$grid$surv, self$tilt_surv$rate(self$data_raw))
      C_new <- self$.draw(self$curves_cens, self$grid$cens, self$tilt_cens$rate(self$data_raw))

      df <- self$data_raw %>%
        select(-time, -status) %>%
        mutate(event_time = T_new, censoring_time = C_new,
               time = pmin(event_time, censoring_time),
               status = as.integer(event_time <= censoring_time))
      if (shuffle) df <- df[sample(n), ]

      data_observed <- df %>% select(time, status, everything(), -censoring_time, -event_time)
      data_oracle   <- df %>% select(time, status, event_time, censoring_time, everything())
      if (return_oracle) list(observed = data_observed, oracle = data_oracle) else data_observed
    },

    ## Report only; nothing here is corrected for.
    tilt_summary = function(t0 = 2, label = NULL) {
      if (!is.null(label)) cat(sprintf("  [%s]\n", label))
      S  <- as.numeric(self$surv_model$predict(self$data_raw, t0)$predictions)
      G  <- as.numeric(self$cens_model$predict(self$data_raw, t0)$predictions)
      rS <- self$tilt_surv$rate(self$data_raw)
      rC <- self$tilt_cens$rate(self$data_raw)
      smp <- self$sample(shuffle = FALSE, return_oracle = TRUE)$oracle
      kS <- max(1L, min(findInterval(t0, self$grid$surv), ncol(self$curves_surv)))
      kC <- max(1L, min(findInterval(t0, self$grid$cens), ncol(self$curves_cens)))
      qq <- function(v) sprintf("%.2f / %.2f / %.2f  (min %.3g, max %.3g)",
                                stats::quantile(v, .1), stats::median(v),
                                stats::quantile(v, .9), min(v), max(v))
      cat(sprintf("  survival tilt : %s\n  censoring tilt: %s\n",
                  self$tilt_surv$describe(), self$tilt_cens$describe()))
      cat(sprintf("  hazard multiplier p10/med/p90\n    survival : %s\n    censoring: %s\n",
                  qq(rS), qq(rC)))
      cat(sprintf("  S(t0|x): mean %.3f sd %.4f   |   G(t0|x): mean %.3f sd %.4f\n",
                  mean(S), stats::sd(S), mean(G), stats::sd(G)))
      cat(sprintf("  P(T<=t0): untilted %.3f -> tilted %.3f    P(C<=t0): untilted %.3f -> tilted %.3f\n",
                  1 - mean(self$curves_surv[, kS]), 1 - mean(S),
                  1 - mean(self$curves_cens[, kC]), 1 - mean(G)))
      cat(sprintf("  P(T <= t0) = %.3f    censored before t0 = %.0f%%\n",
                  mean(smp$event_time <= t0), 100 * mean(smp$status == 0 & smp$time < t0)))
      cat(sprintf("  at grid max: %.1f%% of T, %.1f%% of C\n",
                  100 * mean(smp$event_time     >= max(self$grid$surv)),
                  100 * mean(smp$censoring_time >= max(self$grid$cens))))
      invisible(list(S = S, G = G, sample = smp))
    }
  )
)


## ============================================================
## Side-by-side view of two generators at one horizon.
##
## Under setup v0t two distributions are in play and it is easy to lose
## track of which number came from which: the ANALYSIS distribution
## supplies the train / calibration / test data and defines the estimand,
## while the CENSORING-MODEL distribution supplies only the sample Ghat is
## fitted on. The rows that matter are the censoring ones -- Ghat is
## correctly specified exactly when they agree.
## ============================================================
compare_generators <- function(gen_analysis, gen_cens, t0 = 2,
                               labels = c("ANALYSIS (defines the estimand)",
                                          "CENSORING-MODEL FIT")) {
    grab <- function(g) {
        S  <- as.numeric(g$surv_model$predict(g$data_raw, t0)$predictions)
        G  <- as.numeric(g$cens_model$predict(g$data_raw, t0)$predictions)
        rS <- g$tilt_surv$rate(g$data_raw); rC <- g$tilt_cens$rate(g$data_raw)
        smp <- g$sample(shuffle = FALSE, return_oracle = TRUE)$oracle
        list(fS = g$tilt_surv$describe(), fC = g$tilt_cens$describe(),
             rS = rS, rC = rC, S = S, G = G,
             pT = mean(smp$event_time <= t0),
             pC = mean(smp$censoring_time <= t0),
             cens = mean(smp$status == 0 & smp$time < t0),
             pileT = mean(smp$event_time >= max(g$grid$surv)),
             pileC = mean(smp$censoring_time >= max(g$grid$cens)))
    }
    a <- grab(gen_analysis); b <- grab(gen_cens)
    med3 <- function(v) sprintf("%.2f / %.2f / %.2f",
                                stats::quantile(v, .1), stats::median(v), stats::quantile(v, .9))
    row <- function(nm, x, y, star = FALSE)
        cat(sprintf("  %-30s %24s %24s%s\n", nm, x, y, if (star) "   <--" else ""))

    cat(sprintf("\n  %-30s %24s %24s\n", sprintf("at t0 = %g", t0), labels[1], labels[2]))
    cat("  ", strrep("-", 78), "\n", sep = "")
    row("survival tilt f",            a$fS, b$fS)
    row("censoring tilt f",           a$fC, b$fC)
    row("surv multiplier p10/50/90",  med3(a$rS), med3(b$rS))
    row("cens multiplier p10/50/90",  med3(a$rC), med3(b$rC), TRUE)
    row("S(t0|x) mean (sd)",  sprintf("%.3f (%.3f)", mean(a$S), stats::sd(a$S)),
                              sprintf("%.3f (%.3f)", mean(b$S), stats::sd(b$S)))
    row("G(t0|x) mean (sd)",  sprintf("%.3f (%.3f)", mean(a$G), stats::sd(a$G)),
                              sprintf("%.3f (%.3f)", mean(b$G), stats::sd(b$G)), TRUE)
    row("P(T <= t0)",         sprintf("%.3f", a$pT), sprintf("%.3f", b$pT))
    row("P(C <= t0)",         sprintf("%.3f", a$pC), sprintf("%.3f", b$pC), TRUE)
    row("censored before t0", sprintf("%.0f%%", 100*a$cens), sprintf("%.0f%%", 100*b$cens), TRUE)
    row("at grid max: T / C",
        sprintf("%.1f%% / %.1f%%", 100*a$pileT, 100*a$pileC),
        sprintf("%.1f%% / %.1f%%", 100*b$pileT, 100*b$pileC))

    dG <- mean(abs(a$G - b$G))
    cat(sprintf("\n  mean |G_analysis(t0|x) - G_fit(t0|x)| = %.4f", dG))
    cat(sprintf("   (sd of G_analysis = %.4f)\n", stats::sd(a$G)))
    cat(sprintf("  -> Ghat is %s\n\n",
                if (dG < 0.02) "essentially CORRECTLY SPECIFIED (the two censoring laws agree)"
                else sprintf("MISSPECIFIED by %.0f%% of the analysis spread",
                             100*dG/max(stats::sd(a$G), 1e-12))))
    invisible(list(analysis = a, cens_fit = b, dG = dG))
}
