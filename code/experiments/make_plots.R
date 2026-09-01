options(width = 300)

library(tidyverse)
library(kableExtra)
library(aod)           ## betabin(), for the exceedance table

## ============================================================
## WHAT THIS SCRIPT MAKES
##
## Everything is driven by the FIGURES and TABLES manifests below: one entry
## per output. Comment an entry out to skip it, add one to get another. No
## nested loops over every possible configuration.
##
## The pseudo-outcome is chosen PER ENTRY with `flavour`, so one run can
## produce the IPCW and the AIPCW version of the same figure. FLAVOUR below
## is only the default for entries that do not name one.
## ============================================================

MAKE_FIGURES <- TRUE
MAKE_TABLES  <- TRUE

FLAVOUR <- "dr"                       ## DEFAULT flavour; override per entry
DELTA   <- 0.10                       ## confidence level for the HP methods
N_TRAIN <- 5000
PROB    <- 0.9                        ## target survival level


MIN_NUM_SEL_SHOW <- 10

## Exceedance is withheld for a method/horizon whose share of non-empty
## replicates falls below this.
MIN_NONEMPTY_EXC <- 0.1

## ---- FIGURES -----------------------------------------------------------
## One entry per output actually used in the paper.
##
##   file     exact output basename -- must match the \includegraphics path
##            in the manuscript, so do not rename these casually
##   kind     "boxplot" | "delta" | "lines" | "paired" | "paired-lines"
##            | "flavour-lines"
##   flavour  "et" (IPCW), "ft" (fixed-time IPCW) or "dr" (AIPCW). Optional;
##            defaults to FLAVOUR. Ignored by the "paired", "paired-lines"
##            and "flavour-lines" kinds, which always contrast et with dr.
##            The file name gains FLAVOUR_SUFFIX[[flavour]], so the same
##            entry listed twice with different flavours writes two files.
##   small    TRUE -> Yield + Survival rows only; FALSE -> adds conditional rows
##   xlog     line kinds only: log10 x axis (default TRUE). Set FALSE for a
##            linear axis -- right for a count, wrong for a sample size
##            spanning orders of magnitude.
##   time     screening horizon, for the line kinds only (default 2). The
##            boxplot/delta kinds show every horizon as facet columns instead.
##            The horizon must exist in the results: TIME_LIST in the submit
##            script currently gives v0c and v0t horizon 2 only.
##   w, h     figure size in inches, passed straight to ggsave()
##
## Comment an entry out to skip it.
FIGURES <- list(
  ## ---- main paper ----
  list(id = "Fig 1",  file = "fig1_boxplots_0_gen_grf_surv_grf_ntrain5000_ncal1000",
       kind = "boxplot", setup = "v0", real = 0, surv = "grf", n_cal = 1000, small = FALSE,
       flavour = "dr", w = 7, h = 4.5, y_floor=0.6, note = "semi-synthetic, well-specified GRF"),
  list(id = "Fig 2",  file = "fig2_boxplots_0_gen_grf_surv_grf_ntrain5000_ncal1000_delta",
       kind = "delta",   setup = "v0", real = 0, surv = "grf", n_cal = 1000,
       flavour = "dr", w = 7, h = 3, y_floor=0.6, note = "HP-LTT across delta"),

  ## the same two, with the AIPCW pseudo-outcome. Same `file`; the flavour
  ## suffix keeps the outputs apart.
  ## , list(id = "Fig 1 (dr)", file = "boxplots_0_gen_grf_surv_grf_ntrain5000_ncal1000",
  ##        kind = "boxplot", setup = "v0", real = 0, surv = "grf", n_cal = 1000, small = FALSE,
  ##        flavour = "dr", w = 7, h = 5.2, note = "semi-synthetic GRF, AIPCW")
  ## , list(id = "Fig 2 (dr)", file = "boxplots_0_gen_grf_surv_grf_ntrain5000_ncal1000_delta",
  ##        kind = "delta",   setup = "v0", real = 0, surv = "grf", n_cal = 1000,
  ##        flavour = "dr", w = 7, h = 3.2, note = "HP-LTT across delta, AIPCW")

  list(id = "Fig 3",  file = "fig3_boxplots_0_gen_grf_surv_xgb_ntrain5000_ncal1000_small",
       kind = "boxplot", setup = "v0", real = 0, surv = "xgb", n_cal = 1000, small = TRUE,
       w = 7, h = 3, note = "misspecified XGB survival model"),
  ## IPCW vs AIPCW as the calibration sample grows -- the Figure 5 layout
  ## applied to the Figure 4 experiments. n_cal is the swept variable, so it
  ## is deliberately not held fixed; the log axis suits a range spanning
  ## orders of magnitude.

  list(id = "Fig 4",  file = "fig4_lines_0_gen_grf_ntrain5000_small",
       kind = "lines",   setup = "v0c", real = 0, x = "n_cal", time = 2,
       xlab = "Calibration sample size",
       w = 7, h = 3, note = "yield vs calibration sample size"),

  list(id = "Fig 5", file = "fig5_lines_flavour_0_gen_grf_surv_grf_ntrain5000_mix",
       kind = "flavour-lines", setup = "v0t", real = 0, surv = "grf", n_cal = 1000,
       x = "cens_mix", xlog = FALSE, time = 2, xlab = "Censoring-model misspecification",
       w = 7, h = 3, note = "IPCW vs AIPCW per method, GRF, t0=2"),

  list(id = "Fig 6", file = "fig6_lines_gamma_0_gen_grf_ntrain5000",
       kind = "lines", setup = "v0c", real = 0, x = "n_cal", time = 2, n_cal = 1000,
       xlab = "Calibration sample size",
       methods = c(setNames(as.list(sprintf("conformal-dr-g%g (stable)", c(0, 1, 10, 100, 1000))),
                            sprintf("gamma = %g", c(0, 1, 10, 100, 1000))),
                   list("tuned" = "conformal-dr (stable)")),
       color_value = function(lab) if (lab == "tuned") NA_real_ else
                                                                    as.numeric(sub("^gamma = ", "", lab)),
       clab = "Time shift parameter", errorbars=FALSE,
       w = 7, h = 3, note = "conformal across fixed gamma, plus tuned"),

  
  ## ## ---- supplement, simulated data ----
  list(id = "Fig S1", file = "boxplots_0_gen_grf_surv_cox_ntrain5000_ncal1000",
       kind = "boxplot", setup = "v0", real = 0, surv = "cox", n_cal = 1000, small = FALSE,
       w = 7, h = 4.5, note = "misspecified Cox survival model"),
  list(id = "Fig S2", file = "boxplots_0_gen_grf_surv_grf_ntrain5000_ncal100",
       kind = "boxplot", setup = "v0", real = 0, surv = "grf", n_cal = 100, small = FALSE,
       w = 7, h = 4.5, note = "calibration sample size 100"),
  list(id = "Fig S3", file = "boxplots_0_gen_grf_surv_grf_ntrain5000_ncal100_delta",
       kind = "delta",   setup = "v0", real = 0, surv = "grf", n_cal = 100,
       w = 7, h = 3, note = "delta sweep, calibration sample size 100"),
  list(id = "Fig S4", file = "lines_flavour_0_gen_grf_surv_grf_ntrain5000_ncal",
       kind = "flavour-lines", setup = "v0c", real = 0, surv = "grf", n_cal = 1000,
       x = "n_cal", xlog = TRUE, time = 2, xlab = "Calibration sample size",
       w = 9, h = 3.4, note = "IPCW vs AIPCW per method, vs n_cal"),

  ## ---- supplement, real data ----
  list(id = "Fig S5", file = "boxplots_1_gen_grf_surv_xgb_ntrain5000_ncal1000",
       kind = "boxplot", setup = "v1", real = 1, surv = "xgb", n_cal = 1000, small = FALSE,
       w = 7, h = 8, note = "real data, gradient boosting"),
  list(id = "Fig S6", file = "boxplots_1_gen_grf_surv_cox_ntrain5000_ncal1000",
       kind = "boxplot", setup = "v1", real = 1, surv = "cox", n_cal = 1000, small = FALSE,
       w = 7, h = 8, note = "real data, Cox"),
  list(id = "Fig S7", file = "boxplots_1_gen_grf_surv_grf_ntrain5000_ncal1000",
       kind = "boxplot", setup = "v1", real = 1, surv = "grf", n_cal = 1000, small = FALSE,
       w = 7, h = 8, note = "real data, GRF"),
  list(id = "Fig S8", file = "boxplots_1_gen_grf_surv_xgb_ntrain5000_ncal1000_delta",
       kind = "delta",   setup = "v1", real = 1, surv = "xgb", n_cal = 1000,
       w = 7, h = 4, note = "real data, delta sweep, gradient boosting"),
  list(id = "Fig S9", file = "boxplots_1_gen_grf_surv_cox_ntrain5000_ncal1000_delta",
       kind = "delta",   setup = "v1", real = 1, surv = "cox", n_cal = 1000,
       w = 7, h = 4, note = "real data, delta sweep, Cox"),
  list(id = "Fig S10", file = "boxplots_1_gen_grf_surv_grf_ntrain5000_ncal1000_delta",
       kind = "delta",   setup = "v1", real = 1, surv = "grf", n_cal = 1000,
       w = 7, h = 4, note = "real data, delta sweep, GRF")


  ## ---- double robustness: sweep the CENSORING model misspecification (v0t) ----
  ## The censoring model is fitted on a mixture of the analysis censoring law
  ## and an alternative one, so Ghat goes from correctly specified (0) to
  ## maximally wrong (1) while the analysis data is fixed. One facet per rule,
  ## with IPCW and AIPCW overlaid: Prop. 3.1 says AIPCW should stay flat while
  ## IPCW drifts. Oracle and Model are omitted -- they use no pseudo-outcome.
  ## These kinds show both flavours, so `flavour` does not apply.
  ## list(id = "Fig T2", file = "lines_diff_0_gen_grf_surv_grf_ntrain5000_mix",
  ##      kind = "paired-lines", setup = "v0t", real = 0, surv = "grf", n_cal = 1000,
  ##      x = "cens_mix", xlog = FALSE, time = 2, xlab = "Censoring-model misspecification",
  ##      w = 9, h = 3.4, note = "AIPCW - IPCW per method, t0=2"),

  ## ## ---- IPCW vs AIPCW, paired within replicate (not in the paper yet) ----
  ## list(id = "Diff-sim", file = "boxplots_diff_0_gen_grf_surv_grf_ntrain5000_ncal1000",
  ##      kind = "paired",  setup = "v0", real = 0, surv = "grf", n_cal = 1000,
  ##      w = 7, h = 3, note = "AIPCW - IPCW, semi-synthetic GRF"),
  ## list(id = "Diff-real", file = "boxplots_diff_1_gen_grf_surv_xgb_ntrain5000_ncal1000",
  ##      kind = "paired",  setup = "v1", real = 1, surv = "xgb", n_cal = 1000,
  ##      w = 7, h = 4.5, note = "AIPCW - IPCW, real data")
)
if(!MAKE_FIGURES) FIGURES <- list()

## ---- TABLES ------------------------------------------------------------
## Same `flavour` field as the figures.
##   kind  "summary" (default) | "exceedance" | "exceedance-delta"
TABLES <- list(
    list(id = "Table 2",  file = "data_1_ntrain_5000_ncal_1000_xgb",
         setup = "v1", real = 1, surv = "xgb", n_cal = 1000, note = "real data, gradient boosting"),

    list(id = "Table S1", kind = "exceedance",
         file = "exceedance_0_ntrain_5000_ncal_1000",
         setup = "v0", real = 0, n_cal = 1000,
         note = "mean survival and exceedance"),
    list(id = "Table S2", kind = "exceedance-delta",
         file = "exceedance_delta_0_ntrain_5000_ncal_1000_t23",
         setup = "v0", real = 0, n_cal = 1000, times = c(2, 3),
         note = "HP rules across delta, t0 = 2 and 3"),
    list(id = "Table S3", kind = "exceedance-delta",
         file = "exceedance_delta_0_ntrain_5000_ncal_1000_t69",
         setup = "v0", real = 0, n_cal = 1000, times = c(6, 9),
         note = "HP rules across delta, t0 = 6 and 9"),
    list(id = "Table S4", kind = "exceedance-delta",
         file = "exceedance_data_0_ntrain_5000_ncal_1000_xgb",
         setup = "v0", real = 0, surv = "xgb", n_cal = 1000, deltas = 0.1, note = "matches Fig 3"),
    list(id = "Table S5", kind = "exceedance-delta",
         file = "exceedance_data_0_ntrain_5000_ncal_1000_cox",
         setup = "v0", real = 0, surv = "cox", n_cal = 1000, deltas = 0.1, note = "matches Fig S1"),
    list(id = "Table S6", kind = "exceedance-delta",
         file = "exceedance_data_0_ntrain_5000_ncal_100_grf",
         setup = "v0", real = 0, surv = "grf", n_cal = 100,  deltas = 0.1, note = "matches Fig S2"),
    list(id = "Table S7", file = "data_1_ntrain_5000_ncal_1000_cox",
         setup = "v1", real = 1, surv = "cox", n_cal = 1000, note = "matches Fig S5"),
    list(id = "Table S8", file = "data_1_ntrain_5000_ncal_1000_grf",
         setup = "v1", real = 1, surv = "grf", n_cal = 1000, note = "matches Fig S6"),
    list(id = "Table S9",  file = "data_1_ntrain_5000_ncal_1000_xgb",
         flavour = "et", setup = "v1", real = 1, surv = "xgb", n_cal = 1000, note = "real data, gradient boosting")
)
if(!MAKE_TABLES) TABLES <- list()

## Suffix appended to `file` per flavour. The "et" entry is EMPTY so the IPCW
## outputs keep the exact names the manuscript \includegraphics/\input
## already reference; changing it means editing the LaTeX too.
FLAVOUR_SUFFIX <- c(et = "", ft = "_ipcwft", dr = "_aipcw")

OUTDIR_FIG <- "figures"
OUTDIR_TAB <- "tables"


## ============================================================
## Method lists
## ============================================================

## Method keys are "<flavour>|<delta>|<rule>"; the conformal arm has its own
## key per flavour. Only "dr" has a separate conformal arm in the results;
## "et" and "ft" both read the plain one.
method_list <- function(flavour = FLAVOUR, delta = DELTA) {
    list(
        "Oracle"        = "oracle (group)",
        "Model"         = "model (group)",
        "FDR-Conformal" = if (flavour == "dr") "conformal-dr (stable)" else "conformal (stable)",
        "HP-Pointwise"  = sprintf("%s|%.2f|pt_delta",  flavour, delta),
        "HP-Uniform"    = sprintf("%s|%.2f|uniform",   flavour, delta),
        "HP-LTT"        = sprintf("%s|%.2f|ltt_delta", flavour, delta)
    )
}

## Delta sweep: HP-LTT at several confidence levels.
method_list_delta <- function(flavour = FLAVOUR) {
    ds <- c(0.05, 0.1, 0.2, 0.5)
    c(list("Oracle" = "oracle (group)",
           "FDR-Conformal" = if (flavour == "dr") "conformal-dr (stable)" else "conformal (stable)"),
      setNames(as.list(sprintf("%s|%.2f|ltt_delta", flavour, ds)), sprintf("HP-LTT-%g", ds)))
}

## Delta sweep for the exceedance TABLE, which can carry more rows than the
## Figure 2 boxplot comfortably shows. Extra rules cost nothing here, so
## every HP rule appears at every delta.
method_list_delta_table <- function(flavour = FLAVOUR,
                                    rules  = c("pt_delta", "uniform", "ltt_delta"),
                                    labels = c("HP-Pointwise", "HP-Uniform", "HP-LTT"),
                                    ds     = c(0.05, 0.1, 0.2, 0.5)) {
    out <- list("Oracle" = "oracle (group)", "Model"  = "model (group)",
                "FDR-Conformal" = if (flavour == "dr") "conformal-dr (stable)"
                                  else "conformal (stable)")
    for (i in seq_along(rules))
        out <- c(out, setNames(as.list(sprintf("%s|%.2f|%s", flavour, ds, rules[i])),
                               sprintf("%s-%g", labels[i], ds)))
    out
}

## IPCW / AIPCW contrast, and survival models (used by the v0t figures).
COLORS_FLAVOUR <- c("IPCW" = "grey30", "AIPCW" = "red2")
COLORS_MODEL   <- c("Cox" = "darkgreen", "GRF" = "blue2", "XGB" = "orange2")

COLORS <- c("Oracle" = "orange1", "Model" = "grey50", "FDR-Conformal" = "red2",
            "FDR-Conformal (tuned)" = "red2", "FDR-Conformal (untuned)" = "grey30",
            "HP-Pointwise" = "darkgreen", "HP-Uniform" = "green", "HP-LTT" = "blue2",
            "HP-LTT-0.05" = "blue2", "HP-LTT-0.1" = "dodgerblue3",
            "HP-LTT-0.2" = "deepskyblue2", "HP-LTT-0.5" = "skyblue1")


## ============================================================
## Loading
## ============================================================

## Returns NULL with a warning rather than failing, so one missing setup
## does not abort the run.
load_results <- function(setup) {
    idir <- sprintf("results_hpc/%s", setup)
    if (!dir.exists(idir)) {
        warning(sprintf("results directory not found: %s -- its figures are skipped.", idir),
                call. = FALSE, immediate. = TRUE)
        return(NULL)
    }
    files <- list.files(idir, full.names = TRUE)
    sz <- file.info(files)$size
    files <- files[!is.na(sz) & sz > 0]
    if (!length(files)) {
        warning(sprintf("no usable files in %s -- skipped.", idir), call. = FALSE, immediate. = TRUE)
        return(NULL)
    }
    parts <- lapply(files, function(f)
        tryCatch(suppressWarnings(read_delim(f, delim = ",", col_types = cols(), guess_max = 2)),
                 error = function(e) NULL))
    parts <- parts[!vapply(parts, function(d) is.null(d) || nrow(d) == 0, logical(1))]
    if (!length(parts)) return(NULL)
    bind_rows(parts)
}

## Each setup is read once, however many outputs use it.
.cache <- new.env(parent = emptyenv())
get_results <- function(setup) {
    if (!exists(setup, envir = .cache)) assign(setup, load_results(setup), envir = .cache)
    get(setup, envir = .cache)
}

## Which Method keys exist. Run this first if a figure comes out empty.
report_methods <- function(results, label = "") {
    if (is.null(results) || !nrow(results)) {
        cat(sprintf("[%s] no results loaded\n", label)); return(invisible(NULL))
    }
    have <- sort(unique(results$Method))
    for (fl in c("et", "ft", "dr")) {
        miss <- setdiff(unlist(method_list(fl), use.names = FALSE), have)
        cat(sprintf("[%s] flavour '%s': %s\n", label, fl,
                    if (length(miss)) paste("MISSING", paste(miss, collapse = ", ")) else "complete"))
    }
    invisible(have)
}

## Optional spec field with a default.
`%or%` <- function(a, b) if (is.null(a)) b else a

## The flavour an entry asks for, falling back to the global default.
spec_flavour <- function(spec) spec$flavour %or% FLAVOUR

## The methods an entry compares. `methods` in the spec overrides the
## standard list, which is how a figure can contrast two variants of one
## method rather than the usual six. Names are display labels and need
## matching entries in COLORS.
spec_methods <- function(spec) spec$methods %or% method_list(spec_flavour(spec))

## x axis for the line kinds. A log axis suits a sample size spanning
## orders of magnitude; a fraction or a count reads better linearly, with a
## break at each value actually run.
x_scale <- function(spec, xvals) {
    if (isTRUE(spec$xlog %or% TRUE)) scale_x_continuous(trans = "log10")
    else scale_x_continuous(breaks = sort(unique(xvals)))
}

## The line kinds pick a single horizon. Asking for one the results do not
## contain otherwise yields an empty figure with no explanation.
check_horizon <- function(results, spec) {
    want <- spec$time %or% 2
    have <- sort(unique(results$Time))
    if (want %in% have) return(TRUE)
    warning(sprintf("%s: horizon t0 = %s not in these results (available: %s) -- skipped.",
                    spec$id, want, paste(have, collapse = ", ")),
            call. = FALSE, immediate. = TRUE)
    FALSE
}

## An empty data frame makes facet_grid() error; skip that output instead.
empty_skip <- function(df, what) {
    if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) {
        warning(sprintf("no matching rows for %s -- skipped.", what), call. = FALSE, immediate. = TRUE)
        TRUE
    } else FALSE
}


## ============================================================
## Shared reshaping
## ============================================================

## Which facet rows appear, by figure type.
row_vars <- function(real, small) {
    if (real == 1) {
        ## Because of censoring the event rate is only bounded, so the two
        ## bounds swap when moving from the survival to the event scale:
        ## an upper bound on survival is a LOWER bound on the event rate.
        if (small) c("Screened" = "Yield",
                     "FDP_upper" = "FDP (UB)", "FDP_lower" = "FDP (LB)")
        else c("Screened" = "Yield",
               "FDP_upper" = "FDP (UB)", "FDP_lower" = "FDP (LB)",
               "Event rate (UB) | S>0" = "Event rate (UB)",
               "Event rate (LB) | S>0" = "Event rate (LB)",
               "P(S>0)" = "P(|S| > 0)")
    } else {
        if (small) c("Screened" = "Yield", "FDP" = "FDP")
        else c("Screened" = "Yield", "FDP" = "FDP",
               "Event rate | S>0" = "Event rate",
               "P(S>0)" = "P(|S| > 0)")
    }
}

prepare_long <- function(results, vars, spec, methods) {
    ## Variables whose dashed reference line sits at the target survival
    ## PROB, and those whose line sits at the target event rate 1 - PROB.
    surv_scale  <- c("Survival", "Survival | S>0", "Survival_lower", "Survival_upper",
                     "Surv. (LB) | S>0", "Surv. (UB) | S>0")
    event_scale <- c("FDP", "Event rate | S>0", "FDP_lower", "FDP_upper",
                     "Event rate (LB) | S>0", "Event rate (UB) | S>0")

    results %>%
        filter(Method %in% unname(methods), gen_model_type == "grf",
               surv_model_type == spec$surv,
               cens_model_type == "grf", real_data == spec$real, n_train == N_TRAIN,
               n_cal == spec$n_cal, Probability == PROB, Criterion == "low risk") %>%
        mutate(`Survival | S>0`   = ifelse(Screened > 0, Survival, NA),
               `Surv. (UB) | S>0` = ifelse(Screened > 0, Survival_upper, NA),
               `Surv. (LB) | S>0` = ifelse(Screened > 0, Survival_lower, NA),
               FDP                = 1 - Survival,
               `Event rate | S>0` = ifelse(Screened > 0, 1 - Survival, NA),
               ## bounds swap on the event scale: an upper bound on survival
               ## is a LOWER bound on the event rate
               FDP_lower          = 1 - Survival_upper,
               FDP_upper          = 1 - Survival_lower,
               `Event rate (LB) | S>0` = ifelse(Screened > 0, 1 - Survival_upper, NA),
               `Event rate (UB) | S>0` = ifelse(Screened > 0, 1 - Survival_lower, NA),
               `P(S>0)` = as.numeric(Screened > 0)) %>%
        group_by(Method, Time) %>%
        mutate(`P(S>0)` = mean(`P(S>0)`), num_sel = sum(Screened > 0),
               across(c(`Survival | S>0`, `Surv. (UB) | S>0`, `Surv. (LB) | S>0`,
                        `Event rate | S>0`,
                        `Event rate (LB) | S>0`, `Event rate (UB) | S>0`),
                      ~ ifelse(num_sel >= MIN_NUM_SEL_SHOW, .x, NA))) %>%
        ungroup() %>%
        mutate(Method = factor(Method, levels = unname(methods), labels = names(methods))) %>%
        pivot_longer(cols = any_of(names(vars)), names_to = "Variable", values_to = "Value") %>%
        filter(Variable %in% names(vars)) %>%
        mutate(hbar_y = case_when(Variable %in% surv_scale  ~ PROB,
                                  Variable %in% event_scale ~ 1 - PROB,
                                  TRUE                      ~ NA_real_),
               show_hbar = !is.na(hbar_y),
               ## y_floor is given on the survival scale; on the event scale
               ## the corresponding bound is a ceiling at 1 - y_floor
               floor_flip = Variable %in% event_scale,
               Variable = factor(Variable, levels = names(vars), labels = unname(vars)))
}

boxplot_base <- function(df, y_floor = NULL) {
    pp <- ggplot(df, aes(x = Method, y = Value, color = Method)) +
        geom_boxplot(alpha = 0.5) +
        facet_grid(Variable ~ Time, scales = "free_y",
                   labeller = labeller(Time = function(x) sprintf("Horizon: %d months", as.numeric(x)))) +
        geom_hline(data = df %>% filter(show_hbar) %>% distinct(Variable, hbar_y),
                   aes(yintercept = hbar_y), linetype = 2) +
        labs(x = "Method", y = "") +
        scale_color_manual(values = COLORS) + theme_bw() +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), legend.position = "right")

    ## Extend the survival / event-rate panels to y_floor without clipping:
    ## a blank point extends the scale, and facet_grid(scales="free_y")
    ## leaves the other rows alone.
    if (!is.null(y_floor)) {
        anchor <- df %>%
            filter(show_hbar) %>%
            distinct(Variable, Time, floor_flip) %>%
            mutate(Method = df$Method[1],
                   Value = ifelse(floor_flip, 1 - y_floor, y_floor))
        pp <- pp + geom_blank(data = anchor, aes(x = Method, y = Value), inherit.aes = FALSE)
    }
    pp
}



## ============================================================
## The figure kinds. Each returns a ggplot; the driver saves it.
## ============================================================

fig_boxplot <- function(results, spec) {
    if (empty_skip(results, sprintf("%s [%s] (no results)", spec$id, "boxplot"))) return(NULL)
    vars <- row_vars(spec$real, isTRUE(spec$small))
    methods <- spec_methods(spec)
    if (spec$real == 1) methods$Oracle <- NULL
    df <- prepare_long(results, vars, spec, methods)
    if (empty_skip(df, spec$id)) return(NULL)
    boxplot_base(df, y_floor = spec$y_floor)
}

fig_delta <- function(results, spec) {
    if (empty_skip(results, sprintf("%s [%s] (no results)", spec$id, "delta"))) return(NULL)
    vars <- row_vars(spec$real, small = TRUE)
    methods <- method_list_delta(spec_flavour(spec))
    if (spec$real == 1) methods$Oracle <- NULL
    df <- prepare_long(results, vars, spec, methods)
    if (empty_skip(df, spec$id)) return(NULL)
    boxplot_base(df, y_floor = spec$y_floor)
}

## Colour palette for an entry whose "methods" are one method at a sequence
## of parameter values. `color_value` maps each display label to a number;
## labels it returns NA for sit outside the sequence (e.g. the tuned arm)
## and are given a single contrasting colour instead, so they cannot be
## mistaken for another point on the scale.
spec_palette <- function(spec, labels) {
    v <- vapply(labels, spec$color_value, numeric(1))
    cols <- setNames(rep(spec$color_other %or% "red2", length(labels)), labels)
    idx <- which(is.finite(v))
    if (length(idx)) {
        ## Viridis runs yellow -> purple; the yellow end is too pale to read
        ## as a line and too close to the background. Take the green-to-purple
        ## stretch instead, which keeps the ordering legible throughout.
        ramp <- grDevices::hcl.colors(length(idx) + 2L, "viridis", rev = TRUE)
        cols[idx[order(v[idx])]] <- ramp[-(1:2)]
    }
    cols
}

## Line figure. `x` names the swept column: "n_cal" for the v0c setup,
## "cens_mix" for v0t. Every other dimension in the spec is held fixed;
## the swept one is deliberately NOT filtered.
##
## `color_value` in the spec replaces the standard COLORS lookup with a
## graded palette built from the labels themselves -- discrete, but ordered
## along a sequential scale, which suits a numeric parameter sweep. Labels
## mapped to NA are drawn in `color_other`, at full opacity and with their
## own plotting symbol, so the reference series stands out against the
## sweep, which is drawn faded and with a single shared symbol.
##
## `errorbars = FALSE` in the spec omits the +/- 2 SE whiskers, which are
## worth dropping when several series overlap and the bars obscure them.
fig_lines <- function(results, spec) {
    if (empty_skip(results, sprintf("%s [%s] (no results)", spec$id, "lines"))) return(NULL)
    xvar <- spec$x    %or% "n_cal"
    xlab <- spec$xlab %or% "Calibration sample size"
    vars <- row_vars(spec$real, small = TRUE)
    methods <- spec_methods(spec)
    if (spec$real == 1) methods$Oracle <- NULL
    graded <- !is.null(spec$color_value)

    if (!xvar %in% names(results)) {
        warning(sprintf("%s: column '%s' not in the results -- skipped. Were they produced by the updated experiment?", spec$id, xvar),
                call. = FALSE, immediate. = TRUE)
        return(NULL)
    }
    if (!check_horizon(results, spec)) return(NULL)

    df <- results %>%
        filter(Time == (spec$time %or% 2), Method %in% unname(methods), gen_model_type == "grf",
               surv_model_type %in% c("cox", "grf", "xgb"), n_train == N_TRAIN,
               real_data == spec$real, Probability == PROB, Criterion == "low risk")
    ## hold n_cal fixed unless it is the variable being swept
    if (xvar != "n_cal") df <- df %>% filter(n_cal == spec$n_cal)
    df <- df %>%
        mutate(xval = .data[[xvar]],
               FDP = 1 - Survival,
               `Event rate | S>0` = ifelse(Screened > 0, 1 - Survival, NA),
               FDP_lower = 1 - Survival_upper,
               FDP_upper = 1 - Survival_lower,
               `P(S>0)` = as.numeric(Screened > 0),
               Model = factor(surv_model_type, c("cox", "grf", "xgb"), c("Cox", "GRF", "XGB")),
               Method = factor(Method, levels = unname(methods), labels = names(methods))) %>%
        pivot_longer(cols = any_of(names(vars)), names_to = "Variable", values_to = "Value") %>%
        filter(Variable %in% names(vars)) %>%
        mutate(show_hbar = Variable != "Screened",
               Variable = factor(Variable, levels = names(vars), labels = unname(vars))) %>%
        group_by(Model, Method, xval, Variable, show_hbar) %>%
        summarise(SE = sd(Value)/sqrt(n()), Value = mean(Value), .groups = "drop")
    if (empty_skip(df, spec$id)) return(NULL)

    ebar <- aes(ymin = pmax(0, Value - 2*SE),
                ymax = ifelse(Variable == "Yield", Value + 2*SE, pmin(1, Value + 2*SE)))
    hbar <- df %>% filter(show_hbar) %>% distinct(Variable) %>% mutate(yintercept = 1 - PROB)
    show_eb <- isTRUE(spec$errorbars %or% TRUE)

    if (graded) {
        ## Series = "sweep" for the graded labels, "ref" for the rest. The
        ## sweep is faded and shares one symbol, so the eye reads it as a
        ## family; the reference is opaque with its own symbol.
        v  <- vapply(names(methods), spec$color_value, numeric(1))
        rf <- names(methods)[!is.finite(v)]
        df <- df %>% mutate(Series = factor(ifelse(as.character(Method) %in% rf, "ref", "sweep"),
                                            levels = c("sweep", "ref")))
        a_sweep <- spec$alpha_sweep %or% 0.45
        pp <- ggplot(df, aes(x = xval, y = Value, color = Method,
                             shape = Series, alpha = Series)) +
            geom_line(aes(group = Method)) + geom_point() +
            {if (show_eb) geom_errorbar(ebar, width = 0.1)} +
            scale_color_manual(values = spec_palette(spec, names(methods))) +
            scale_shape_manual(values = c(sweep = 16, ref = 17), guide = "none") +
            scale_alpha_manual(values = c(sweep = a_sweep, ref = 1), guide = "none") +
            guides(color = guide_legend(override.aes = list(
                shape = ifelse(names(methods) %in% rf, 17, 16),
                alpha = ifelse(names(methods) %in% rf, 1, a_sweep))))
    } else {
        pp <- ggplot(df, aes(x = xval, y = Value, color = Method, shape = Method)) +
            geom_line() + geom_point() +
            {if (show_eb) geom_errorbar(ebar, alpha = 0.5, width = 0.1)} +
            scale_color_manual(values = COLORS)
    }

    pp +
        facet_grid(Variable ~ Model, scales = "free_y",
                   labeller = labeller(Model = function(x) sprintf("Survival model: %s", x))) +
        geom_hline(data = hbar, aes(yintercept = yintercept), linetype = 2,
                   inherit.aes = FALSE) +
        labs(x = xlab, y = "", color = spec$clab %or% "Method") +
        x_scale(spec, df$xval) + theme_bw() +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), legend.position = "right")
}

## IPCW vs AIPCW, one facet per calibration rule. Always contrasts et with
## dr, so `flavour` in the spec does not apply.
##
## For the v0t setup: the censoring model is fitted on a mixture of the
## analysis censoring law and an alternative one, so moving right along the
## x axis makes Ghat progressively wrong while the analysis data is fixed.
## This is the direct read on Proposition 3.1: IPCW should drift as Ghat
## degrades, AIPCW should stay flat.
##
## Oracle and Model are omitted -- they do not use a pseudo-outcome, so they
## have no IPCW/AIPCW distinction to show.
fig_lines_flavour <- function(results, spec) {
    if (empty_skip(results, sprintf("%s [%s] (no results)", spec$id, "flavour lines"))) return(NULL)
    xvar <- spec$x    %or% "n_cal"
    xlab <- spec$xlab %or% "Calibration sample size"
    if (!xvar %in% names(results)) {
        warning(sprintf("%s: column '%s' not in the results -- skipped.", spec$id, xvar),
                call. = FALSE, immediate. = TRUE)
        return(NULL)
    }
    if (!check_horizon(results, spec)) return(NULL)
    et <- method_list("et"); dr <- method_list("dr")
    rules <- setdiff(names(et), c("Oracle", "Model"))
    map <- bind_rows(
        tibble(Method = unlist(et[rules], use.names = FALSE), Rule = rules, Flavour = "IPCW"),
        tibble(Method = unlist(dr[rules], use.names = FALSE), Rule = rules, Flavour = "AIPCW"))
    vars <- row_vars(spec$real, small = TRUE)

    df <- results %>%
        filter(Time == (spec$time %or% 2), Method %in% map$Method, gen_model_type == "grf",
               surv_model_type == spec$surv, n_train == N_TRAIN, real_data == spec$real,
               Probability == PROB, Criterion == "low risk")
    if (xvar != "n_cal") df <- df %>% filter(n_cal == spec$n_cal)
    df <- df %>%
        left_join(map, by = "Method") %>%
        mutate(xval = .data[[xvar]],
               Rule = factor(Rule, levels = rules),
               Flavour = factor(Flavour, levels = c("IPCW", "AIPCW")),
               FDP = 1 - Survival,
               FDP_lower = 1 - Survival_upper,
               FDP_upper = 1 - Survival_lower,
               `Event rate | S>0` = ifelse(Screened > 0, 1 - Survival, NA),
               `P(S>0)` = as.numeric(Screened > 0)) %>%
        pivot_longer(cols = any_of(names(vars)), names_to = "Variable", values_to = "Value") %>%
        filter(Variable %in% names(vars)) %>%
        mutate(show_hbar = Variable != "Screened",
               Variable = factor(Variable, levels = names(vars), labels = unname(vars))) %>%
        group_by(Rule, Flavour, xval, Variable, show_hbar) %>%
        summarise(SE = sd(Value)/sqrt(n()), Value = mean(Value), .groups = "drop")
    if (empty_skip(df, spec$id)) return(NULL)

    ggplot(df, aes(x = xval, y = Value, color = Flavour, shape = Flavour)) +
        geom_point() + geom_line() +
        geom_errorbar(aes(ymin = pmax(0, Value - 2*SE),
                          ymax = ifelse(Variable == "Yield", Value + 2*SE, pmin(1, Value + 2*SE))),
                      alpha = 0.5, width = 0.1) +
        facet_grid(Variable ~ Rule, scales = "free_y") +
        geom_hline(data = df %>% filter(show_hbar) %>% distinct(Variable) %>%
                       mutate(yintercept = 1 - PROB),
                   aes(yintercept = yintercept), linetype = 2) +
        labs(x = xlab, y = "", color = "Calibration", shape = "Calibration") +
        x_scale(spec, df$xval) +
        scale_color_manual(values = COLORS_FLAVOUR) + theme_bw() +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), legend.position = "right")
}

## Paired AIPCW - IPCW differences as a function of the swept variable.
## Always contrasts et with dr, so `flavour` does not apply. The SE shown is
## that of the paired difference, which is smaller than the difference of
## the two unpaired SEs whenever the flavours are positively correlated
## across replicates -- as they are here, since they share the data.
fig_lines_paired <- function(results, spec) {
    if (empty_skip(results, sprintf("%s [%s] (no results)", spec$id, "paired lines"))) return(NULL)
    xvar <- spec$x    %or% "n_cal"
    xlab <- spec$xlab %or% "Calibration sample size"
    if (!xvar %in% names(results)) {
        warning(sprintf("%s: column '%s' not in the results -- skipped.", spec$id, xvar),
                call. = FALSE, immediate. = TRUE)
        return(NULL)
    }
    if (!check_horizon(results, spec)) return(NULL)
    et <- method_list("et"); dr <- method_list("dr")
    rules <- setdiff(names(et), c("Oracle", "Model"))
    map <- bind_rows(
        tibble(Method = unlist(et[rules], use.names = FALSE), Rule = rules, Flavour = "IPCW"),
        tibble(Method = unlist(dr[rules], use.names = FALSE), Rule = rules, Flavour = "AIPCW"))
    vars <- if (spec$real == 1) c("Yield", "FDP (LB)", "FDP (UB)") else c("Yield", "FDP")
    
    ## One survival model if the spec names one, otherwise all three.
    keep_models <- if (is.null(spec$surv)) c("cox", "grf", "xgb") else spec$surv
    df <- results %>%
        filter(Time == (spec$time %or% 2), Method %in% map$Method, gen_model_type == "grf",
               surv_model_type %in% keep_models, n_train == N_TRAIN,
               real_data == spec$real, Probability == PROB, Criterion == "low risk")
    if (xvar != "n_cal") df <- df %>% filter(n_cal == spec$n_cal)
    df <- df %>%
        mutate(xval = .data[[xvar]],
               Model = factor(surv_model_type, c("cox", "grf", "xgb"), c("Cox", "GRF", "XGB"))) %>%
        left_join(map, by = "Method") %>%
        select(Seed, Model, xval, Rule, Flavour, Screened, Survival, Survival_lower, Survival_upper) %>%
        pivot_wider(id_cols = c(Seed, Model, xval, Rule), names_from = Flavour,
                    values_from = c(Screened, Survival, Survival_lower, Survival_upper)) %>%
        filter(!is.na(Screened_IPCW), !is.na(Screened_AIPCW)) %>%
        mutate(both = Screened_IPCW > 0 & Screened_AIPCW > 0,
               Yield       = Screened_AIPCW - Screened_IPCW,
               Survival    = ifelse(both, Survival_AIPCW - Survival_IPCW, NA),
               FDP         = ifelse(both, Survival_IPCW - Survival_AIPCW, NA),
               `FDP (LB)`  = ifelse(both, Survival_upper_IPCW - Survival_upper_AIPCW, NA),
               `FDP (UB)`  = ifelse(both, Survival_lower_IPCW - Survival_lower_AIPCW, NA),
               `Surv.(LB)` = ifelse(both, Survival_lower_AIPCW - Survival_lower_IPCW, NA),
               `Surv.(UB)` = ifelse(both, Survival_upper_AIPCW - Survival_upper_IPCW, NA),
               Rule = factor(Rule, levels = rules)) %>%
        select(Seed, Model, xval, Rule, any_of(vars)) %>%
        pivot_longer(cols = any_of(vars), names_to = "Variable", values_to = "Value") %>%
        filter(!is.na(Value)) %>%
        mutate(Variable = factor(Variable, levels = vars)) %>%
        mutate(Rule = factor(Rule, levels = rules)) %>%
        group_by(Model, Rule, xval, Variable) %>%
        summarise(SE = sd(Value)/sqrt(n()), Value = mean(Value), .groups = "drop")
    if (empty_skip(df, spec$id)) return(NULL)

    ## With a single survival model the Model colour carries no information,
    ## so colour by rule instead and drop the legend (the facet already
    ## names the rule).
    one_model <- length(unique(df$Model)) == 1
    pp <- ggplot(df, aes(x = xval, y = Value)) +
        geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
        facet_grid(Variable ~ Rule, scales = "free_y") +
        labs(x = xlab, y = "AIPCW - IPCW (paired by replicate)") +
        x_scale(spec, df$xval) + theme_bw() +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
    if (one_model) {
        pp + geom_point(aes(color = Rule)) + geom_line(aes(color = Rule)) +
            geom_errorbar(aes(ymin = Value - 2*SE, ymax = Value + 2*SE, color = Rule),
                          alpha = 0.5, width = 0.1) +
            scale_color_manual(values = COLORS) + theme(legend.position = "none")
    } else {
        pp + geom_point(aes(color = Model, shape = Model)) +
            geom_line(aes(color = Model)) +
            geom_errorbar(aes(ymin = Value - 2*SE, ymax = Value + 2*SE, color = Model),
                          alpha = 0.5, width = 0.1) +
            labs(color = "Survival model", shape = "Survival model") +
            scale_color_manual(values = COLORS_MODEL) + theme(legend.position = "right")
    }
}

## IPCW vs AIPCW, differenced WITHIN replicate. Always contrasts et with dr,
## so `flavour` does not apply.
##
## Both flavours run on the same replicate (same Seed, same split, same
## fitted models), so differencing removes the between-replicate variability
## that dominates side-by-side boxplots.
##
## Survival differences are defined only when BOTH flavours selected someone:
## the recorded survival of an empty selection is the convention 1, so
## differencing it is meaningless. Those replicates drop out of the survival
## panel but stay in the yield panel, where an empty selection is a real 0.
fig_paired <- function(results, spec) {
    if (empty_skip(results, sprintf("%s [%s] (no results)", spec$id, "paired diff"))) return(NULL)
    et <- method_list("et"); dr <- method_list("dr")
    rules <- setdiff(names(et), c("Oracle", "Model"))
    map <- bind_rows(
        tibble(Method = unlist(et[rules], use.names = FALSE), Rule = rules, Flavour = "IPCW"),
        tibble(Method = unlist(dr[rules], use.names = FALSE), Rule = rules, Flavour = "AIPCW"))
    vars <- if (spec$real == 1) c("Yield", "Surv.(LB)", "Surv.(UB)") else c("Yield", "Survival")

    df <- results %>%
        filter(Time > 1, Method %in% map$Method, gen_model_type == "grf",
               surv_model_type == spec$surv, n_train == N_TRAIN, real_data == spec$real,
               n_cal == spec$n_cal, Probability == PROB, Criterion == "low risk") %>%
        left_join(map, by = "Method") %>%
        select(Seed, Time, Rule, Flavour, Screened, Survival, Survival_lower, Survival_upper) %>%
        pivot_wider(id_cols = c(Seed, Time, Rule), names_from = Flavour,
                    values_from = c(Screened, Survival, Survival_lower, Survival_upper)) %>%
        filter(!is.na(Screened_IPCW), !is.na(Screened_AIPCW)) %>%
        mutate(both = Screened_IPCW > 0 & Screened_AIPCW > 0,
               Yield       = Screened_AIPCW - Screened_IPCW,
               Survival    = ifelse(both, Survival_AIPCW - Survival_IPCW, NA),
               `Surv.(LB)` = ifelse(both, Survival_lower_AIPCW - Survival_lower_IPCW, NA),
               `Surv.(UB)` = ifelse(both, Survival_upper_AIPCW - Survival_upper_IPCW, NA),
               Rule = factor(Rule, levels = rules)) %>%
        select(Seed, Time, Rule, any_of(vars)) %>%
        pivot_longer(cols = any_of(vars), names_to = "Variable", values_to = "Value") %>%
        filter(!is.na(Value)) %>%
        mutate(Variable = factor(Variable, levels = vars))
    if (empty_skip(df, spec$id)) return(NULL)

    ggplot(df, aes(x = Rule, y = Value, color = Rule)) +
        geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
        geom_boxplot(alpha = 0.5) +
        facet_grid(Variable ~ Time, scales = "free_y",
                   labeller = labeller(Time = function(x) sprintf("Horizon: %d months", as.numeric(x)))) +
        labs(x = "Method", y = "AIPCW - IPCW (paired by replicate)") +
        scale_color_manual(values = COLORS) + theme_bw() +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), legend.position = "right")
}


## ============================================================
## Driver
## ============================================================

## Kinds that show both flavours at once; their file names carry no suffix.
BOTH_FLAVOUR_KINDS <- c("paired", "paired-lines", "flavour-lines")

cat("=== figures ===\n")
for (spec in FIGURES) {
    res <- get_results(spec$setup)
    if (is.null(res)) next
    pp <- switch(spec$kind,
                 boxplot = fig_boxplot(res, spec),
                 delta   = fig_delta(res, spec),
                 lines   = fig_lines(res, spec),
                 paired  = fig_paired(res, spec),
                 "paired-lines"  = fig_lines_paired(res, spec),
                 "flavour-lines" = fig_lines_flavour(res, spec),
                 stop("unknown figure kind: ", spec$kind))
    if (is.null(pp)) next
    dir.create(OUTDIR_FIG, showWarnings = FALSE, recursive = TRUE)
    fl  <- spec_flavour(spec)
    sfx <- if (spec$kind %in% BOTH_FLAVOUR_KINDS) "" else FLAVOUR_SUFFIX[[fl]]
    fname <- file.path(OUTDIR_FIG, sprintf("%s%s.pdf", spec$file, sfx))
    ggsave(fname, pp, width = spec$w, height = spec$h, device = cairo_pdf)
    cat(sprintf("  %-12s %-38s [%s] -> %s\n", spec$id, spec$note,
                if (spec$kind %in% BOTH_FLAVOUR_KINDS) "et+dr" else fl, basename(fname)))
}


## ============================================================
## Tables
## ============================================================

summarize_table <- function(results, spec) {
    if (empty_skip(results, sprintf("%s [table] (no results)", spec$id))) return(NULL)
    methods <- method_list(spec_flavour(spec))
    if (spec$real == 1) {
        methods$Oracle <- NULL
        ## Bounds swap on the event scale: an upper bound on survival is a
        ## LOWER bound on the FDP / event rate.
        vars <- c("Screened" = "Yield", "FDP_lower" = "FDP (LB)", "FDP_upper" = "FDP (UB)",
                  "Event rate (LB) | S>0" = "Ev.rate (LB)",
                  "Event rate (UB) | S>0" = "Ev.rate (UB)",
                  "P(S>0)" = "P($|S|>0$)")
    } else {
        vars <- c("Screened" = "Yield", "FDP" = "FDP",
                  "Event rate | S>0" = "Event rate", "P(S>0)" = "P($|S|>0$)")
    }

    df.long <- results %>%
        filter(Method %in% unname(methods), gen_model_type == "grf",
               surv_model_type == spec$surv,
               cens_model_type == "grf", real_data == spec$real, n_train == N_TRAIN,
               n_cal == spec$n_cal, Probability == PROB, Criterion == "low risk") %>%
        mutate(FDP                     = 1 - Survival,
               `Event rate | S>0`      = ifelse(Screened > 0, 1 - Survival, NA),
               FDP_lower               = 1 - Survival_upper,
               FDP_upper               = 1 - Survival_lower,
               `Event rate (LB) | S>0` = ifelse(Screened > 0, 1 - Survival_upper, NA),
               `Event rate (UB) | S>0` = ifelse(Screened > 0, 1 - Survival_lower, NA),
               `P(S>0)` = as.numeric(Screened > 0),
               Method = factor(Method, levels = unname(methods), labels = names(methods))) %>%
        group_by(Time, Method) %>%
        mutate(prop_sel = mean(Screened > 0),
               across(c(`Event rate | S>0`,
                        `Event rate (LB) | S>0`, `Event rate (UB) | S>0`),
                      ~ ifelse(prop_sel >= 0.09, .x, NA))) %>%
        summarise(across(all_of(names(vars)),
                         list(mean = ~mean(., na.rm = TRUE),
                              se   = ~sd(., na.rm = TRUE)/sqrt(sum(!is.na(.)))),
                         .names = "{.col}_{.fn}"), .groups = "drop") %>%        
        mutate(across(where(is.numeric), ~ ifelse(is.nan(.), NA, .)))
    if (empty_skip(df.long, spec$id)) return(NULL)
   
    ## A method fails if its LOWER bound is significantly ABOVE the target
    ## event rate: on the survival scale this was the upper bound falling
    ## below PROB, and 1 - Survival_upper is the FDP lower bound.
    fdp_col  <- if (spec$real == 1) "FDP_lower" else "FDP"
    cond_col <- if (spec$real == 1) "Event rate (LB) | S>0" else "Event rate | S>0"
    TARGET   <- 1 - PROB

    df.long %>%
        group_by(Time) %>%
        mutate(too_high = .data[[paste0(fdp_col, "_mean")]] -
                          2 * .data[[paste0(fdp_col, "_se")]] > TARGET,
               cond_high = coalesce(.data[[paste0(cond_col, "_mean")]] -
                                    2 * .data[[paste0(cond_col, "_se")]] > TARGET, FALSE),
               best = { ok <- Method != "Oracle" & !too_high & !cond_high
                        v <- Screened_mean[ok]
                        if (!length(v) || all(is.na(v))) NA_real_ else max(v, na.rm = TRUE) }) %>%
        ungroup() %>%
        mutate(best = Screened_mean == best & Method != "Oracle") %>%
        pivot_longer(cols = ends_with("_mean") | ends_with("_se"),
                     names_to = c("Variable", ".value"), names_pattern = "(.*)_(mean|se)") %>%
        filter(Variable %in% names(vars)) %>%
        mutate(Variable = factor(Variable, levels = names(vars), labels = unname(vars)),
               formatted = case_when(
                   is.na(mean) ~ "--",
                   Method == "Oracle" ~ sprintf("%.0f (%.2f)", mean, 2*se),
                   Variable == "Yield" & !is.na(best) & best & round(mean, 2) > 0 ~
                       sprintf("\\textbf{%.0f} (%.2f)", mean, 2*se),
                   Variable == "Yield" ~ sprintf("%.0f (%.2f)", mean, 2*se),
                   (Variable %in% c("FDP (LB)", "FDP") & too_high) |
                   (Variable %in% c("Ev.rate (LB)", "Event rate") & cond_high) ~
                       sprintf("\\textcolor{red}{%.2f (%.2f)}", mean, 2*se),
                   TRUE ~ sprintf("%.2f (%.2f)", mean, 2*se))) %>%
        select(Time, Method, Variable, formatted) %>%
        pivot_wider(names_from = Variable, values_from = formatted) %>%
        arrange(Time, Method)
}

## Model-based exceedance for one group of replicates.
##
##     X_i | N_i ~ Binom(N_i, theta_i),   theta_i ~ Beta(a, b)
##     reported:  P(theta < target) = pbeta(target, a, b)
##
## aod parametrises by (mu, phi), mu = a/(a+b), phi = 1/(a+b+1), so
## a + b = 1/phi - 1. See ?aod::betabin and Prentice (1986, JASA 81, 321).
##
## No interval is reported. As phi -> 0 the fitted Beta collapses to a point
## mass at mu and the exceedance becomes a STEP function of mu at the target,
## so every way of propagating uncertainty in (mu, phi) -- delta method,
## marginal rectangle, percentile bootstrap -- returns something either
## degenerate or vacuous. The point estimate is a useful summary; a
## companion interval would not be.
##
## NA is returned when the fit fails or the data cannot identify the two
## parameters (fewer than five non-empty replicates, or no variation in the
## observed proportions).
estimate.under.survival <- function(res, target = 0.9) {
    res <- res %>% filter(Screened > 0, !is.na(Survival))
    if (nrow(res) < 1) return(NA_real_)
    d <- res %>% mutate(N = as.integer(Screened), X = as.integer(round(Screened * Survival))) %>%
                 filter(N > 0, X >= 0, X <= N)
    if (nrow(d) < 5 || length(unique(d$X / d$N)) < 2) return(NA_real_)

    f <- try(suppressWarnings(aod::betabin(cbind(X, N - X) ~ 1, ~ 1, data = d)), silent = TRUE)
    if (inherits(f, "try-error")) return(NA_real_)
    s <- try(summary(f), silent = TRUE)
    if (inherits(s, "try-error")) return(NA_real_)
    mu  <- tryCatch(stats::plogis(s@Coef[1, "Estimate"]), error = function(e) NA_real_)
    phi <- tryCatch(s@Phi[1, "Estimate"],                 error = function(e) NA_real_)

    if (!is.finite(mu) || !is.finite(phi) || phi >= 1) return(NA_real_)
    if (phi <= 0) return(as.numeric(mu < target))
    k <- 1/phi - 1
    stats::pbeta(target, mu * k, (1 - mu) * k)
}

## Shared body of both exceedance tables.
##
##   methods    named method list, as method_list() or method_list_delta()
##   delta_for  function of the DISPLAY label returning the delta that row
##              promises, or NA for a row that makes no such claim
summarize_exceedance_core <- function(results, spec, methods, delta_for,
                                      split_delta = FALSE) {
    if (empty_skip(results, sprintf("%s [exceedance] (no results)", spec$id))) return(NULL)
    if (spec$real == 1) {
        warning(sprintf("%s: needs the oracle outcome, which real-data runs do not have -- skipped.",
                        spec$id), call. = FALSE, immediate. = TRUE)
        return(NULL)
    }

    ## ONE survival model -- three made the table too wide. `surv` in the
    ## spec picks it; GRF by default.
    keep_models <- spec$surv %or% "grf"

    ## Same filter as prepare_long(), so this table and the corresponding
    ## boxplot figure describe exactly the same replicates.
    df <- results %>%
        filter(Time > 1, Method %in% unname(methods), gen_model_type == "grf",
               surv_model_type %in% keep_models,
               n_train == N_TRAIN, real_data == spec$real, n_cal == spec$n_cal,
               Probability == PROB, Criterion == "low risk") %>%
        mutate(Method = factor(Method, levels = unname(methods), labels = names(methods)))
    ## optional horizon subset, so a long table can be split across files
    if (!is.null(spec$times)) df <- df %>% filter(Time %in% spec$times)
    if (empty_skip(df, spec$id)) return(NULL)
    if (all(is.na(df$Survival))) {
        warning(sprintf("%s: Survival is all NA -- skipped.", spec$id),
                call. = FALSE, immediate. = TRUE)
        return(NULL)
    }

    stats <- df %>%
        group_by(Time, Method) %>%
        group_modify(~ {
            fdp_v  <- 1 - .x$Survival
            ev_v   <- 1 - .x$Survival[.x$Screened > 0]
            se     <- function(v) { v <- v[!is.na(v)]
                                    if (length(v) < 2) NA_real_ else sd(v)/sqrt(length(v)) }
            tibble(
                fdp      = mean(fdp_v, na.rm = TRUE),  fdp_se    = se(fdp_v),
                evrate   = if (mean(.x$Screened > 0) < MIN_NONEMPTY_EXC) NA_real_
                           else mean(ev_v,  na.rm = TRUE),  evrate_se = se(ev_v),
                yield    = mean(.x$Screened),          yield_se  = se(.x$Screened),
                nonempty = mean(.x$Screened > 0),      nonempty_se = se(as.numeric(.x$Screened > 0)),
                ## no standard error for the exceedance: it comes from a
                ## beta-binomial fit whose uncertainty cannot be propagated
                ## sensibly (see estimate.under.survival)
                exc      = if (mean(.x$Screened > 0) < MIN_NONEMPTY_EXC) NA_real_
                           else estimate.under.survival(.x, target = PROB))
        }) %>%
        ungroup()

    failed <- sum(is.na(stats$exc))
    if (failed)
        cat(sprintf("    note: exceedance not identified in %d of %d cells\n",
                    failed, nrow(stats)))

    stats %>%
        ## Two decimals: on a few dozen replicates of a few dozen patients
        ## the third decimal is noise.
        ##
        ## RED marks a quantity above the level it is supposed to respect:
        ## the event rate above the target alpha, and the estimated
        ## exceedance above the delta that row promises. Both are point
        ## comparisons, not tests -- with no interval they should be read as
        ## flags rather than as evidence of a violation. Rows making no
        ## high-probability claim -- Oracle, Model, and FDR-Conformal, which
        ## controls E[FDP] instead -- get NA from delta_for() and are never
        ## flagged on the exceedance.
        mutate(dlt = vapply(as.character(Method), delta_for, numeric(1)),
               hl  = is.finite(dlt) & !is.na(exc) & round(exc, 2) > dlt,
               hl_ev = is.finite(evrate) & round(evrate, 3) > (1 - PROB) & as.character(Method) != "Oracle") %>%
        ## bold the largest yield among methods not flagged on either the
        ## event rate or the exceedance, excluding the Oracle
        group_by(Time) %>%
        mutate(best = { ok <- as.character(Method) != "Oracle" & !hl & !hl_ev
                        v <- yield[ok]
                        if (!length(v) || all(is.na(v))) NA_real_ else max(v, na.rm = TRUE) },
               best = is.finite(best) & yield > 0 & yield == best & as.character(Method) != "Oracle") %>%
        ungroup() %>%
        mutate(Yield = ifelse(is.nan(yield), "--",
                       ifelse(best & round(yield, 2) > 0,
                              sprintf("\\textbf{%.0f} (%.2f)", yield, 2*yield_se),
                              sprintf("%.0f (%.2f)", yield, 2*yield_se))),
               FDP     = ifelse(is.nan(fdp),   "--", sprintf("%.3f (%.3f)", fdp, 2*fdp_se)),
               EvRate  = case_when(
                   !is.finite(evrate) ~ "--",
                   hl_ev              ~ sprintf("\\textcolor{red}{%.3f (%.3f)}", evrate, 2*evrate_se),
                   TRUE               ~ sprintf("%.3f (%.3f)", evrate, 2*evrate_se)),
               Exc     = case_when(
                   is.na(exc) ~ "--",
                   hl         ~ sprintf("\\textcolor{red}{%.2f}", exc),
                   TRUE       ~ sprintf("%.2f", exc)),
               Nonempty = sprintf("%.2f (%.2f)", nonempty, 2*nonempty_se),
               ## short rule name and its own delta, taken from the RAW
               ## label before any markup is added to it
               Rule  = sub("-[0-9.]+$", "", as.character(Method)),
               Delta = ifelse(is.finite(dlt), sprintf("%.2f", dlt), "--"),
               ord   = ifelse(is.finite(dlt), dlt, -1),
               ## the method name is coloured too, so the left-hand column
               ## alone says which methods ever failed
               Method = as.character(Method)) %>%
        { if (!split_delta) {
              arrange(., Time) %>%
                  select(Time, Method = Rule, Yield, FDP, EvRate, Exc, Nonempty)
          } else {
              ## Ordered by delta so the three rules sit together at each
              ## level; rows with no delta claim (Oracle, FDR-Conformal)
              ## come first, once per horizon.
              ##
              ## delta leads the row and is printed ONLY on the first of
              ## each block. Repeating it beside all three rules turned the
              ## column into a wall of the same four numbers and hid where
              ## one level ended and the next began.
              arrange(., Time, ord, Rule) %>%
                  group_by(Time, Delta) %>%
                  mutate(Delta = ifelse(dplyr::row_number() == 1L, Delta, "")) %>%
                  ungroup() %>%
                  select(Time, Delta, Method = Rule, Yield, FDP, EvRate, Exc, Nonempty)
          } }
}

## Companion to Figure 1: the standard method list, all at the single
## confidence level DELTA.
summarize_exceedance <- function(results, spec) {
    summarize_exceedance_core(results, spec, method_list(spec_flavour(spec)),
                              function(m) if (grepl("^HP-", m)) DELTA else NA_real_)
}

## Companion to Figure 2: HP-LTT at several confidence levels. Each row is
## judged against ITS OWN delta, read back out of the label the method list
## built with sprintf("HP-LTT-%g", ds).
summarize_exceedance_delta <- function(results, spec) {
    ds <- spec$deltas %or% c(0.05, 0.1, 0.2, 0.5)
    summarize_exceedance_core(results, spec,
                              method_list_delta_table(spec_flavour(spec), ds = ds),
                              function(m) suppressWarnings(as.numeric(sub("^HP-[^-]+-", "", m))),
                              ## a lone delta makes the leading column a
                              ## single repeated value, so use the plain layout
                              split_delta = length(ds) > 1)
}

## LaTeX for the exceedance table. One survival model, five value columns:
##   Yield        mean number selected, over ALL replicates
##   FDP          mean FDP over ALL replicates (matches Fig 1)
##   Event rate   mean event rate over non-empty replicates only
##   P(r>alpha)   beta-binomial exceedance estimate; red when it exceeds
##                that row's delta
##   P(|S|>0)     share of replicates that selected anyone
write_exceedance_latex <- function(tab, path) {
    body  <- tab %>% select(-Time)
    sizes <- tab %>% count(Time, name = "n")
    ## five columns for the single-delta table, six when a Delta column is
    ## present (the delta sweep), which keeps the method labels short
    names(body) <- if ("Delta" %in% names(body))
        c("$\\delta$", "Method", "Yield", "FDP", "Event rate", "$P(r>\\alpha)$", "$P(|S|>0)$")
    else
        c("Method", "Yield", "FDP", "Event rate", "$P(r>\\alpha)$", "$P(|S|>0)$")
    n_left <- if ("$\\delta$" %in% names(body)) 2L else 1L
    kb <- body %>%
        kbl(escape = FALSE, booktabs = TRUE, format = "latex",
            align = c(rep("l", n_left), rep("r", ncol(body) - n_left)))
    start <- 1L
    for (i in seq_len(nrow(sizes))) {
        end <- start + sizes$n[i] - 1L
        kb <- kb %>% group_rows(group_label = sprintf("Horizon: %s months", sizes$Time[i]),
                                start_row = start, end_row = end)
        start <- end + 1L
    }
    dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
    cat(as.character(kb), file = path)
}

write_table_latex <- function(tab, path) {
    body  <- tab %>% select(-Time)
    sizes <- tab %>% count(Time, name = "n")
    kb <- body %>% kbl(escape = FALSE, booktabs = TRUE, format = "latex",
                       align = c("l", rep("r", ncol(body) - 1)))
    start <- 1L
    for (i in seq_len(nrow(sizes))) {
        end <- start + sizes$n[i] - 1L
        kb <- kb %>% group_rows(group_label = sprintf("Horizon: %s months", sizes$Time[i]),
                                start_row = start, end_row = end)
        start <- end + 1L
    }
    dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
    cat(as.character(kb), file = path)
}

cat("\n=== tables ===\n")
for (spec in TABLES) {
    res <- get_results(spec$setup)
    if (is.null(res)) next
    kind <- spec$kind %or% "summary"
    tab <- switch(kind,
                  summary            = summarize_table(res, spec),
                  exceedance         = summarize_exceedance(res, spec),
                  "exceedance-delta" = summarize_exceedance_delta(res, spec),
                  stop("unknown table kind: ", kind))
    if (is.null(tab)) next
    fl   <- spec_flavour(spec)
    path <- file.path(OUTDIR_TAB, sprintf("%s%s.txt", spec$file, FLAVOUR_SUFFIX[[fl]]))
    if (kind %in% c("exceedance", "exceedance-delta")) write_exceedance_latex(tab, path)
    else write_table_latex(tab, path)
    cat(sprintf("  %-12s %-38s [%s] -> %s\n", spec$id, spec$note, fl, basename(path)))
}

cat("\nDefault FLAVOUR =", FLAVOUR,
    "-- set `flavour = \"dr\"` (or \"ft\") on any entry to override it.\n")
cat("If an output is missing, run  report_methods(get_results(\"v0\"), \"v0\")\n")
