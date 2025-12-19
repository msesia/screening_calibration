options(width = 300)

library(tidyverse)
library(kableExtra)

interpret_screening_rule <- function(prob, crit, latex=FALSE) {
    stopifnot(crit %in% c("low risk", "high risk"))
    stopifnot(prob>=0 & prob<=1)
    if(crit=="high risk") {
        if(latex) {
            interpretation <- sprintf("high-risk patients with P($T > t$)$ < %.2f$", prob)
        } else {
            interpretation <- sprintf("high-risk patients with P(T > t) < %.2f", prob)
        }
    } else {
        if(latex) {
            interpretation <- sprintf("low-risk patients with P($T > t$)$ > %.2f$", prob)
        } else {
            interpretation <- sprintf("low-risk patients with P(T> t) > %.2f", prob)
        }    }
    return(interpretation)
}

get_custom_colors_model <- function() {
    custom_colors <- c(
        "Model" = "black",
        "KM" = "forestgreen",
        "CSB" = "red",
        "Boot (LTT)" = "purple",
        "Boot (greedy)" = "steelblue"
    )
}
get_custom_shapes_model <- function() {
    custom_shapes <- c(
        "Model" = 16,
        "KM" = 18,
        "CSB" = 8,
        "Boot (LTT)" = 17,
        "Boot (greedy)" = 2
    )
}

load_results <- function(setup) {
    idir <- sprintf("results_hpc/%s", setup)
    ifile.list <- list.files(idir)
    results <- do.call("rbind", lapply(ifile.list, function(ifile) {
        df <- read_delim(sprintf("%s/%s", idir, ifile), delim=",", col_types=cols(), guess_max=2)
    }))
    return(results)
}

####################
## SEMI-SYNTHETIC ##
####################

make_figures_boxplots <- function(results, method.list, colors.list,
                                  plot.real = 1, plot.gen_model = "grf", plot.surv_model = "grf",
                                  plot.n_train = 1000, plot.n_cal = 1000, prob.screen = 0.9,
                                  small.figure = FALSE, save_fig = FALSE) {
  stopifnot(is.list(method.list), !is.null(names(method.list)))

  ## custom_colors <- get_custom_colors_model()
  ## custom_shapes <- get_custom_shapes_model()

  if (plot.real == 1) {
      method.list$Oracle <- NULL
      if(small.figure) {
          variable.values <- c("Screened", "Survival_lower", "Survival_upper")
          variable.labels <- c("Yield", "Surv.(LB)", "Surv.(UB)")
      } else {
          variable.values <- c("Screened", "Survival_lower", "Survival_upper", "Surv. (LB) | S>0", "Surv. (UB) | S>0", "P(S>0)")
          variable.labels <- c("Yield", "Surv.(LB)", "Surv.(UB)", "Cond.Surv.(LB)", "Cond.Surv.(UB)", "P(non-empty)")
      }          
  } else {
      if(small.figure) {
          variable.values <- c("Screened", "Survival")
          variable.labels <- c("Yield", "Survival")
      } else {
          variable.values <- c("Screened", "Survival", "Survival | S>0", "P(S>0)")
          variable.labels <- c("Yield", "Survival", "Survival (cond.)", "P(non-empty)")
      }
  }
  method.values <- unname(method.list)
  method.labels <- names(method.list)

  if (plot.real == 1) {
      if(small.figure) {
          df.range <- tibble(Value=c(0,1,0,1), Method="Model",
                             Variable=c("Survival_lower", "Survival_lower", "Survival_upper", "Survival_upper"))
      } else {
          df.range <- tibble(Value=c(0,1,0,1,0,1,0,1), Method="Model",
                             Variable=c("Survival_lower", "Survival_lower", "Survival_upper", "Survival_upper",
                                        "Surv. (UB) | S>0", "Surv. (UB) | S>0", "Surv. (LB) | S>0", "Surv. (LB) | S>0"))
      }
  } else {
      if(small.figure) {
          df.range <- rbind(tibble(Method="Model", Value=c(0,1), Variable="Survival"))
      } else {
          df.range <- rbind(tibble(Method="Model", Value=c(0,1), Variable="Survival"),
                            tibble(Method="Model", Value=c(0,1), Variable="Survival | S>0"))
      }
  }
  df.range <- df.range %>% mutate(Variable = factor(Variable, levels = variable.values, labels = variable.labels))     
  
  df <- results %>%
      filter(
          Time > 1,
          Method %in% method.values,
          gen_model_type == plot.gen_model,
          surv_model_type == plot.surv_model,
          n_train==plot.n_train,
          real_data == plot.real,
          n_cal == plot.n_cal,
          Probability == prob.screen,
          Criterion == "low risk"
      ) %>%
      mutate(
          `Survival | S>0` = ifelse(Screened>0, Survival, NA),
          `Surv. (UB) | S>0` = ifelse(Screened>0, Survival_upper, NA),
          `Surv. (LB) | S>0` = ifelse(Screened>0, Survival_lower, NA),
          `P(S>0)` = Screened>0) %>%
      group_by(Method, Time) %>%
      mutate(`P(S>0)` = mean(`P(S>0)`),
             num_sel = sum(Screened>0),
             `Survival | S>0` = ifelse(num_sel>=9, `Survival | S>0`, NA),
             `Surv. (UB) | S>0` = ifelse(num_sel>=9, `Surv. (UB) | S>0`, NA),
             `Surv. (LB) | S>0` = ifelse(num_sel>=9, `Surv. (LB) | S>0`, NA)
             ) %>%
      ungroup() %>%
      mutate(Method = factor(Method, levels = method.values, labels = method.labels)) %>%
      select(Method, Time, Screened, Survival, `Survival | S>0`, Survival_lower, Survival_upper, `Surv. (LB) | S>0`, `Surv. (UB) | S>0`, `P(S>0)`) %>%
      pivot_longer(
          cols = c(Screened, Survival, `Survival | S>0`, Survival_lower, Survival_upper, `Surv. (LB) | S>0`, `Surv. (UB) | S>0`, `P(S>0)`),
          names_to = c("Variable"),
          values_to = c("Value")
      ) %>%
      filter(Variable %in% variable.values) %>%
      mutate(
          show_hbar = Variable %in% c("Survival", "Surv. (UB) | S>0", "Surv. (LB) | S>0", "Survival | S>0", "Survival_lower", "Survival_upper"),
          Variable = factor(Variable, levels = variable.values, labels = variable.labels)
      )         
      
  title.str <- sprintf("Experiments with %s data. Screening for: %s.\nNum train: %d, num cal: %d, num test: 1000.\nGen model: %s, Surv model: %s",
                       ifelse(plot.real, "REAL", "SEMI-SYNTHETIC"),
                       interpret_screening_rule(prob.screen, "low risk"),
                       plot.n_train, plot.n_cal, plot.gen_model, plot.surv_model)

  pp <- ggplot(df, aes(x = Method, y = Value, color = Method, shape = Method)) +
      geom_boxplot(alpha=0.5) +
      facet_grid(
          Variable ~ Time,
          scales = "free_y",
          labeller = labeller(
              Time = function(x) sprintf("Horizon: %d months", as.numeric(x))
          )
      ) +
      geom_point(data=df.range, alpha=0) +
      geom_hline(
          data = df %>% filter(show_hbar) %>% distinct(Variable) %>% mutate(yintercept = prob.screen),
          aes(yintercept = yintercept),
          linetype = 2
      ) +
      labs(
          x = "Method",
          y = "",
          #subtitle = title.str
      ) +
      scale_color_manual(values = colors.list) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position="right")      

  if (save_fig) {
      if(small.figure) {
          fig.height <- ifelse(plot.real, 8, 3.2)
          fname <- sprintf("figures/boxplots_%s_gen_%s_surv_%s_ntrain%s_ncal%s_small.pdf",
                           plot.real, plot.gen_model, plot.surv_model, plot.n_train, plot.n_cal)
      } else {
      fig.height <- ifelse(plot.real, 8, 5.2)
          fname <- sprintf("figures/boxplots_%s_gen_%s_surv_%s_ntrain%s_ncal%s.pdf",
                           plot.real, plot.gen_model, plot.surv_model, plot.n_train, plot.n_cal)
      }          
      ggsave(fname, pp, width = 7, height = fig.height)
      cat(sprintf("Saved figure in: %s.\n", fname))      
  } else {
      return(pp)
  }
}

ipcw_delta <- 0.1
ipcw_style <- "ft"
method.list <- list(
  "Oracle" = "oracle (group)",
  "Model" = "model (group)",
#  "Model (CI)" = "model (ci, group)",
#  "Conformal" = "conformal",
  "FDR-Conformal" = "conformal (stable)",
  "HP-Greedy" = sprintf("%s|%.2f|pt_delta", ipcw_style, ipcw_delta),
  "HP-Uniform" = sprintf("%s|%.2f|uniform", ipcw_style, ipcw_delta),
  "HP-LTT" = sprintf("%s|%.2f|ltt_delta", ipcw_style, ipcw_delta)
#  "Delta (LTT, 0.5)" = sprintf("LTT-0.5-ipcw%d-ipcw_delta", ipcw_style),
#  "Delta (greedy, 0.5)" = sprintf("LTT-0.5-greedy-ipcw%d-ipcw_delta", ipcw_style)
)
colors.list <- c(
  "Oracle" = "orange1",
  "Model" = "grey50",
#  "Model (CI)" = "grey20",
#  "Conformal" = "red",
  "FDR-Conformal" = "red2",
  "HP-Greedy" = "darkgreen",
  "HP-Uniform" = "green",
  "HP-LTT" = "blue2"
#  "Delta (LTT, 0.5)" = "steelblue2",
#  "Delta (greedy, 0.5)" = "seagreen1"
)

results <- load_results("v0")
save_fig <- TRUE

for(small.figure in c(FALSE,TRUE)) {
    for(plot.surv_model in c("grf", "cox", "xgb")) {
        for(plot.n_cal in c(100, 1000)) {
            make_figures_boxplots(results, method.list, colors.list, plot.real=0, plot.gen_model = "grf", plot.surv_model = plot.surv_model,
                                  plot.n_train = 5000, plot.n_cal=plot.n_cal, small.figure=small.figure, save_fig=save_fig)
        }
    }
}


results <- load_results("v1")
for(plot.surv_model in c("grf", "cox", "xgb")) {
    for(plot.n_cal in c(100, 1000)) {
        make_figures_boxplots(results, method.list, colors.list, plot.real=1, plot.gen_model = "grf", plot.surv_model = plot.surv_model,
                              plot.n_train = 5000, plot.n_cal=plot.n_cal, save_fig=save_fig)
    }
}

summarize_table <- function(results, method.list, colors.list,
                            plot.real = 1, plot.cens_model = "grf", plot.surv_model = "grf",
                            plot.n_train = 2000, plot.n_cal = 1000, prob.screen = 0.9,
                            save_latex = TRUE) {
    stopifnot(is.list(method.list), !is.null(names(method.list)))

    if (plot.real == 1) {
        method.list$Oracle <- NULL
        variable.values <- c("Screened", "Survival_lower", "Survival_upper", "Surv. (LB) | S>0", "Surv. (UB) | S>0", "P(S>0)")
        variable.labels <- c("Yield", "Surv.(LB)", "Surv.(UB)", "C.Surv.(LB)", "C.Surv.(UB)", "P($|S|>0$)")
    } else {
        variable.values <- c("Screened", "Survival", "Survival | S>0", "P(S>0)")
        variable.labels <- c("Yield", "Survival", "Survival (cond.)", "P($|S|>0$)")
    }
    method.values <- unname(method.list)
    method.labels <- names(method.list)

    df.long <- results %>%
        filter(
            Method %in% method.values,
            surv_model_type == plot.surv_model,
            cens_model_type == plot.cens_model,
            real_data == plot.real,
            n_train == plot.n_train,
            n_cal == plot.n_cal,
            Probability == prob.screen,
            Criterion == "low risk"
        ) %>%
        mutate(
            `Survival | S>0` = ifelse(Screened>0, Survival, NA),
            `Surv. (UB) | S>0` = ifelse(Screened>0, Survival_upper, NA),
            `Surv. (LB) | S>0` = ifelse(Screened>0, Survival_lower, NA),
            `P(S>0)` = Screened>0) %>%
        mutate(Method = factor(Method, levels = method.values, labels = method.labels)) %>%
        group_by(Time, Method, Criterion, Probability) %>%
        mutate(
            prop_sel = mean(Screened>0),
            `Survival | S>0` = ifelse(prop_sel>=0.09, `Survival | S>0`, NA),
            `Surv. (UB) | S>0` = ifelse(prop_sel>=0.09, `Surv. (UB) | S>0`, NA),
            `Surv. (LB) | S>0` = ifelse(prop_sel>=0.09, `Surv. (LB) | S>0`, NA)
        ) %>%
        summarise(
            across(
                c(Screened, Survival, "Survival | S>0", Survival_lower, Survival_upper, "Surv. (LB) | S>0", "Surv. (UB) | S>0", "P(S>0)"),
                list(mean = ~mean(., na.rm = TRUE), se = ~sd(., na.rm = TRUE) / sqrt(n())),
                .names = "{.col}_{.fn}"
            ),
            .groups = "drop"
        ) %>%
        ungroup() %>%
        mutate(across(where(is.numeric), ~ ifelse(is.nan(.), NA, .)))


    df <- df.long %>%
        pivot_longer(
            cols = ends_with("_mean") | ends_with("_se"),
            names_to = c("Variable", ".value"),
            names_pattern = "(.*)_(mean|se)"
        ) %>%
        filter(Variable %in% variable.values) %>%
        mutate(
            Variable = factor(Variable, levels = variable.values, labels = variable.labels),
            formatted = sprintf("%.2f (%.2f)", mean, 2 * se)
        ) %>%
        select(Time, Method, Variable, formatted) %>%
        pivot_wider(names_from = Variable, values_from = formatted) %>%
        arrange(Time, Method)

    ## ------- FORMATTED TABLE (bold max yield; red too-low survival & cond. surv. UB) -------
    if (plot.real) {
        df.formatted <- df.long %>%
            group_by(Time) %>%
            mutate(
                Survival_too_low   = Survival_upper_mean + 2 * Survival_upper_se < prob.screen,
                CondSurv_too_low   = (`Surv. (UB) | S>0_mean` + 2 * `Surv. (UB) | S>0_se`) < prob.screen,
                CondSurv_too_low = ifelse(is.na(CondSurv_too_low), FALSE, CondSurv_too_low),
                Screened_max = {
                    ok   <- Method != "Oracle" & !Survival_too_low & !CondSurv_too_low 
                    vals <- Screened_mean[ok]
                    if (length(vals) == 0 || all(is.na(vals))) NA_real_ else max(vals, na.rm = TRUE)
                }
            ) %>%
            ungroup()
    } else {
        df.formatted <- df.long %>%
            group_by(Time) %>%
            mutate(
                Survival_too_low   = Survival_mean + 2 * Survival_se < prob.screen,
                CondSurv_too_low   = (`Survival | S>0_mean` + 2 * `Survival | S>0_se`) < prob.screen,
                CondSurv_too_low = ifelse(is.na(CondSurv_too_low), FALSE, CondSurv_too_low),
                Screened_max = {
                    ok   <- Method != "Oracle" & !Survival_too_low & !CondSurv_too_low 
                    vals <- Screened_mean[ok]
                    if (length(vals) == 0 || all(is.na(vals))) NA_real_ else max(vals, na.rm = TRUE)
                }                
            ) %>%
            ungroup()
    }
    
    df.formatted <- df.formatted %>%
        mutate(Screened_max = (Screened_mean == Screened_max & Method != "Oracle")) %>%
        pivot_longer(
            cols = ends_with("_mean") | ends_with("_se"),
            names_to = c("Variable", ".value"),
            names_pattern = "(.*)_(mean|se)"
        ) %>%
        filter(Variable %in% variable.values) %>%
        mutate(
            Variable = factor(Variable, levels = variable.values, labels = variable.labels),
            formatted = ifelse(is.na(mean),
                               "NA",
                        ifelse(
                            ## Oracle: always plain
                            Method == "Oracle", sprintf("%.0f (%.2f)", mean, 2 * se),
                            ## Non-Oracle
                        ifelse(
                            Variable == "Yield",
                        ifelse(!is.na(Screened_max) & (Screened_max == 1)==TRUE & round(mean,2)>0,
                               sprintf("\\textbf{%.0f} (%.2f)", mean, 2 * se),
                               sprintf("%.0f (%.2f)", mean, 2 * se)),
                        ifelse(
                        (Variable %in% c("Surv.(UB)", "Survival")        & Survival_too_low) |
                        (Variable %in% c("C.Surv.(UB)", "Survival (cond.)") & CondSurv_too_low),
                        sprintf("\\textcolor{red}{%.2f (%.2f)}", mean, 2 * se),
                        sprintf("%.2f (%.2f)", mean, 2 * se)
                        )
                        )
                        )
                        )
        ) %>%
        select(Time, Method, Variable, formatted) %>%
        pivot_wider(names_from = Variable, values_from = formatted) %>%
        arrange(Time, Method)        
    

                                        # After you build df.formatted and do arrange(Time, Method):

    if (save_latex) {
        dir.create("tables", showWarnings = FALSE, recursive = TRUE)
        file_path <- sprintf("tables/data_%s_ntrain_%s_ncal_%s_%s.txt",
                             plot.real, plot.n_train, plot.n_cal, plot.surv_model)

        if (!requireNamespace("kableExtra", quietly = TRUE)) {
            stop("Please install the 'kableExtra' package to use LaTeX output.")
        }
        library(kableExtra)
        library(dplyr)

                                        # 1) Table body without the Time column
        df_body <- df.formatted %>%
            arrange(Time, Method) %>%
            select(-Time)

                                        # 2) Group sizes in the same order as the table body
        time_sizes <- df.formatted %>%
            arrange(Time, Method) %>%
            count(Time, name = "n")

                                        # 3) Build the base kable (header appears once)
        kb <- df_body %>%
            kbl(
                escape = FALSE, booktabs = TRUE, format = "latex",
                align = c("l", rep("r", ncol(df_body) - 1))
            )

                                        # 4) Add bold group headers per Time without repeating table headers
        start <- 1L
        for (i in seq_len(nrow(time_sizes))) {
            n_i <- time_sizes$n[i]
            end <- start + n_i - 1L
            lbl <- sprintf("Horizon: %s months", time_sizes$Time[i])
            kb <- kb %>% group_rows(group_label = lbl, start_row = start, end_row = end)
            start <- end + 1L
        }

                                        # 5) Write LaTeX to file
        cat(as.character(kb), file = file_path)
        cat(sprintf("Saved grouped table (by Time) in: %s.\n", file_path))
    }

    return(list(table_raw = df, table_formatted = df.formatted))
}

results <- load_results("v0")
for(plot.surv_model in c("grf", "cox", "xgb")) {
    for(plot.n_cal in c(100,1000)) {
        summarize_table(results, method.list, colors.list, plot.real=0, plot.cens_model = "grf", plot.surv_model = plot.surv_model,
                        plot.n_train = 5000, plot.n_cal=plot.n_cal, save_latex=TRUE)
    }
}


results <- load_results("v1")
for(plot.surv_model in c("grf", "cox", "xgb")) {
    for(plot.n_cal in c(100,1000)) {
        summarize_table(results, method.list, colors.list, plot.real=1, plot.cens_model = "grf", plot.surv_model = plot.surv_model,
                        plot.n_train = 5000, plot.n_cal=plot.n_cal, save_latex=TRUE)
    }
}


#############################################
## Plot results vs calibration sample size ##
#############################################

make_figures_lines <- function(results, method.list, colors.list,
                               plot.real = 1, plot.gen_model = "grf",
                               plot.n_train = 1000, prob.screen = 0.9,
                               small.figure = TRUE, save_fig = FALSE) {
    stopifnot(is.list(method.list), !is.null(names(method.list)))

    ## custom_colors <- get_custom_colors_model()
    ## custom_shapes <- get_custom_shapes_model()

    if (plot.real == 1) {
        method.list$Oracle <- NULL
        if(small.figure) {
            variable.values <- c("Screened", "Survival_lower", "Survival_upper")
            variable.labels <- c("Yield", "Surv.(LB)", "Surv.(UB)")
        } else {
            variable.values <- c("Screened", "Survival_lower", "Survival_upper", "Surv. (LB) | S>0", "Surv. (UB) | S>0", "P(S>0)")
            variable.labels <- c("Yield", "Surv.(LB)", "Surv.(UB)", "Cond.Surv.(LB)", "Cond.Surv.(UB)", "P(non-empty)")
        }          
    } else {
        if(small.figure) {
            variable.values <- c("Screened", "Survival")
            variable.labels <- c("Yield", "Survival")
        } else {
            variable.values <- c("Screened", "Survival", "Survival | S>0", "P(S>0)")
            variable.labels <- c("Yield", "Survival", "Survival (cond.)", "P(non-empty)")
        }
    }
    method.values <- unname(method.list)
    method.labels <- names(method.list)

    range.ymin <- 0.8
    
    if (plot.real == 1) {
        if(small.figure) {
            df.range <- tibble(Value=c(range.ymin,1,range.ymin,1), Method="Model", n_cal=1000,
                               Variable=c("Survival_lower", "Survival_lower", "Survival_upper", "Survival_upper"))
        } else {
            df.range <- tibble(Value=c(range.ymin,1,range.ymin,1,range.ymin,1,range.ymin,1), Method="Model",n_cal=1000,
                               Variable=c("Survival_lower", "Survival_lower", "Survival_upper", "Survival_upper",
                                          "Surv. (UB) | S>0", "Surv. (UB) | S>0", "Surv. (LB) | S>0", "Surv. (LB) | S>0"))
        }
    } else {
        if(small.figure) {
            df.range <- rbind(tibble(Method="Model", Value=c(range.ymin,1), Variable="Survival", n_cal=1000))
        } else {
            df.range <- rbind(tibble(Method="Model", Value=c(range.ymin,1), Variable="Survival", n_cal=1000),
                              tibble(Method="Model", Value=c(range.ymin,1), Variable="Survival | S>0",n_cal=1000))
        }
    }
    df.range <- df.range %>% mutate(Variable = factor(Variable, levels = variable.values, labels = variable.labels))     
    
    df <- results %>%
        filter(
            Time %in% c(2),
            Method %in% method.values,
            gen_model_type == plot.gen_model,
            surv_model_type %in% c("cox", "grf", "xgb"),
            n_train==plot.n_train,
            real_data == plot.real,
            Probability == prob.screen,
            Criterion == "low risk"
        ) %>%
        mutate(
            `Survival | S>0` = ifelse(Screened>0, Survival, NA),
            `Surv. (UB) | S>0` = ifelse(Screened>0, Survival_upper, NA),
            `Surv. (LB) | S>0` = ifelse(Screened>0, Survival_lower, NA),
            `P(S>0)` = Screened>0,
             Model = factor(surv_model_type, c("cox", "grf", "xgb"), c("Cox", "GRF", "XGB")) ) %>%
        group_by(Model, Method, n_cal) %>%
        mutate(`P(S>0)` = mean(`P(S>0)`),
               num_sel = sum(Screened>0),
               `Survival | S>0` = ifelse(num_sel>=9, `Survival | S>0`, NA),
               `Surv. (UB) | S>0` = ifelse(num_sel>=9, `Surv. (UB) | S>0`, NA),
               `Surv. (LB) | S>0` = ifelse(num_sel>=9, `Surv. (LB) | S>0`, NA)
               ) %>%
        ungroup() %>%
        mutate(Method = factor(Method, levels = method.values, labels = method.labels)) %>%
        select(Model, Method, n_cal, Seed, Screened, Survival, `Survival | S>0`, Survival_lower, Survival_upper, `Surv. (UB) | S>0`, `Surv. (LB) | S>0`, `P(S>0)`) %>%
        pivot_longer(
            cols = c(Screened, Survival, `Survival | S>0`, Survival_upper, Survival_lower, `Surv. (UB) | S>0`, `Surv. (LB) | S>0`, `P(S>0)`),
            names_to = c("Variable"),
            values_to = c("Value")
        ) %>%
        filter(Variable %in% variable.values) %>%
        mutate(
            show_hbar = Variable %in% c("Survival", "Surv. (UB) | S>0", "Surv. (LB) | S>0", "Survival | S>0", "Survival_lower", "Survival_upper"),
            Variable = factor(Variable, levels = variable.values, labels = variable.labels)
        ) %>%
        group_by(Model, Method, n_cal, Variable, show_hbar) %>%
        summarise(SE = sd(Value)/sqrt(n()), Value = mean(Value))
    
    pp <- ggplot(df, aes(x = n_cal, y = Value, color = Method, shape = Method)) +
        geom_point(alpha=1) +
        geom_line(alpha=1) +
        geom_errorbar(aes(ymin=pmax(0,Value-2*SE), ymax=ifelse(Variable=="Yield", pmin(1000,Value+2*SE), pmin(1,Value+2*SE))),
                      alpha=0.5, width=0.1) +
        facet_grid(
            Variable~Model,
            scales = "free_y",
            labeller = labeller(
                Time = function(x) sprintf("Horizon: %d months", as.numeric(x)),
                Model = function(x) sprintf("Survival model: %s", x)                
            )
        ) +
        geom_point(data=df.range, alpha=0) +
        geom_hline(
            data = df %>% filter(show_hbar) %>% distinct(Variable) %>% mutate(yintercept = prob.screen),
            aes(yintercept = yintercept),
            linetype = 2
        ) +
        labs(
            x = "Calibration sample size",
            y = ""
            ) +
        scale_x_continuous(trans='log10') +        
        scale_color_manual(values = colors.list) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position="right")      

    if (save_fig) {
        if(small.figure) {
            fig.height <- ifelse(plot.real, 8, 3.2)
            fname <- sprintf("figures/lines_%s_gen_%s_ntrain%s_small.pdf",
                             plot.real, plot.gen_model, plot.n_train)
        } else {
            fig.height <- ifelse(plot.real, 8, 5.2)
            fname <- sprintf("figures/lines_%s_gen_%s_ntrain%s.pdf",
                             plot.real, plot.gen_model, plot.n_train)
        }          
        ggsave(fname, pp, width = 7, height = fig.height)
        cat(sprintf("Saved figure in: %s.\n", fname))      
    } else {
        return(pp)
    }
}

ipcw_delta <- 0.1
ipcw_style <- "ft"
method.list <- list(
  "Oracle" = "oracle (group)",
  "Model" = "model (group)",
#  "Model (CI)" = "model (ci, group)",
#  "Conformal" = "conformal",
  "FDR-Conformal" = "conformal (stable)",
  "HP-Greedy" = sprintf("%s|%.2f|pt_delta", ipcw_style, ipcw_delta),
  "HP-Uniform" = sprintf("%s|%.2f|uniform", ipcw_style, ipcw_delta),
  "HP-LTT" = sprintf("%s|%.2f|ltt_delta", ipcw_style, ipcw_delta)
#  "Delta (LTT, 0.5)" = sprintf("LTT-0.5-ipcw%d-ipcw_delta", ipcw_style),
#  "Delta (greedy, 0.5)" = sprintf("LTT-0.5-greedy-ipcw%d-ipcw_delta", ipcw_style)
)
colors.list <- c(
  "Oracle" = "orange1",
  "Model" = "grey50",
#  "Model (CI)" = "grey20",
#  "Conformal" = "red",
  "FDR-Conformal" = "red2",
  "HP-Greedy" = "darkgreen",
  "HP-Uniform" = "green",
  "HP-LTT" = "blue2"
#  "Delta (LTT, 0.5)" = "steelblue2",
#  "Delta (greedy, 0.5)" = "seagreen1"
)

results <- load_results("v0c")
save_fig <- TRUE

for(small.figure in c(TRUE)) {
    for(plot.surv_model in c("grf", "xgb")) {
        make_figures_lines(results, method.list, colors.list, plot.real=0, plot.gen_model = "grf", plot.n_train = 5000, small.figure=small.figure, save_fig=save_fig)
    }
}



############################################
## Plot LTT for different values of delta ##
############################################

make_figures_boxplots_delta <- function(results, method.list, colors.list,
                                        plot.real = 1, plot.gen_model = "grf", plot.surv_model = "grf",
                                        plot.n_train = 1000, plot.n_cal = 1000, prob.screen = 0.9,
                                        save_fig = FALSE) {
    stopifnot(is.list(method.list), !is.null(names(method.list)))

  ## custom_colors <- get_custom_colors_model()
  ## custom_shapes <- get_custom_shapes_model()

  if (plot.real == 1) {
      method.list$Oracle <- NULL
    variable.values <- c("Screened", "Survival_lower", "Survival_upper")
    variable.labels <- c("Yield", "Surv.(LB)", "Surv.(UB)")
  } else {
    variable.values <- c("Screened", "Survival")
    variable.labels <- c("Yield", "Survival")
  }
  method.values <- unname(method.list)
  method.labels <- names(method.list)

  if (plot.real == 1) {
      df.range <- tibble(Value=c(0,1,0,1), Method="FDR-Conformal",
                         Variable=c("Survival_lower", "Survival_lower", "Survival_upper", "Survival_upper"))
  } else {
      df.range <- rbind(tibble(Method="Oracle", Value=c(0,1), Variable="Survival"))
  }
  df.range <- df.range %>% mutate(Variable = factor(Variable, levels = variable.values, labels = variable.labels))     
  
  df <- results %>%
      filter(
          Time > 1,
          Method %in% method.values,
          gen_model_type == plot.gen_model,
          surv_model_type == plot.surv_model,
          n_train==plot.n_train,
          real_data == plot.real,
          n_cal == plot.n_cal,
          Probability == prob.screen,
          Criterion == "low risk"
      ) %>%
      mutate(
          `Survival | S>0` = ifelse(Screened>0, Survival, NA),
          `Surv. (UB) | S>0` = ifelse(Screened>0, Survival_upper, NA),
          `Surv. (LB) | S>0` = ifelse(Screened>0, Survival_lower, NA),
          `P(S>0)` = Screened>0) %>%
      group_by(Method, Time) %>%
      mutate(`P(S>0)` = mean(`P(S>0)`),
             num_sel = sum(Screened>0),
             `Survival | S>0` = ifelse(num_sel>=9, `Survival | S>0`, NA),
             `Surv. (UB) | S>0` = ifelse(num_sel>=9, `Surv. (UB) | S>0`, NA),
             `Surv. (LB) | S>0` = ifelse(num_sel>=9, `Surv. (LB) | S>0`, NA)
             ) %>%
      ungroup() %>%
      mutate(Method = factor(Method, levels = method.values, labels = method.labels)) %>%
      select(Method, Time, Screened, Survival, `Survival | S>0`, Survival_lower, Survival_upper, `Surv. (UB) | S>0`, `Surv. (LB) | S>0`, `P(S>0)`) %>%
      pivot_longer(
          cols = c(Screened, Survival, `Survival | S>0`, Survival_upper, Survival_lower, `Surv. (UB) | S>0`, `Surv. (LB) | S>0`, `P(S>0)`),
          names_to = c("Variable"),
          values_to = c("Value")
      ) %>%
      filter(Variable %in% variable.values) %>%
      mutate(
          show_hbar = Variable %in% c("Survival", "Surv. (UB) | S>0", "Surv. (LB) | S>0", "Survival | S>0", "Survival_lower", "Survival_upper"),
          Variable = factor(Variable, levels = variable.values, labels = variable.labels)
      )         
      
  title.str <- sprintf("Experiments with %s data. Screening for: %s.\nNum train: %d, num cal: %d, num test: 1000.\nGen model: %s, Surv model: %s",
                       ifelse(plot.real, "REAL", "SEMI-SYNTHETIC"),
                       interpret_screening_rule(prob.screen, "low risk"),
                       plot.n_train, plot.n_cal, plot.gen_model, plot.surv_model)

  pp <- ggplot(df, aes(x = Method, y = Value, color = Method, shape = Method)) +
      geom_boxplot(alpha=0.5) +
      facet_grid(
          Variable ~ Time,
          scales = "free_y",
          labeller = labeller(
              Time = function(x) sprintf("Horizon: %d months", as.numeric(x))
          )
      ) +
      geom_point(data=df.range, alpha=0) +
      geom_hline(
          data = df %>% filter(show_hbar) %>% distinct(Variable) %>% mutate(yintercept = prob.screen),
          aes(yintercept = yintercept),
          linetype = 2
      ) +
      labs(
          x = "Method",
          y = "",
          #subtitle = title.str
      ) +
      scale_color_manual(values = colors.list) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position="right")      

  if (save_fig) {
      fig.height <- ifelse(plot.real, 4, 3.2)
      fname <- sprintf("figures/boxplots_%s_gen_%s_surv_%s_ntrain%s_ncal%s_delta.pdf", plot.real, plot.gen_model, plot.surv_model, plot.n_train, plot.n_cal)
      ggsave(fname, pp, width = 7, height = fig.height)
      cat(sprintf("Saved figure in: %s.\n", fname))      
  } else {
      return(pp)
  }
}


ipcw_delta <- 0.1
ipcw_style <- "ft"
method.list <- list(
  "Oracle" = "oracle (group)",
  "FDR-Conformal" = "conformal (stable)",
  "HP-LTT-0.05" = sprintf("%s|%.2f|ltt_delta", ipcw_style, 0.05),
  "HP-LTT-0.1" = sprintf("%s|%.2f|ltt_delta", ipcw_style, 0.1),
  "HP-LTT-0.2" = sprintf("%s|%.2f|ltt_delta", ipcw_style, 0.2),
  "HP-LTT-0.5" = sprintf("%s|%.2f|ltt_delta", ipcw_style, 0.5)
)
colors.list <- c(
  "Oracle" = "orange1",
  "FDR-Conformal" = "red2",
  "HP-LTT-0.05" = "blue2",
  "HP-LTT-0.1" = "dodgerblue3",
  "HP-LTT-0.2" = "deepskyblue2",
  "HP-LTT-0.5" = "skyblue1"
)

results <- load_results("v0")
save_fig <- TRUE

for(plot.surv_model in c("grf", "cox", "xgb")) {
    for(plot.n_cal in c(100, 1000)) {
        make_figures_boxplots_delta(results, method.list, colors.list, plot.real=0, plot.gen_model = "grf", plot.surv_model = plot.surv_model,
                                    plot.n_train = 5000, plot.n_cal=plot.n_cal, save_fig=save_fig)
    }
}


results <- load_results("v1")

for(plot.surv_model in c("grf", "cox", "xgb")) {
    for(plot.n_cal in c(100, 1000)) {
        make_figures_boxplots_delta(results, method.list, colors.list, plot.real=1, plot.gen_model = "grf", plot.surv_model = plot.surv_model,
                                    plot.n_train = 5000, plot.n_cal=plot.n_cal, save_fig=save_fig)
    }
}
