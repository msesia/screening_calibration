#'Evaluate Screening Selections Without Oracle Event Times
#'
#' Computes conservative bounds on the survival rate of selected patients
#' without access to true event times. This function evaluates selection performance
#' under two extreme assumptions for censored data:
#' - lower bound assumes censored patients failed at the censoring time
#' - upper bound assumes censored patients never fail (event time = ∞)
#'
#' @param data.test A data frame containing the test dataset. Must include columns `time` and `status`.
#' @param idx.selected A vector of indices for selected patients (row indices into `data.test`).
#' @param screening_time Numeric, the time horizon for screening (e.g., 5 months).
#' @param screening_prob Numeric between 0 and 1; the survival threshold used for classification.
#' @param screening_crit A string, either `"low risk"` or `"high risk"`, indicating the direction of the screening rule.
#'
#' @return A data frame with the following columns:
#' \describe{
#'   \item{Screening.time}{The time horizon used for screening}
#'   \item{Screening.rule}{A textual description of the rule applied}
#'   \item{Screened}{Number of selected patients}
#'   \item{Proportion.survived.lower}{Conservative estimate assuming censored events occurred at censoring time}
#'   \item{Proportion.survived.upper}{Liberal estimate assuming censored events never occurred}
#' }
#'
#' @seealso \code{\link{evaluate_selections}}, \code{\link{select_patients_band}}
#' @export
evaluate_selections_without_oracle <- function(data.test, idx.selected, screening_time, screening_prob, screening_crit) {
  stopifnot(screening_crit %in% c("low risk", "high risk"))

  # Create oracle versions of the data
  data.test.oracle.lower <- data.test
  data.test.oracle.lower$event_time <- data.test$time  # always use time, regardless of censoring

  data.test.oracle.upper <- data.test
  data.test.oracle.upper$event_time <- ifelse(data.test$status == 1, data.test$time, Inf)

  # Evaluate selections using lower oracle
  res.lower <- evaluate_selections(data.test.oracle.lower, idx.selected, screening_time, screening_prob, screening_crit)
  res.lower$Proportion.survived.lower <- res.lower$Proportion.survived
  res.lower$Proportion.survived <- NULL
  res.lower$Valid <- NULL

  # Evaluate selections using upper oracle
  res.upper <- evaluate_selections(data.test.oracle.upper, idx.selected, screening_time, screening_prob, screening_crit)
  res.upper$Proportion.survived.upper <- res.upper$Proportion.survived
  res.upper$Proportion.survived <- NULL
  res.upper$Valid <- NULL

  # Merge the two results by common keys
  merge(res.lower, res.upper, by = c("Screening.time", "Screening.rule", "Screened"))
}


#' Evaluate Screening Rule Against True or Imputed Event Times
#'
#' Given a set of selected patients, compute the proportion who survive beyond a given time,
#' and determine whether this satisfies the conservative screening criterion.
#'
#' @param data.test.oracle A data frame that must contain a column `event_time`, representing the known or imputed event time.
#' @param idx.selected Integer vector of row indices indicating selected patients.
#' @param screening_time Numeric value representing the time horizon for evaluation.
#' @param screening_prob Numeric value between 0 and 1, used as a threshold for the screening rule.
#' @param screening_crit Character string; either `"low risk"` (expecting high survival) or `"high risk"` (expecting low survival).
#'
#' @return A named list with the screening time, rule description, number of patients screened,
#'         observed survival proportion, and whether the result satisfies the conservative screening rule.
#' @export
evaluate_selections <- function(data.test.oracle, idx.selected, screening_time, screening_prob, screening_crit) {
  stopifnot(screening_crit %in% c("low risk", "high risk"))

  is_alive <- data.test.oracle$event_time >= screening_time
  n <- nrow(data.test.oracle)

  num_selections <- length(idx.selected)
  if (num_selections > 0) {
    selected_alive <- sum(is_alive[idx.selected])
    selected_dead <- num_selections - selected_alive
    prop_selected_alive <- selected_alive / num_selections
  } else {
    selected_alive <- 0
    selected_dead <- 0
    prop_selected_alive <- if (screening_crit == "low risk") 1 else 0
  }

  if (screening_crit == "low risk") {
    conservative <- prop_selected_alive > screening_prob
  } else {
    conservative <- prop_selected_alive < screening_prob
  }

  # Fallback in degenerate case
  if (is.na(conservative)) conservative <- TRUE

  rule_desc <- interpret_screening_rule(screening_time, screening_prob, screening_crit)

  result <- as.data.frame(list(
    Screening.time = screening_time,
    Screening.rule = rule_desc,
    Screened = num_selections,
    Proportion.survived = prop_selected_alive,
    Valid = conservative
  ))

  return(result)
}


#' Format Screening Rule as String
#'
#' @param time Screening time.
#' @param prob Survival probability threshold.
#' @param crit Screening direction (`"low risk"` or `"high risk"`).
#'
#' @return A formatted string describing the rule.
#' @export
interpret_screening_rule <- function(time, prob, crit) {
    stopifnot(crit %in% c("low risk", "high risk"))
    stopifnot(time>=0)
    stopifnot(prob>=0 & prob<=1)
    if(crit=="high risk") {
        interpretation <- sprintf("high-risk patients with P(survival at time %.2f) < %.2f", time, prob)
    } else {
        interpretation <- sprintf("low-risk patients with P(survival at time %.2f) > %.2f", time, prob)
    }
    return(interpretation)
}
