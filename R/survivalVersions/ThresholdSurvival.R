
#' Estimation and inference for the threshold-response cumulative incidence.
#' Let `T` be the time-until-event variable and `C` be the time-until-censoring. The observed data-structure is `(W, A, Ttilde, Delta)` where `W` is a vector of covariates,
#' `A` is a continuous treatment or biomarker, `Ttilde := min(T, C)` is the minimum of the censoring and event times,
#' `Delta := `(T <= C)` is the indicator of whether the event is observed.
#' For a given threshold/cutoff `v` of the treatment `A` and reference time point `t`,
#' the estimand of interest is given by `E_W[P(T <= t | A >= v, W)]`.
#' The estimand can be estimated using the below functions for multiple values of `t` and `v`.
#' The argument `cutoffs_A` specifies the threshold of `A` at which to estimate the estimand.
#' The argument `target_times` specifies the time-points at which to estimate the estimand.
#'
#' NOTE: This function only supports discrete time-to-event variables that begin counting at `t=1`.
#' For continuous times, one can discretize the time into bins and then apply this method to the discretized time-to-event dataset.
#' Having a very large grid of times (e.g. <=20-30) may make this function run slow.
#'
#' @param data A data.frame containing the observed data.
#' Its columns should include the covariates for which to adjust, the treatment variable,
#' the time-to-event-or-censoring variable `Ttilde := min(T, C)`, the indicator of whether the survival event is observed `Delta := 1(C <= T)`.
#' @param covariates A character vector of covariate names for which to adjust.
#' @param trt A string giving the column name of the continuous treatment variable.
#' @param Ttilde A string giving the column name of the time-to-event-or-censoring variable `Ttilde := min(T, C)`.
#' This time-to-event variable should be discrete and start counting at `t=1`.
#' @param Delta A string giving the column name of the indicator of observing the event variable `Delta := 1(T<=C)`.
#' @param biased_sampling_indicator A string giving the column name of the indicator variable of whether or not the observation's treatment is observed.
#' The indicator should take the value `1` if the treatment is observed and `0` otherwise.
#' This parameter combined with the weight parameter \code{weights_var} allows this method to be used with biased sampling designs.
#' @param weights_var (Optional) A string giving the column name of the weight variable.
#' For instance, this variable could correspond with IPW weights used to adjust for biased sampling or treatment missingness.
#' @param cutoffs_A A numeric vector of treatment values at which the threshold-response cumulative incidence should be estimated.
#' @param target_times A numeric vector of time-points at which to estimate the cumulative incidence.
#' @param lrnr_A A \code{tlverse/sl3} \code{sl3_Learner} object for binary outcomes that is used to estimate the conditional density of the treatment using discretization and pooled logistic regression.
#' @param lrnr_C A \code{tlverse/sl3} \code{sl3_Learner} object for binary outcomes that is used to estimate the conditional censoring hazard function using pooled logistic regression.
#' @param lrnr_T A \code{tlverse/sl3} \code{sl3_Learner} object for binary outcomes that is used to estimate the conditional event hazard function using pooled logistic regression.
#' @param ngrid_A (Internal use) The number of bins used when discretizing the continuous treatment for density estimation.
#' @param monotone_decreasing Whether the monotone corrected estimates of the threshold-response cumulative incidence should be monotone decreasing or monotone increasing.
#' @export
ThresholdSurvival <- function(data, covariates, trt = "A", Ttilde = "Ttilde", Delta = "Delta",   biased_sampling_indicator = NULL, weights_var, cutoffs_A, target_times,  lrnr_A = Lrnr_hal9001$new(smoothness_orders = 0, max_degree = 2, num_knots = c(30,15)) , lrnr_C = Lrnr_hal9001$new(smoothness_orders = 0, max_degree = 2, num_knots = c(30,15)) , lrnr_T = Lrnr_hal9001$new(smoothness_orders = 0, max_degree = 2, num_knots = c(30,15)) ,  ngrid_A = 25,  monotone_decreasing = T) {

  data$J123 <- 1
  suppressWarnings({suppressMessages({survivalThresh(data, covariates, trt = trt, Ttilde = Ttilde, Delta = Delta, J = "J123", biased_sampling_indicator = biased_sampling_indicator,  biased_sampling_group = NULL, weights_var = weights_var, cutoffs_A = cutoffs_A, cutoffs_J = 1, target_times = target_times, lrnr = lrnr, lrnr_A = lrnr_A, lrnr_C = lrnr_C, lrnr_N = lrnr_T, lrnr_J = Lrnr_glm$new(), ngrid_A = 25, type_J = c("above", "below", "equal"), max_eps = 0.25, max_iter = 50,  verbose = FALSE, monotone_decreasing = T)
  })})
}
