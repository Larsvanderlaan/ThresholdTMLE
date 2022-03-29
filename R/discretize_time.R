
#' Takes a dataset with continuous time-to-event variables and target times at which to estimate the threshold-response cumulative incidence
#' and returns a new dataset with the time-to-event variables discretized into integer values and the new target integer times that correspond with the continuous target times.
#' @param data A data.frame containing the observed data.
#' Its columns should include the covariates for which to adjust, the treatment variable,
#' the time-to-event-or-censoring variable `Ttilde := min(T, C)`, the indicator of whether the survival event is observed `Delta := 1(C <= T)`.
#' @param Ttilde A string giving the column name of the time-to-event-or-censoring variable `Ttilde := min(T, C)`.
#' @param biased_sampling_indicator A string giving the column name of the indicator variable of whether or not the observation's treatment is observed.
#' The indicator should take the value `1` if the treatment is observed and `0` otherwise.
#' This parameter combined with the weight parameter \code{weights_var} allows this method to be used with biased sampling designs.
#' @param target_times A numeric vector of time-points at which to estimate the cumulative incidence.
#' @param nbins_t The number of bins to use for discretizing time. The recommended value is 20-30.
#' @export
discretize_time <- function(data, Ttilde, biased_sampling_indicator = NULL, target_times, nbins_t = 25) {
  if(is.null(biased_sampling_indicator)) {
    R <- 1
  } else {
    R <- data[[biased_sampling_indicator]]
  }
  # The maximum reference time of interest
  upper_t <- max(target_times)
  # Discretize continuous times
  # Number of bins to use in time discretization

  # The discrete grid of times (equally spaced on quantile scale)
  time_grid <- unique(sort(quantile(data[["Ttilde"]][data[["Ttilde"]] <= upper_t + 1 & R ==1  ], seq(0,1, length = nbins_t))))
  # Make sure target times appear in time grid
  time_grid <- sort(union(time_grid, target_times))

  Ttilde_discrete <- findInterval(data[["Ttilde"]], time_grid, all.inside = TRUE)
  target_times_discrete <- findInterval(target_times, time_grid, all.inside = TRUE)
  names(target_times_discrete) <- paste0("t=",target_times )
  data$Ttilde <- Ttilde_discrete
  out<- list(data_discrete = data, target_times_discrete = target_times_discrete)
  return(out)
}
