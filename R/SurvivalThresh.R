
#' @param data A data.frame or data.table containing the data.
#' This should contain all data including those with missing "A" variable.
#' @param covariates A vector of covariate names in \code{data} (i.e. the "W" node)
#' @param trt Treatment variable name (i.e. the "A" node)
#' @param Ttilde Time-to-event or censoring variable Ttilde  (i.e. the "Ttilde" node)
#' @param Delta Variable name for indicator of being not censored  (i.e. the "Delta" node)
#' @param J Variable name for mark/type of event  (i.e. the "J" node)
#' @param weights Variable name storing the IPW weights for the missing treatment/"A" values. Can be NULL.
#' If \code{biased_sampling_indicator} and \code{biased_sampling_group} are specifyed then the weights will be estimated using NPMLE.
#' @param biased_sampling_indicator
#' @param biased_sampling_group
#' @param cutoffs_A
#' @param cutoffs_J
#' @param target_times
#' @param lrnr
#' @param lrnr_A
#' @param lrnr_C
#' @param lrnr_N
#' @param lrnr_J
#' @param lrnr_A
#' @param ngrid_A
#' @param type_A
#' @param type_J
#' @param split_by_J
#' @param max_eps
#' @param max_iter
#' @param fast_analysis
#' @param verbose
#' @export
survivalThresh <- function(data, covariates, trt = "A", Ttilde = "Ttilde", Delta = "Delta", J = "J", biased_sampling_indicator = NULL,  biased_sampling_group = NULL, weights = NULL, cutoffs_A, cutoffs_J, target_times, lrnr = Lrnr_glmnet$new(), lrnr_A = lrnr, lrnr_C = lrnr, lrnr_N = lrnr, lrnr_J = lrnr, ngrid_A = 25, type_A = c("above", "below", "equal"), type_J = c("above", "below", "equal"), split_by_J = TRUE, max_eps = 0.05, max_iter = 100, fast_analysis = F, verbose = TRUE) {
  data <- as.data.table(data)
  n_full_sample <- nrow(data)
  data_full <- data
  print(weights)
  if(is.null(weights) & !is.null(biased_sampling_group)) {
    groups <- unique(data_full[[biased_sampling_group]])
    data_full$weights <- 0

    for(grp in groups) {
      keep <- data_full[[biased_sampling_group]] == grp

      set(data_full, which(keep), "weights", rep(1/mean(data_full[keep, biased_sampling_indicator, with = F][[1]]), sum(keep)))
    }
    weights <- "weights"
  }
  if(verbose) {
    print("Processing data...")
  }
  processed <- process_data(data_full, covariates, trt, Ttilde, Delta, J, weights)
  data_full <- processed$data
  data_full[[biased_sampling_indicator]] <- data[[biased_sampling_indicator]]
  data_full[[biased_sampling_group]] <- data[[biased_sampling_group]]


  node_list <- processed$node_list
  if(!is.null(biased_sampling_indicator)) {
    data <- data_full[data_full[[biased_sampling_indicator]]==1,]
  }
  print(data_full)
  print(data)
  type_J <- match.arg(type_J)
  type_A <- match.arg(type_A)

  data_orig <- data
  # Shape data to long format (n*t rows)
  data <- shape_long(data, node_list)
  # Assumes data starts at t =1
  times <- 1:max(data_orig$Ttilde)
  nt <- length(times)
  n <- nrow(data_orig)
  if(type_A=="equal") {
    grid <- sort(unique(cutoffs_A))
  } else {
    grid <-as.vector(quantile(data_orig$A, seq(0,1, length = ngrid_A)))
    donothing <- sapply(cutoffs_A, function(c) {
      grid[which.min(abs(c - grid))[1]] <<-c

    })
  }
  # Fit initial nuisance estimates
  if(verbose) {
    print("Fitting likelihood...")
  }

  fits <- fit_likelihood(data, node_list, grid, cutoffs_A, cutoffs_J, lrnr = lrnr, lrnr_A = lrnr_A, lrnr_C = lrnr_C, lrnr_N = lrnr_N, lrnr_J = lrnr_J, type_A = type_A, type_J = type_J, verbose = verbose)
  # Get observed data likelihoods
  if(verbose) {
    print("Computing likelihoods...")
  }
  likelihoods <- get_likelihoods(fits, data, node_list)
  # Target the J-distribution
  if(verbose) {
    print("Targeting distribution of mark")
  }
  update_list <- target_J(likelihoods, fits, data, target_times, node_list)
  likelihoods <- update_list$likelihoods
  fits <- update_list$fits
  out <- list(likelihoods = likelihoods, fits = fits)
  colMeans(out$fits$EIC_J)
  # Target the overall survival (N) distribution
  if(verbose) {
    print("Targeting survival hazard")
  }
  out_orig <- out
  out_final <- out
  if(split_by_J) {

    res_list <- list()
    for(j in seq_along(cutoffs_J)) {


      out <- out_orig

      out$likelihoods$J <- out$likelihoods$J[,j, with = F]
      out$likelihoods$outcomes_J <- out$likelihoods$outcomes_J[,j, with = F]
      out$fits$cutoffs_J <- cutoffs_J[j]


      out$fits$EIC_J <- as.matrix(out$fits$EIC_J[,((j-1)*length(target_times)*length(cutoffs_A)+1):((j)*length(target_times)*length(cutoffs_A)),drop =  F])
      out$fits$epsilons_J <- out$fits$epsilons_J[j]

      for(i in 1:max_iter){
        print(i)
        out <- target_N(data, out$likelihoods, out$fits, node_list, target_times = target_times, max_eps = max_eps, n_full_sample = n_full_sample)
        if(out$converged) {
          break
        }


        if(i==85) {
          out$fits$max_eps <- out$fits$max_eps/1.5
        }

      }

      if(!out$converged) {
        out <- target_N(data, out$likelihoods, out$fits, node_list, target_times = target_times, max_eps = max_eps, n_full_sample = n_full_sample, force_converge = TRUE)
      }

      if(j==1) {
        out_final$fits <- list(out$fits)
        out_final$likelihoods <- list(out$likelihoods)
      } else {
        out_final$fits[[j]] <- (out$fits)
        out_final$likelihoods[[j]] <- (out$likelihoods)
        }

      likelihoods <- out$likelihoods
      fits <- out$fits
      # Perform sequential regression to account for stochastic intervention targeting
      if(type_A!="equal") {
        if(verbose) {
          print("Computing parameter estimates with sequential regression")
        }

        res <- sequential_targeting(data, data_orig, data_full, fits, likelihoods, target_times, node_list, n_full_sample , biased_sampling_group , biased_sampling_indicator , fast_analysis = fast_analysis)

      } else if (type_A == "equal") {
        if(verbose) {
          print("Computing parameter estimates")
        }
        res <- Fixed_treatment_targeting(data, data_orig, fits, likelihoods, target_times, node_list)

      }
      res_list[[j]] <- res
    }
    res <- res_list
    names(res) = paste0("J=",cutoffs_J)
    estimates_t <- list()
    for(index_t in seq_along(target_times)) {
      estimates_J <-list()
      for(index_j in seq_along(cutoffs_J)) {
        estimates_J[[paste0("J=",cutoffs_J[index_j])]] <- res[[index_j]][[index_t]][[1]]
      }
      estimates_t[[paste0("t=", target_times[index_t])]] <- estimates_J
    }
    res <-estimates_t
    out <- out_final
    names(out$fits) <- paste0("J=", cutoffs_J)
    names(out$likelihoods) <- paste0("J=", cutoffs_J)
  } else {


    for(i in 1:max_iter){
      print(i)
      out <- target_N(data, out$likelihoods, out$fits, node_list, target_times = target_times, max_eps = max_eps, n_full_sample = n_full_sample)
      if(out$converged) {
        break
      }
      if(i==50) {
        out$fits$max_eps <- out$fits$max_eps/1.3
      }

      if(i==85) {
        out$fits$max_eps <- out$fits$max_eps/1.3
      }
    }
    if(!out$converged) {
      out <- target_N(data, out$likelihoods, out$fits, node_list, target_times = target_times, max_eps = max_eps, n_full_sample = n_full_sample, force_converge = TRUE)
    }
    likelihoods <- out$likelihoods
    fits <- out$fits
    # Perform sequential regression to account for stochastic intervention targeting
    if(type_A!="equal") {
      if(verbose) {
        print("Computing parameter estimates with sequential regression")
      }
      res <- sequential_targeting(data, data_orig, data_full, fits, likelihoods, target_times, node_list, n_full_sample = n_full_sample, fast_analysis = fast_analysis)
    } else if (type_A == "equal") {
      if(verbose) {
        print("Computing parameter estimates")
      }
      res <- Fixed_treatment_targeting(data, data_orig, data_full, fits, likelihoods, target_times, node_list, n_full_sample = n_full_sample, fast_analysis = fast_analysis)

    }



  }
  likelihoods <- out$likelihoods
  fits <- out$fits
  return(list(estimates = res, likelihoods = likelihoods, fits = fits, data = data))


}







process_data <- function(data, covariates, trt = "A", Ttilde = "Ttilde", Delta = "Delta", J = "J", weights = NULL) {
  node_list <- list(W = covariates, A = "A", weights = "weights", "J" = "J", Delta = "Delta", Ttilde = "Ttilde")
  data_new <- data.table( A = data[[trt]], Ttilde = data[[Ttilde]], Delta = data[[Delta]], J = data[[J]])
  set(data_new,, covariates, data[, covariates, with = F])

  if(!is.null(weights)) {
    data_new$weights <- data[[weights]]
  } else {
    data_new$weights <- 1
  }
  data <- data_new
  data$id <- as.factor(1:nrow(data))
  rm(data_new)
  return(list(data = data, node_list = node_list))
}

