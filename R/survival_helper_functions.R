
#' Wald-Style Confidence Intervals
#'
#' @importFrom stats qnorm
#'
wald_ci <- function(est, se, level = 0.95, q = NULL) {
  if (is.null(q)) {
    q <- abs(stats::qnorm(p = (1 - level) / 2))
  }

  ci_low <- est - q * se
  ci_high <- est + q * se
  return(cbind(ci_low, ci_high))
}


bound <- function(x, bounds) {
  lower <- bounds[[1]]
  if (length(bounds) > 1) {
    upper <- bounds[[2]]
  } else {
    upper <- 1 - lower
  }
  pmin(pmax(x, lower), upper)
}


shape_long <- function(data, node_list, times = 1:max(data$Ttilde)) {
  times <- 1:max(data$Ttilde)
  t_mat <- matrix(rep(times, nrow(data)), ncol = length(times), byrow = T)
  Tt <- 1*(t_mat >= data$Ttilde)
  Nt <- data$Delta*Tt
  at_risk  <- 1*(t_mat <= data$Ttilde)
  Ct  <- (1-data$Delta)*Tt
  #at_risk_Ct <- at_risk_Nt*(1-Nt)
  data_long <- data.table(t = as.vector(t_mat), Nt = as.vector(Nt), Ct = as.vector(Ct), at_risk = as.vector(at_risk),id = rep(1:nrow(data), length(times)), Event = as.vector(1*(t_mat == data$Ttilde)))
  set(data_long, , node_list$W, data.table(apply(data[,node_list$W, with = F],2, rep, length(times))))
  set(data_long, , "A", rep(data[[node_list$A]],length(times)))
  set(data_long, , "J", rep(data[[node_list$J]],length(times)))
  set(data_long, , node_list$weights, rep(data[[node_list$weights]],length(times)))
  return(data_long)
}




hazard_to_survival <- function(v, nt, left = F) {
  s <- t(apply(matrix(1-v, ncol = nt),1,cumprod))
  if(left) {
    s <- cbind(rep(1, nrow(s)),s[,-ncol(s)])
  }
  return(as.vector(s))
}

compute_survival_functions <- function(likelihoods, nt) {
  Ft_list <- list()
  surv_N <- hazard_to_survival(likelihoods$N, nt = nt, left = T)
  Ft_p1 <- matrix(surv_N, ncol = nt)*matrix(likelihoods$N, ncol = nt)
  for(i in 1:ncol(likelihoods$J)) {
    Ft <- Ft_p1*matrix(likelihoods$J[[i]], ncol = nt)
    Ft <- t(apply(Ft,1, cumsum))
    Ft_list[[i]] <- as.vector(Ft)

  }
  Ft <- as.matrix(do.call(cbind, Ft_list))
  rm(Ft_list)
  return(list(St = surv_N, Ft = Ft))
}

