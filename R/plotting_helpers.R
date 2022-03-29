
#' @export
plot_threshold_response <- function(output, simultaneous_CI = T, monotone = F) {
  estimates <- output
  if (simultaneous_CI) {
    lower_index <- 6
    upper_index <- 7
    title <- "Threshold-response function with simultaneous 95% confidence bands"
  } else {
    lower_index <- 4
    upper_index <- 5
    title <- "Threshold-response function with point-wise 95% confidence intervals"
  }
  lower <- bound(estimates[, lower_index], c(0, 2 * max(estimates[, 2])))
  upper <- bound(estimates[, upper_index], c(0, 2 * max(estimates[, 2])))
  if (T | !(monotone & simultaneous_CI)) {
    no_event <- attr(output, "no_event")
    upper[no_event] <- 0
  }
  subtitle <- NULL

  if (monotone) {
    subtitle <- paste0("Assumes monotonicity of true function")
    mon <- isoreg(estimates[, 1], -estimates[, 2])
    estimates[, 2] <- -mon$yf
    if (simultaneous_CI) {
      for (i in 2:length(upper)) {
        if(is.na(upper[i - 1]) || is.na(lower[i - 1])) {
          next
        }
        if (upper[i] > upper[i - 1]) {
          upper[i] <- upper[i - 1]
        }
        if (lower[i] > lower[i - 1]) {
          lower[i] <- lower[i - 1]
        }
      }
    }
  }

  # vec <- c()
  # remove <- c()
  # for(k in 1:(nrow(est))) {
  #   if(k==1) {
  #     deriv <-  (est[k+1,1] - est[k,1])/(est[k+1,3] - est[k,3])
  #   }
  #   else if(k==nrow(est)) {
  #     deriv <-  (est[k-1,1] - est[k,1])/(est[k-1,3] - est[k,3])
  #   } else {
  #     deriv <- (est[k+1,1] - est[k,1])/(est[k+1,3] - est[k,3])/2 + (est[k-1,1] - est[k,1])/(est[k-1,3] - est[k,3])/2
  #   }
  #   index <- which.min(abs(esttmle[,1] -  est[k,3]))
  #   print(index)
  #   radius <- 1.96*sd(IC[,index]/deriv)/sqrt(nrow(IC))
  #
  #   vec <- c( vec, est[k,3]- radius , est[k,3], radius + est[k,3])
  # }
  # mat <- cbind(est[,1], matrix(vec, ncol = 3, byrow = T))
  # colnames(mat) <- c("risk", "CI_lower", "thresh_est", "CI_upper")
  # mat


  if (!is.null(estimates)) {
    library(ggplot2)
    plot_data <- data.frame(
      cutoffs = round(estimates[, 1], 2), est = estimates[, 2],
      lower = lower,
      upper = upper
    )
    g1 <- ggplot2::ggplot(data = plot_data, aes_string("cutoffs", "est")) +
      scale_x_continuous(breaks = plot_data$cutoffs) +
      #geom_point(aes_string(x = "cutoffs", y = "est"), legend = F, colour = alpha("red")) +
      #geom_point(aes_string(x = "cutoffs", y = "lower"), legend = F, colour = alpha("red", 0.2)) +
      #geom_point(aes_string(x = "cutoffs", y = "upper"), legend = F, colour = alpha("red", 0.2)) +
      geom_line() +
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, color = NA)

    g1 <- g1 +
      theme(
        text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, hjust = 1)
      ) + ylab("Cumulative incidence") + xlab("Threshold")
    g1 <- g1 + ggtitle(title) + labs(subtitle = subtitle) + theme(plot.title = element_text(size = 15), plot.subtitle = element_text(size = 13))
  } else {
    g1 <- NULL
  }



  return(g1)
}
