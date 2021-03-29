
make_bins <- function(x, type = c("equal_range", "equal_mass"), n_bins = NULL) {
  # clean up arguments
  type <- match.arg(type)
  
  # set grid along x
  if (type == "equal_range") {
    bins <- ggplot2::cut_interval(x, n_bins,
                                  right = FALSE,
                                  ordered_result = TRUE, dig.lab = 12
    )
  } else if (type == "equal_mass") {
    bins <- ggplot2::cut_number(x, n_bins,
                                right = FALSE,
                                ordered_result = TRUE, dig.lab = 12
    )
  }
  
  # https://stackoverflow.com/questions/36581075/extract-the-breakpoints-from-cut
  breaks_left <- as.numeric(sub(".(.+),.+", "\\1", levels(bins)))
  breaks_right <- as.numeric(sub(".+,(.+).", "\\1", levels(bins)))
  breaks <- c(breaks_left[1], breaks_right)
  return(breaks)
}
discretize_variable <- function(x, type = c("equal_range", "equal_mass"), n_bins = NULL, breaks = NULL) {
  if (is.null(breaks)) {
    breaks <- make_bins(x, type, n_bins)
  }
  
  # for predict method, only need to assign observations to existing intervals
  # NOTE: findInterval() and cut() might return slightly different results...
  bin_id <- findInterval(x, breaks, all.inside = T)
  x_in_bin <- x - breaks[bin_id]
  list(x_discrete = bin_id, x_in_bin = x_in_bin, breaks = breaks)
}



Lrnr_density_discretize <- R6::R6Class(
  classname = "Lrnr_density_discretize",
  inherit = Lrnr_base, portable = TRUE,
  class = TRUE,
  public = list(
    initialize = function(categorical_learner = NULL, type = "equal_mass", n_bins = 20, predict_type = c("density", "CDF", "RCDF", "probability"), breaks = NULL, ...) {
      predict_type <- match.arg(predict_type)
      if (is.null(categorical_learner)) {
        categorical_learner <- make_learner(Lrnr_glmnet)
      }
      params <- list(
        type = type, n_bins = n_bins, predict_type = predict_type,
        categorical_learner = categorical_learner, breaks = breaks, ...
      )
      super$initialize(params = params, ...)
    }
  ),
  
  private = list(
    .properties = c("density"),
    
    .train = function(task) {
      discretized <- sl3:::discretize_variable(task$Y,
                                               type = self$params$type,
                                               n_bins = self$params$n_bins,
                                               breaks = self$params$breaks
      )
      
      # make discretized task
      new_columns <-
        task$add_columns(data.table(
          discrete_Y =
            factor(discretized$x_discrete)
        ))
      discrete_task <- task$next_in_chain(
        outcome = "discrete_Y",
        column_names = new_columns
      )
      # fit categorical learner to discretized task
      categorical_fit <- self$params$categorical_learner$train(discrete_task)
      
      fit_object <- list(
        categorical_fit = categorical_fit,
        breaks = discretized$breaks,
        levels = levels(factor(discretized$x_discrete))
      )
      return(fit_object)
    },
    
    .predict = function(task) {
      # make discretized task
      discretized <- sl3:::discretize_variable(task$Y,
                                               breaks = self$fit_object$breaks
      )
      data <- task$data
      data$discrete_Y =
        factor(discretized$x_discrete)
      nodes <- task$nodes
      nodes$outcome <- "discrete_Y"
      discrete_task <- sl3_Task$new(data = data, nodes = nodes, outcome_type  = "categorical", outcome_levels = self$fit_object$levels, row_index <- task$row_index, folds = task$folds)
      # new_columns <-
      #   task$add_columns(data.table(
      #     discrete_Y =
      #       factor(discretized$x_discrete)
      #   ))
      # discrete_task <- task$next_in_chain(
      #   outcome = "discrete_Y",
      #   column_names = new_columns,
      #   outcome_levels = self$fit_object$levels
      # )
      
      # predict categorical learner on discretized task
      raw_preds <- self$fit_object$categorical_fit$predict(discrete_task)
      predmat <- unpack_predictions(raw_preds)
      if(self$params$predict_type ==  "probability") {
        obs_pred <- predmat[cbind(seq_len(task$nrow), discretized$x_discrete)]
      } else if(self$params$predict_type != "density") {
        predmat <- t(apply(predmat, 1, function(v) {
          c(0, cumsum(v))
        }))
         bin_lengths <- diff(self$fit_object$breaks)
        obs_pred_left <- predmat[cbind(seq_len(task$nrow), discretized$x_discrete)]
        obs_pred_right <- predmat[cbind(seq_len(task$nrow), discretized$x_discrete+1)]
        x_diff <- discretized$x_in_bin
        obs_pred <- obs_pred_left + x_diff*(obs_pred_right - obs_pred_left)/bin_lengths[discretized$x_discrete]
        obs_pred <- pmin(obs_pred,1)
        obs_pred <- pmax(obs_pred,0)
        
        if(self$params$predict_type == "RCDF") {
          obs_pred <- 1 - obs_pred
        }
      } else {
        bin_lengths <- diff(self$fit_object$breaks)
        scale_mat <- matrix(rep(1 / bin_lengths, each = task$nrow),
                            nrow = task$nrow
        )
        predmat <- predmat * scale_mat
        obs_pred <- predmat[cbind(seq_len(task$nrow), discretized$x_discrete)]
        
      }
      # subset predictions to only those bins relevant
      
      
      return(obs_pred)
    },
    .required_packages = c()
  )
)