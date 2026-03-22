#' Evaluate vine copula-based classifiers on test data
#'
#' @description This function calculates posterior probabilities on an external test
#' dataset for a vine copula-based classifier.
#' Given the test data and the vine models (from the training phase),
#' it transforms the test data into u-scale, calculates marginal
#' densities, applies the fitted vine models, and calculates posterior probabilities
#' using the class priors stored in \code{model_output}.
#' @param test_data A matrix or data frame containing test data on the original scale. The data
#'   should have the same predictor columns as the learning data used for model fitting.
#'   If a class-label column is present, it is ignored.
#' @param model_output A list produced by the vine copula fitting function (e.g., \code{vineclass}),
#'   which contains the fitted vine models, chosen variables, and fitted kernel density
#'   estimates, along with class prior probabilities.
#' @param var_types A character vector indicating the type of each variable (e.g., \code{"c"} for
#'   continuous and \code{"d"} for discrete). Its length should equal the number of variable columns.
#' @return A matrix of posterior probabilities for each test observation across all classes.
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' data <- data.frame(matrix(rnorm(1000), ncol = 5))
#' data$class <- sample(1:2, nrow(data), replace = TRUE)
#'
#' # Run vine-based classification
#' result <- vineclass(
#'   data=data,
#'   learn_ratio = 1,
#'   copfam = "parametric",
#'   trunc_lvl = Inf,
#'   var_types = c("c", "c", "c", "c","c")
#' )
#'
#' test_data <- matrix(rnorm(100),ncol = 5)
#' var_types <- c("c", "c", "c", "c", "c")
#' test_post_prob <- test_vineclass(test_data, result, var_types)
#' head(test_post_prob)
#' }
#'
#' @export

test_vineclass <- function(test_data, model_output, var_types) {
  truncate_u <- function(u) pmin(pmax(u, 1e-12), 1-1e-12)
  total_features <- length(var_types)
  total_obs_test <- nrow(test_data)
  count_c <- sum(var_types == "c")
  count_d <- sum(var_types == "d")
  ncol_udata <- count_c + 2*count_d
  total_cl <- ncol(model_output$post_prob_learn)
  dmar_test <- lapply(1:total_cl, matrix, data=NA, nrow=total_obs_test, ncol=total_features)
  u_test_df <- lapply(1:total_cl, matrix, data=NA, nrow=total_obs_test, ncol=ncol_udata)
  fit_kde1d <- model_output$fitted_kdes
  for(k in 1:total_cl){
    count_disc <- 1
    for(j in 1:total_features){
      if(var_types[j]=="d"){
        # Using the fitted kernels on the learning set, apply PIT on the test set ----
        u_test <- fit_kde1d[[k]][[j]](test_data[,j])
        u_test <- truncate_u(u_test)
        u_test_df[[k]][,j] <- u_test

        u_test_minus <- fit_kde1d[[k]][[j]](test_data[,j] - 1)
        u_test_minus <- truncate_u(u_test_minus)
        u_test_df[[k]][,total_features+count_disc] <- u_test_minus
        # Using the fitted kernels on the learning set, calculate marginal densities on the test set -----
        dmar_test[[k]][,j]<- u_test - u_test_minus
        count_disc <- count_disc + 1
      }
      else{

        # Using the fitted kernels on the learning set, apply PIT on the test set ----
        u_test_df[[k]][,j] <- pkde1d(test_data[,j], fit_kde1d[[k]][[j]])
        u_test_df[[k]][,j]  <- truncate_u(u_test_df[[k]][,j])

        # Using the fitted kernels on the learning set, calculate marginal densities on the test set -----
        dmar_test[[k]][,j] <- dkde1d(test_data[,j], fit_kde1d[[k]][[j]])
        na_indices <- is.na(dmar_test[[k]][, j])
        if (any(na_indices)) {
          # Check if corresponding test data values are too high
          high_values <- test_data[na_indices, j] > max(test_data[,j])
          #dmar_test[[k]][na_indices, j] <- ifelse(high_values, 0.99999, 1e-5)
        }

      }
    }
  }

  class_priors <- if (!is.null(model_output$class_priors)) {
    resolve_class_priors(total_cl = total_cl, class_priors = model_output$class_priors)
  } else if (!is.null(model_output$train_data)) {
    resolve_class_priors(
      labels = model_output$train_data[, ncol(model_output$train_data)],
      total_cl = total_cl
    )
  } else {
    rep(1 / total_cl, total_cl)
  }

  vine_out(
    test = u_test_df,
    dmar_test = dmar_test,
    total_cl = total_cl,
    chosen_cand = model_output$chosen_vars,
    var_types = var_types,
    vine_mdl = model_output$vine,
    class_priors = class_priors
  )$post_prob
}
