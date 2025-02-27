#' Evaluate vine copula-based classifiers on test data
#'
#' @description This function calculates posterior probabilities on an external test
#' dataset for a vine copula-based classifier.
#' Given the test data and the vine models (from the training phase),
#' it transforms the test data into u-scale, calculates marginal
#' densities, applies the fitted vine models, and calculates posterior probabilities.
#' @param test_data A matrix or data frame containing test data on the original scale. The data
#'   should have the same structure (number of columns) as the learning data used for model fitting.
#'   The last column must contain class labels.
#' @param model_output A list produced by the vine copula fitting function (e.g., \code{vineclass}),
#'   which contains the fitted vine models, chosen variables, and fitted kernel density
#'   estimates.
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
  total_features <- ncol(test_data)
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


  chosen_cand <- model_output$chosen_vars
  vars <- chosen_cand
  nvars <- length(var_types)

  # Check if the elements at specified indices are "d"
  is_d <- var_types[chosen_cand] == "d"
  # Extract the indices where "d" is found
  d_indices <- chosen_cand[is_d]
  d_positions <- which(var_types == "d")
  indices <- match(d_indices, d_positions)
  if(length(is_d) >0){
    chosen_cand <- c(vars, nvars+indices)
  }
  else{
    chosen_cand <- vars
  }
  vine_mdl <- model_output$vine

  rvine_densities <- sapply(1:total_cl, function(j) {
    result <-  dvinecop(u_test_df[[j]], vine_mdl[[j]])
    result[is.na(result)] <- 1e-5
    result[which(result==0)] <- 1e-5
    return(result)
  })

  lik_points_test <- sapply(1:total_cl, function(j) {
    marginal_densities <- apply(dmar_test[[j]], 1, prod)
    marginal_densities[is.na(marginal_densities)] <- 1e-5
    marginal_densities[which(marginal_densities==0)] <- 1e-5
    result <- marginal_densities * rvine_densities[, j]
    return(result)
  })


  lik_per_obs_test <- apply(lik_points_test, 1, sum)
  post_prob_test <- lik_points_test/rep(lik_per_obs_test, total_cl)
  post_prob_test
}
