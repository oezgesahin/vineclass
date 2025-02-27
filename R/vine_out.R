#' Test vine-classifier
#' @description This function evaluates the performance of vine copula-based classifiers
#' on a test set by computing class labels and posterior probabilities for each observation internally.
#'
#' @param test A list of test sets on u-scale where each list item corresponds
#' to the copula data of a class. The last column must contain class labels.
#' @param dmar_test A list of marginal densities of the test set, where each
#' list item corresponds to the marginal densities of a class.
#' @param total_cl A positive integer specifying the total number of classes in the data.
#' @param chosen_cand A vector of integers indicating the indices of selected variables.
#' @param var_types A character vector of length \code{ncol(data) - 1}, indicating "c" or "d" for each variable.
#' @param vine_mdl A list of fitted vine copula models for each class.
#'
#' @return A list containing:
#' \itemize{
#' \item \code{vine_lbl}: A vector containing the predicted class labels for each observation.
#' \item \code{post_prob}: A data frame of posterior probabilities with dimensions
#' \code{n_observations x total_cl}.
#' }
#'
#'
#' @noRd
vine_out <- function(test, dmar_test, total_cl, chosen_cand, var_types, vine_mdl) {
  # Validate inputs
  stopifnot(is.list(test), is.list(dmar_test), is.list(vine_mdl))
  stopifnot(length(test) == total_cl, length(dmar_test) == total_cl, length(vine_mdl) == total_cl)
  stopifnot(is.numeric(chosen_cand) && all(chosen_cand > 0))

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


  # Compute R-vine densities for each class
  rvine_densities <- sapply(seq_len(total_cl), function(j) {
    result <- dvinecop(test[[j]][, chosen_cand], vine_mdl[[j]])
    result[is.na(result)] <- 1e-5
    result[result == 0] <- 1e-5
    result
  })


  # Compute likelihood points for each class
  lik_points_test <- sapply(seq_len(total_cl), function(j) {
    margin_dens <-apply(dmar_test[[j]][, vars], 1, prod)
    margin_dens[is.na(margin_dens)] <- 1e-5
    margin_dens[margin_dens == 0] <- 1e-5
    result <- margin_dens * rvine_densities[, j]
    result
  })

  # Compute posterior probabilities
  lik_per_obs_test <- rowSums(lik_points_test)
  post_prob <- lik_points_test / matrix(rep(lik_per_obs_test, total_cl), ncol = total_cl, byrow = FALSE)

  # Determine predicted class labels
  vine_lbl <- apply(post_prob, 1, function(x) which.max(x))

  # Return results
  list(vine_lbl = vine_lbl, post_prob = as.data.frame(post_prob))
}
