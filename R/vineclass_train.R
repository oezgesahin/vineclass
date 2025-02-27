#' Vine copula-based classifiers
#' @description
#' This is the main function for vine copula-based classification. It splits the
#' data into learning and test sets (if applicable), fits vine copula models for each class,
#' and returns posterior probabilities for both the learning and test sets.
#'
#' @param data A data frame or matrix. The last column should contain class labels,
#'   which will be re-labeled internally to 1, 2, ..., K if not already.
#' @param vars A vector of integer indices specifying which columns (excluding the class label)
#'   are used for modeling. If missing or \code{NULL}, all variables except the last (class label)
#'   are used.
#' @param learn_ratio A numeric value in \eqn{(0,1]} specifying the proportion of data for learning.
#'   If \code{learn_ratio = 1}, all data is used for learning (no test set, default).
#' @param copfam A character vector specifying the pair copula families (e.g., "parametric", "nonparametric").
#' @param nCore A positive integer specifying the number of cores to use for parallel computations.
#'   If missing or \code{NULL}, it is set to the minimum of the total number of classes and
#'   \code{parallel::detectCores() - 1}.
#' @param trunc_lvl A numeric value specifying the truncation level for the vine copula models. Inf means no truncation
#' @param var_types A character vector of length \code{ncol(data) - 1}, indicating "c" or "d" for each variable.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{vine}}{List of fitted vine copula models for each class.}
#'   \item{\code{chosen_vars}}{Vector of selected variables used for vine modeling.}
#'   \item{\code{fitted_kdes}}{List of KDEs (marginal densities) for each class in the learning set.}
#'   \item{\code{post_prob_learn}}{Posterior probabilities for the learning data.}
#'   \item{\code{post_prob_test}}{Posterior probabilities for the test data (if a test set exists).}
#'   \item{\code{test_data}}{Data frame containing the test data used for evaluation (if a test set exists).}
#' }
#'
#' @references
#' Sahin O and Joe H (2024), Vine Copula-Based Classifiers with Applications.
#' Journal of Classification, 1-29.
#'
#' Brechmann E C and Joe H (2015), Truncation of vine copulas using fit indices.
#' Journal of Multivariate Analysis, 138, 19-33.
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' data <- data.frame(matrix(rnorm(1000), ncol = 5))
#' data$class <- sample(1:2, nrow(data), replace = TRUE)
#'
#' # Run vine-based classification using 80% of the data for learning
#' result <- vineclass(
#'   data=data,
#'   learn_ratio = 0.8,
#'   copfam = "parametric",
#'   trunc_lvl = Inf,
#'   var_types = c("c", "c", "c", "c","c")
#' )
#' head(result$post_prob_learn)
#' head(result$post_prob_test)
#' }
#'
#' @export
vineclass <- function(data,
                      vars = NULL,
                      learn_ratio = 1.0,
                      copfam = "parametric",
                      nCore = NULL,
                      trunc_lvl = Inf,
                      var_types)
{
  # Validate 'data'
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("'data' must be a data frame or matrix.")
  }

  # Re-label the class column to be consecutive integers starting at 1
  class_col_index <- ncol(data)
  original_classes <- unique(data[, class_col_index])
  # Sort classes to ensure consistent ordering
  sorted_classes <- sort(original_classes)
  # Convert to factor with levels in ascending order, then to integer
  data[, class_col_index] <- as.integer(
    factor(
      data[, class_col_index],
      levels = sorted_classes,
      labels = seq_along(sorted_classes)
    )
  )

  # If 'vars' is NULL, use all columns except the last
  if (is.null(vars)) {
    vars <- seq_len(ncol(data) - 1)
  }

  # Validate 'vars'
  if (!all(vars %in% seq_len(ncol(data) - 1))) {
    stop("Elements of 'vars' must be valid column indices of 'data', excluding the class label.")
  }

  # Validate 'learn_ratio'
  if (!is.numeric(learn_ratio) || learn_ratio <= 0 || learn_ratio > 1) {
    stop("learn_ratio must be between 0 and 1.")
  }

  # Validate 'var_types'
  if (length(var_types) != (ncol(data) - 1)) {
    stop("'var_types' must have length equal to the number of variables (columns minus the class label).")
  }

  # Determine the total number of classes after re-labeling
  total_classes <- length(unique(data[, class_col_index]))

  # Default 'nCore' if not specified: min(total_classes, detectCores()-1)
  if (is.null(nCore)) {
    nCore <- min(total_classes, max(1, parallel::detectCores() - 1))
  }

  # Split data into learning and test sets
  if (learn_ratio < 1) {
    split <- split_data(data, learn_ratio)
    learn <- split$train
    test <- split$test
  } else {
    learn <- data
    test <- NULL
  }

  # Get copula data (u-scale) and marginal densities from the learning (and optional test) sets
  copula_data <- get_u_data_dens(learn=learn, test=test, total_cl = total_classes, var_types= var_types)

  # Extract necessary components
  learn_u <- copula_data$train_u
  dmar_learn <- copula_data$dmar_train

  test_u <- if (!is.null(test)) copula_data$test_u else NULL
  dmar_test <- if (!is.null(test)) copula_data$dmar_test else NULL
  # Fit vine copulas for each class and compute learning posterior probabilities
  vine_fit <- vine_all(
    learn = learn_u,
    vars = vars,
    total_cl = total_classes,
    dmar_learn = dmar_learn,
    copfam = copfam,
    nCore = nCore,
    trunc_lvl = trunc_lvl,
    var_types = var_types
  )
  # Retrieve posterior probabilities for the learning data
  post_prob_learn <- vine_fit$post_prob_learn

  # Compute posterior probabilities for the test set (if it exists)
  if (!is.null(test_u)) {
    post_prob_test <- vine_out(test_u, dmar_test, total_classes, vine_fit$chosen_cand, var_types, vine_fit$vine)$post_prob
  } else {
    post_prob_test <- NULL
  }

  # Return final results
  list(
    vine = vine_fit$vine,
    chosen_vars = vine_fit$chosen_cand,
    fitted_kdes = copula_data$margins,
    post_prob_learn = post_prob_learn,
    post_prob_test = post_prob_test,
    train_data = learn,
    test_data = test
  )
}
