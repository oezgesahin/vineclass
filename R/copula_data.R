#' Get copula data and densities
#
#' @description Converts data into u-scale using non-parametric kernel density
#' estimation (unbounded for continuous variables) based on the learning set.
#' Computes marginal densities for learning, and optionally validation and test sets.
#'
#' @param learn A matrix containing learning data on the x-scale. The last column
#' must contain class labels (1 to total_cl).
#' @param test Optional matrix containing test data on the x-scale. The last column must
#' contain class labels (1 to total_cl). If NULL, test data is not processed.
#' @param valid Optional validation matrix on the x-scale. The last column must
#' contain class labels (1 to total_cl). If NULL, validation data is not processed.
#' @param total_cl A positive integer specifying the total number of classes in
#' the data.
#' @param var_types A character vector specifying variable types for each feature
#' ("c" for continuous, "d" for discrete).
#' @return A list containing:
#' \itemize{
#' \item \code{train_u}: Learning set on the u-scale (list of matrices by class).
#' \item \code{dmar_train}: Marginal densities of the learning set (list of matrices by class).
#' \item \code{test_u}: Test set on the u-scale (optional; list of matrices by class).
#' \item \code{dmar_test}: Marginal densities of the test set (optional; list of matrices by class).
#' \item \code{valid_u}: Validation set on the u-scale (optional; list of matrices by class).
#' \item \code{dmar_valid}: Marginal densities of the validation set (optional; list of matrices by class).
#' \item \code{margins}: List of fitted kernel density models for each class and variable.
#' }
#'
#' @examples
#' # Simulate data from two classes
#' set.seed(1)
#' df1 <- matrix(rnorm(900), 300, 3) # Class 1 features
#' df2 <- matrix(rnorm(900, 2, 2), 300, 3) # Class 2 features
#'
#' # Assign class labels
#' df1 <- cbind(df1, class = 1)
#' df2 <- cbind(df2, class = 2)
#'
#' # Combine into a single dataset
#' data <- rbind(df1, df2)
#'
#' # Split data into training (60%), validation (20%), and test (20%) sets
#' train_indices <- 1:360
#' valid_indices <- 361:480
#' test_indices <- 481:600
#'
#' df_train <- data[train_indices, ]
#' df_valid <- data[valid_indices, ]
#' df_test <- data[test_indices, ]
#'
#' # Apply the function with all sets
#' copuladata <- get_u_data_dens(
#'   learn = df_train,
#'   test = df_test,
#'   valid = df_valid,
#'   total_cl = 2,
#'   var_types = c("c", "c", "c")
#' )
#'
#' # Access copula data for training set
#' copuladata_train <- copuladata$train_u
#' copuladata_train_cl1var1 <- copuladata_train[[1]][, 1] # First class, first variable
#'
#' # Access copula data for test set
#' copuladata_test <- copuladata$test_u
#' copuladata_test_cl1var1 <- copuladata_test[[1]][, 1] # First class, first variable
#'
#' # Access copula data for validation set
#' copuladata_valid <- copuladata$valid_u
#' copuladata_valid_cl1var1 <- copuladata_valid[[1]][, 1] # First class, first variable
#'
#' # Access marginal densities of training set
#' univ_kde_train <- copuladata$dmar_train
#' univ_kde_train_cl1var1 <- univ_kde_train[[1]][, 1] # First class, first variable
#'
#' # Access fitted kde1d objects
#' fitted_kdes <- copuladata$margins
#' fitted_kdes_cl1var1 <- fitted_kdes[[1]][[1]] # First class, first variable
#'
#' @import kde1d
#' @importFrom stats ecdf
#' @export

get_u_data_dens <- function(learn, test = NULL, valid = NULL, total_cl, var_types) {
  stopifnot(is.character(var_types), length(var_types) == ncol(learn) - 1)

  # Setup parameters
  total_features <- ncol(learn) - 1
  ncol_udata <- sum(var_types == "c") + 2 * sum(var_types == "d")
  truncate_u <- function(u) pmin(pmax(u, 1e-12), 1-1e-12)

  # Initialize outputs
  initialize_list <- function(rows, cols) replicate(total_cl, matrix(NA, rows, cols), simplify = FALSE)
  dmar_learn <- initialize_list(nrow(learn), total_features)
  u_learn_df <- initialize_list(nrow(learn), ncol_udata + 1)
  fit_kde1d <- vector("list", total_cl)

  if (!is.null(test)) {
    dmar_test <- initialize_list(nrow(test), total_features)
    u_test_df <- initialize_list(nrow(test), ncol_udata + 1)
  }

  if (!is.null(valid)) {
    dmar_valid <- initialize_list(nrow(valid), total_features)
    u_valid_df <- initialize_list(nrow(valid), ncol_udata + 1)
  }

  # Process each class and variable
  for (k in seq_len(total_cl)) {
    count_disc <- 1
    fit_kde1d[[k]] <- vector("list", total_features)
    for (j in seq_len(total_features)) {
      class_indices <- learn[, total_features + 1] == k
      var_learn <- learn[class_indices, j]

      if (var_types[j] == "d") {
        # Discrete variable handling
        fit_kde1d[[k]][[j]] <- ecdf(var_learn)
        u_var <- truncate_u(fit_kde1d[[k]][[j]](learn[, j]))
        uminus <- truncate_u(fit_kde1d[[k]][[j]](learn[, j] - 1))
        pmf <- u_var - uminus

        u_learn_df[[k]][, j] <- u_var
        u_learn_df[[k]][, total_features + count_disc] <- uminus
        dmar_learn[[k]][, j] <- pmf

        apply_to_set <- function(dataset, ecdf_obj) {
          u_val <- truncate_u(ecdf_obj(dataset[, j]))
          uminus_val <- truncate_u(ecdf_obj(dataset[, j] - 1))
          list(u = u_val, pmf = u_val - uminus_val, uminus = uminus_val)
        }

        if (!is.null(valid)) {
          valid_data <- apply_to_set(valid, fit_kde1d[[k]][[j]])
          u_valid_df[[k]][, j] <- valid_data$u
          u_valid_df[[k]][, total_features  + count_disc] <- valid_data$uminus
          dmar_valid[[k]][, j] <- valid_data$pmf
        }

        if (!is.null(test)) {
          test_data <- apply_to_set(test, fit_kde1d[[k]][[j]])
          u_test_df[[k]][, j] <- test_data$u
          u_test_df[[k]][, total_features  + count_disc] <- test_data$uminus
          dmar_test[[k]][, j] <- test_data$pmf
        }
        count_disc <- count_disc + 1

      } else {
        # Continuous variable handling
        fit_kde1d[[k]][[j]] <- kde1d(var_learn)
        apply_continuous <- function(dataset) {
          u_vals <- truncate_u(pkde1d(dataset[, j], fit_kde1d[[k]][[j]]))
          densities <- dkde1d(dataset[, j], fit_kde1d[[k]][[j]])
          na_indices <- is.na(densities)
          max_indx <- which.max(var_learn)
          min_indx <- which.min(var_learn)
          densities[na_indices] <- ifelse(dataset[na_indices, j] > max(var_learn), densities[max_indx], densities[min_indx])
          list(u = u_vals, density = densities)
        }

        train_data <- apply_continuous(learn)
        u_learn_df[[k]][, j] <- train_data$u
        dmar_learn[[k]][, j] <- train_data$density

        if (!is.null(valid)) {
          valid_data <- apply_continuous(valid)
          u_valid_df[[k]][, j] <- valid_data$u
          dmar_valid[[k]][, j] <- valid_data$density
        }

        if (!is.null(test)) {
          test_data <- apply_continuous(test)
          u_test_df[[k]][, j] <- test_data$u
          dmar_test[[k]][, j] <- test_data$density
        }
      }
    }

    # Assign class labels
    u_learn_df[[k]][, ncol_udata + 1] <- learn[, total_features + 1]
    if (!is.null(test)) u_test_df[[k]][, ncol_udata + 1] <- test[, total_features + 1]
    if (!is.null(valid)) u_valid_df[[k]][, ncol_udata + 1] <- valid[, total_features + 1]
  }

  # Return structured results
  result <- list(
    train_u = u_learn_df,
    dmar_train = dmar_learn,
    margins = fit_kde1d
  )
  if (!is.null(valid)) {
    result$valid_u <- u_valid_df
    result$dmar_valid <- dmar_valid
  }
  if (!is.null(test)) {
    result$test_u <- u_test_df
    result$dmar_test <- dmar_test
  }

  return(result)
}
