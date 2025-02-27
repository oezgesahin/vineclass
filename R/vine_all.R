#' Fit vine models to all classes
#' @description
#' Fits vine copula models to the given data for each class, based on selected variables.
#' The function uses the `rvinecopulib` package for vine copula estimation.
#'
#' @param learn A list where each item corresponds to the copula data of a class on the u-scale.
#'   The last column contains class labels (1 to \code{total_cl}).
#' @param vars An integer vector specifying the indices of selected variables for modeling.
#' @param total_cl A positive integer indicating the total number of classes in the data.
#' @param dmar_learn A list of marginal densities for the data, where each item corresponds to
#'   the marginal densities of a class.
#' @param copfam A character vector specifying the pair copula families to be used.
#' @param nCore A positive integer specifying the number of cores to use for vine copula fitting.
#' @param trunc_lvl A numeric value specifying the truncation level for the vine copula models. Inf means no truncation
#' @param var_types A character vector of variable types (\code{"c"} for continuous,
#'   \code{"d"} for discrete).
#'
#' @return A list with the following components:
#' \describe{
#'   \item{\code{vine}}{A list of fitted vine copula models for each class.}
#'   \item{\code{chosen_cand}}{A vector of selected variables used in the vine copula models.}
#'   \item{\code{post_prob_learn}}{A data frame of posterior probabilities for each observation in the learning set.}
#' }
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' df1 <- matrix(rnorm(900), 300, 3)
#' df2 <- matrix(rnorm(900, 2, 2), 300, 3)
#' df1 <- cbind(df1, class = 1)
#' df2 <- cbind(df2, class = 2)
#' df <- rbind(df1, df2)
#' df_splitted <- split_data(df, 0.8)
#' df_train <- df_splitted$train
#' df_test <- df_splitted$test
#' copuladata <- get_u_data_dens(learn=df_train, test=df_test, total_cl=2, var_types= c("c", "c", "c"))
#' vinefit <- vine_all(copuladata$train_u, 1:3, 2, copuladata$dmar_train,"parametric", 1, 2, c("c", "c", "c"))
#' }
#'
#' @import rvinecopulib
#' @importFrom parallel mclapply
#' @noRd
vine_all <- function(learn, vars, total_cl, dmar_learn, copfam, nCore, trunc_lvl, var_types) {
  cols <- ncol(learn[[1]])
  vine_mdl <- list()
  nvars <- length(var_types)
  disc_vine <- list()

  # Identify discrete variables and their positions
  is_d <- var_types[vars] == "d"
  d_indices <- vars[is_d]
  d_positions <- which(var_types == "d")
  indices <- match(d_indices, d_positions)

  # Combine selected variables with discrete variable indices
  if (length(is_d) > 0) {
    chosen_cand <- c(vars, nvars + indices)
  } else {
    chosen_cand <- vars
  }
  # Fit vine models with discrete handling
  if (any(var_types == "d") && length(nvars) > 2) {
    for (k in 1:total_cl) {
      input_df <- learn[[k]][which(learn[[k]][, cols] == k), vars]
      n1 <- nrow(input_df)
      rr1 <- vWcor(input_df, var_types[vars])
      disc_vine[[k]] <- GaussVineMST(rr1, n1, 1.01)$VineA[1:length(vars), length(vars):1]
    }
    vine_mdl <- parallel::mclapply(1:total_cl, function(k) {
      vinecop(
        learn[[k]][which(learn[[k]][, cols] == k), chosen_cand],
        presel = FALSE,
        structure = disc_vine[[k]],
        family_set = copfam,
        trunc_lvl = trunc_lvl,
        var_types = var_types[vars]
      )
    }, mc.cores = nCore)
  } else {
    vine_mdl <- parallel::mclapply(1:total_cl, function(k) {
      vinecop(
        learn[[k]][which(learn[[k]][, cols] == k), chosen_cand],
        presel = FALSE,
        family_set = copfam,
        trunc_lvl = trunc_lvl,
        var_types = var_types[vars]
      )
    }, mc.cores = nCore)
  }

  # Compute densities and posterior probabilities
  dvine_learn <- sapply(1:total_cl, function(j) {
    result <- dvinecop(learn[[j]][, chosen_cand], vine_mdl[[j]])
    result[is.na(result)] <- 1e-5
    result[result == 0] <- 1e-5
    result
  })

  lik_points_learn <- sapply(1:total_cl, function(j) {
    margin_dens <-apply(dmar_learn[[j]][, vars], 1, prod)
    margin_dens[is.na(margin_dens)] <- 1e-5
    margin_dens[margin_dens == 0] <- 1e-5
    result <- margin_dens * dvine_learn[, j]
    result
  })

  lik_per_obs_learn <- apply(lik_points_learn, 1, sum)
  post_prob_learn <- lik_points_learn / rep(lik_per_obs_learn, total_cl)

  list(
    vine = vine_mdl,
    chosen_cand = vars,
    post_prob_learn = post_prob_learn
  )
}
