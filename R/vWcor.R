#' Vander Waerden correlation matrix
#' @description
#' Computes a correlation matrix for mixed continuous and discrete data types:
#' - Continuous-continuous: Vander Waerden correlation of normal scores.
#' - Continuous-discrete: Polyserial correlation using \code{\link[polycor]{polyserial}}.
#' - Discrete-discrete: Polychoric correlation using \code{\link[psych]{polychoric}}.
#'
#' @param data A dataframe or matrix containing the data.
#' @param var_types A character vector of length equal to the number of columns in \code{data},
#'   where each element is "c" for continuous or "d" for discrete.
#'
#' @return A symmetric correlation matrix.
#'
#' @importFrom psych polychoric
#' @importFrom polycor polyserial
#' @importFrom Matrix nearPD
#' @importFrom stats cor
#'
#' @examples
#' set.seed(123)
#' data <- data.frame(
#'   cont1 = rnorm(100),
#'   cont2 = rnorm(100),
#'   disc1 = sample(0:1, 100, replace = TRUE),
#'   disc2 = sample(0:2, 100, replace = TRUE)
#' )
#' var_types <- c("c", "c", "d", "d")
#' vWcor(data, var_types)
#'
#' @export
vWcor <- function(data, var_types) {
  d <- ncol(data)
  stopifnot(d == length(var_types))

  rmat <- matrix(0, d, d)
  is_disc <- ifelse(var_types == "d", TRUE, FALSE)
  at_least_two_c <- sum(var_types == "c") >= 2

  if (at_least_two_c) {
    nscore_continuous <- norm_score(data[, !is_disc])
    rmat[!is_disc, !is_disc] <- cor(nscore_continuous)
  }

  for (i in which(is_disc)) {
    for (j in 1:d) {
      if (i == j) next

      if (!is_disc[j]) {
        # i is discrete, j is continuous
        rmat[i, j] <- rmat[j, i] <- polycor::polyserial(data[, j], data[, i])
      } else {
        # both i and j are discrete
        rmat[i, j] <- rmat[j, i] <- psych::polychoric(table(data[, i], data[, j]))$rho
      }
    }
  }

  # Ensure positive definiteness
  diag(rmat) <- 1
  if (!isPosDef(rmat)) {
    rmat <- as.matrix(Matrix::nearPD(rmat, corr = TRUE)$mat)
  }

  rmat
}
