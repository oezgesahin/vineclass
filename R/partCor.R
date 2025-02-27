#' Partial correlation estimation
#' @description
#' Computes the partial correlation between two variables, \code{j} and \code{k},
#' given a set of conditioning variables using the correlation matrix.
#'
#' @param S A correlation matrix.
#' @param given A vector of indices for the conditioning variables.
#' @param j An index for the first variable.
#' @param k An index for the second variable.
#'
#' @return A numeric value representing the partial correlation between \code{j} and \code{k}
#'   given the variables in \code{given}.
#'
#' @examples
#' # Example correlation matrix
#' S <- matrix(c(1, 0.5, 0.3, 0.5, 1, 0.2, 0.3, 0.2, 1), ncol = 3)
#' # Partial correlation of variables 1 and 3 given variable 2
#' partCor(S, given = 2, j = 1, k = 3)
#'
#' @export
partCor <- function(S, given, j, k) {
  S11 <- S[given, given]
  jk <- c(j, k)
  S12 <- S[given, jk]
  S21 <- S[jk, given]
  S22 <- S[jk, jk]

  if (length(given) > 1) {
    tem <- solve(S11, S12)
    Om212 <- S21 %*% tem
  } else {
    tem <- S12 / S11
    Om212 <- outer(S21, tem)
  }

  om11 <- 1 - Om212[1, 1]
  om22 <- 1 - Om212[2, 2]
  om12 <- S[j, k] - Om212[1, 2]

  om12 / sqrt(om11 * om22)
}
