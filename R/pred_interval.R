#' Prediction intervals for classification
#' @description This function computes prediction intervals for categorical classes based on
#' a posterior probability matrix, given a specified prediction interval level (\code{alpha}).
#'
#' @param ProbMatrix A data frame or matrix of estimated posterior probabilities with
#' dimensions \code{n x k}, where \code{k} is the number of categories and \code{n} is
#' the number of observations.
#' @param labels A vector of class names of length \code{k}.
#' @param true_class A vector of true class labels of length \code{n}.
#' @param alpha A numeric value between 0 and 1 specifying the prediction interval level.
#'
#' @return A contingency table showing the true class labels against the predicted intervals.
#'
#' @examples
#' set.seed(1)
#' # Simulate a posterior probability matrix for 3 classes and 10 observations
#' ProbMatrix <- matrix(runif(30, 0, 1), nrow = 10, ncol = 3)
#' ProbMatrix <- ProbMatrix / rowSums(ProbMatrix) # Normalize to probabilities
#' labels <- c("A", "B", "C")
#' true_class <- sample(labels, 10, replace = TRUE)
#' alpha <- 0.8
#'
#' # Compute prediction intervals at alpha=0.8
#' pred_table <- pred_interval(ProbMatrix, labels, true_class, alpha)
#' print(pred_table)
#'
#' @export
pred_interval <- function(ProbMatrix, labels, true_class, alpha) {
  # Validate inputs
  stopifnot(is.matrix(ProbMatrix) || is.data.frame(ProbMatrix))
  stopifnot(length(labels) == ncol(ProbMatrix))
  stopifnot(length(true_class) == nrow(ProbMatrix))
  stopifnot(is.numeric(alpha) && alpha > 0 && alpha <= 1)

  # Initialize prediction intervals
  ncases <- nrow(ProbMatrix)
  predInt <- character(ncases)

  # Compute prediction intervals for each observation
  for (i in seq_len(ncases)) {
    p <- ProbMatrix[i, ]
    ip <- order(p, decreasing = TRUE)
    pOrdered <- p[ip] # Probabilities in decreasing order
    labelsOrdered <- labels[ip] # Corresponding labels
    G <- cumsum(pOrdered) # Cumulative sum of probabilities
    k1 <- min(which(G >= alpha)) # Find smallest k such that sum >= alpha
    pred1 <- labelsOrdered[1:k1] # Select top-k labels
    predInt[i] <- paste(pred1, collapse = ",") # Combine into a single string
  }

  # Create a contingency table
  table(true_class, predInt)
}
