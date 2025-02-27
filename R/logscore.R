#' Negative log-likelihood score for classification
#' @description Calculate the negative log-likelihood score for classification.
#' @param predmatrix A numeric matrix where each row contains predicted probabilities for each class.
#' @param labels A numeric vector of true class indices corresponding to each row of predmatrix.
#' @return The negative log-likelihood score (numeric).
#' @export
logscore <- function(predmatrix, labels) {
  # Check inputs
  if (length(labels) != nrow(predmatrix)) {
    stop("Number of labels must match the number of rows in predmatrix.")
  }
  eps <- 1e-5

  # Extract probabilities corresponding to true labels
  probs <- predmatrix[cbind(1:nrow(predmatrix), labels)]

  # Replace any probabilities below the threshold with the threshold
  probs <- pmax(probs, eps)

  # Compute the negative log-likelihood
  score <- -sum(log(probs))

  return(score)
}
