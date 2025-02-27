#' Categorical misclassification rates
#' @description This function computes the misclassification rates and average lengths
#' of prediction intervals for each category based on the provided prediction interval matrix.
#'
#' @param predIntervalMatrix A contingency table (matrix or data frame) where rows
#' represent true classes and columns represent predicted intervals.
#' @param iprint Logical; if \code{TRUE}, prints details for each class during computation.
#'
#' @return A matrix with misclassification rates and average prediction interval lengths for
#' each category, along with the overall average for each metric.
#'
#' @examples
#' set.seed(1)
#' # Example data
#' predIntervalMatrix <- matrix(sample(0:5, 12, replace = TRUE), nrow = 3, ncol = 4)
#' rownames(predIntervalMatrix) <- c("A", "B", "C")
#' colnames(predIntervalMatrix) <- c("A", "A,B", "B,C", "C")
#'
#' # Compute misclassification rates and average interval lengths
#' result <- categMisclassRates(predIntervalMatrix, iprint = FALSE)
#' print(result)
#'
#' @import stringr
#' @export

categMisclassRates <- function(predIntervalMatrix, iprint = FALSE) {
  # Validate inputs
  stopifnot(is.matrix(predIntervalMatrix) || is.data.frame(predIntervalMatrix))
  cols <- colnames(predIntervalMatrix)
  rows <- rownames(predIntervalMatrix)

  misclass <- numeric(length(rows))
  avglen <- numeric(length(rows))

  # Compute metrics for each true class
  for (i in seq_along(rows)) {
    row <- predIntervalMatrix[i, ]
    total <- sum(row)
    class <- rows[i]
    correct <- 0
    total_length <- 0

    for (j in seq_along(cols)) {
      col <- cols[j]
      count <- row[col]
      total_length <- total_length + stringr::str_length(col) * count
      if (grepl(class, col, fixed = TRUE)) {
        correct <- correct + count
      }
    }

    avglen[i] <- total_length / total
    misclass[i] <- 1 - correct / total

    if (iprint) {
      cat("Class:", class, "\n")
      cat("Number of misclassifications :", total - correct, "\n")
      cat("Misclassification rate:", round(misclass[i], 3), "\n")
      cat("Average interval length:", round(avglen[i], 2), "\n")
    }
  }

  # Compile results
  result <- rbind(misclass, avglen)
  colnames(result) <- rows
  result <- cbind(result, average = rowMeans(result))
  rownames(result) <- c("Misclassification Rate", "Average Interval Length")

  result
}
