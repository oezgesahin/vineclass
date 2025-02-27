#' Data splitting
#' @description Splits the input dataset into learning (training), validation (optional),
#' and test sets based on the specified proportions.
#'
#' @param data A data frame or matrix to be split.
#' @param learn_rat A numeric value specifying the proportion of data to allocate to the learning (training) set.
#' @param valid_rat Optional numeric value specifying the proportion of data to allocate to the validation set.
#' If not provided, no validation set will be created.
#'
#' @return A list containing:
#' \itemize{
#' \item \code{learning}: The learning (training) set.
#' \item \code{valid}: The validation set (if \code{valid_rat} is provided).
#' \item \code{test}: The test set.
#' }
#'
#' @examples
#' set.seed(1)
#' data <- matrix(rnorm(900), ncol = 3)
#'
#' # Split into training (60%), validation (20%), and test (20%)
#' split <- split_data(data, learn_rat = 0.6, valid_rat = 0.2)
#' train <- split$learning
#' valid <- split$valid
#' test <- split$test
#'
#' # Split into training (70%) and test (30%) without validation
#' split <- split_data(data, learn_rat = 0.7)
#' train <- split$train
#' test <- split$test
#' @export

split_data <- function(data, learn_rat, valid_rat = NA) {
  # Validate inputs
  stopifnot(is.data.frame(data) || is.matrix(data))
  stopifnot(is.numeric(learn_rat) && learn_rat > 0 && learn_rat < 1)
  if (!is.na(valid_rat)) stopifnot(is.numeric(valid_rat) && valid_rat > 0 && valid_rat < 1)
  stopifnot(is.na(valid_rat) || (learn_rat + valid_rat < 1))

  total_obs <- nrow(data)
  test_rat <- if (is.na(valid_rat)) 1 - learn_rat else 1 - learn_rat - valid_rat

  # Split data
  learn_obs <- sample(total_obs, round(total_obs * learn_rat))
  remaining_obs <- setdiff(1:total_obs, learn_obs)

  colnames(data)[ncol(data)] <- "class"

  if (!is.na(valid_rat)) {
    val_obs <- sample(remaining_obs, round(total_obs * valid_rat))
    test_obs <- setdiff(remaining_obs, val_obs)
    list(
      learning = data[learn_obs, ],
      valid = data[val_obs, ],
      test = data[test_obs, ]
    )
  } else {
    list(
      train = data[learn_obs, ],
      test = data[remaining_obs, ]
    )
  }
}
