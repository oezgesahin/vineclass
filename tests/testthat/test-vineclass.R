test_that("vineclass handles invalid data inputs", {
  # Non-data.frame or non-matrix 'data'
  non_df_data <- list(x = 1:10, y = 1:10)
  expect_error(
    vineclass(non_df_data, var_types = c("c", "c")),
    "'data' must be a data frame or matrix."
  )
})

test_that("vineclass re-labels class column to 1, 2, ..., K", {
  # Classes are re-labeled
  set.seed(1)
  data <- data.frame(matrix(rnorm(900), ncol = 3))
  data$class <- sample(c(5, 10), nrow(data), replace = TRUE) # Non-sequential labels

  # Provide minimal var_types
  # - ncol(data) - 1 = 3 - 1 = 2
  # - so we need length 2 for var_types
  result <- vineclass(
    data        = data,
    var_types   = c("c", "c","c"),
    learn_ratio = 0.8
  )

  # Check if the returned train_data and test_data have sequential classes
  train_classes <- unique(result$train_data[, ncol(result$train_data)])
  test_classes  <- if (!is.null(result$test_data)) {
    unique(result$test_data[, ncol(result$test_data)])
  } else NULL

  # Expect classes are in {1, 2} because we used only 2 distinct original classes (5, 10).
  expect_true(all(train_classes %in% c(1, 2)))
  if (!is.null(test_classes)) {
    expect_true(all(test_classes %in% c(1, 2)))
  }
})

test_that("vineclass uses all columns except last if vars is NULL", {
  # If vars is NULL, it uses all columns except last
  set.seed(1)
  data <- data.frame(matrix(rnorm(900), ncol = 5))
  data$class <- sample(1:2, nrow(data), replace = TRUE)

  # var_types: length = ncol(data) - 1 = 5
  var_types <- rep("c", 5) # 5 numeric columns, last is class
  result <- vineclass(
    data        = data,
    learn_ratio = 1,  # no test set
    var_types   = var_types
  )
  # chosen_vars should be 1:5-1 = 1:4
  expect_equal(result$chosen_vars, c(1, 2, 3, 4, 5))
})

test_that("vineclass stops if vars indices are invalid", {
  # An invalid 'vars' triggers an error
  set.seed(1)
  data <- data.frame(matrix(rnorm(900), ncol = 5))
  data$class <- sample(1:2, nrow(data), replace = TRUE)
  var_types <- rep("c", 5)

  # Provide invalid 'vars' (6 is out of range)
  expect_error(
    vineclass(data, vars = c(2, 6), var_types = var_types),
    "Elements of 'vars' must be valid column indices of 'data', excluding the class label."
  )
})

test_that("vineclass stops if learn_ratio is out of range", {
  # provide invalid learn_ratio
  set.seed(1)
  data <- data.frame(matrix(rnorm(900), ncol = 5))
  data$class <- sample(1:2, nrow(data), replace = TRUE)
  var_types <- rep("c", 5)

  expect_error(
    vineclass(data, learn_ratio = 1.1, var_types = var_types),
    "learn_ratio must be between 0 and 1."
  )
})

test_that("vineclass stops if var_types length is incorrect", {
  # var_types that's too short/long
  set.seed(1)
  data <- data.frame(matrix(rnorm(900), ncol = 5))
  data$class <- sample(1:2, nrow(data), replace = TRUE)

  # We need length = ncol(data) - 1 = 5 - 1 = 4
  # Provide length 3
  var_types <- c("c", "c", "c")

  expect_error(
    vineclass(data, var_types = var_types),
    "'var_types' must have length equal to the number of variables"
  )
})

test_that("vineclass splits data when learn_ratio < 1 and returns testdata", {
  # A ratio less than 1 and check if test_data is non-null
  set.seed(1)
  data <- data.frame(matrix(rnorm(900), ncol = 5))
  data$class <- sample(1:3, nrow(data), replace = TRUE)

  var_types <- rep("c", 5)  # 4 numeric columns + 1 class
  result <- vineclass(data, learn_ratio = 0.5, var_types = var_types, nCore = 1 )

  # We expect some split between train_data and test_data
  expect_true(nrow(result$train_data) > 0)
  expect_true(nrow(result$test_data) > 0)
})

test_that("vineclass estimates class priors from the learning data", {
  set.seed(11)
  data <- data.frame(matrix(rnorm(400), ncol = 2))
  data$class <- c(rep(1, 150), rep(2, 50))

  result <- vineclass(
    data = data,
    learn_ratio = 0.7,
    var_types = c("c", "c"),
    nCore = 1
  )

  expected_priors <- prop.table(table(result$train_data[, ncol(result$train_data)]))
  expect_equal(unname(result$class_priors), unname(as.numeric(expected_priors)))
})

test_that("test_vineclass uses stored class priors", {
  set.seed(12)
  data1 <- data.frame(x1 = rnorm(60), x2 = rnorm(60))
  data2 <- data.frame(x1 = rnorm(60), x2 = rnorm(60))
  data1$class <- 1
  data2$class <- 2
  data <- rbind(data1, data2)

  result <- vineclass(
    data = data,
    learn_ratio = 1,
    var_types = c("c", "c"),
    nCore = 1
  )

  test_data <- data.frame(x1 = c(-0.5, 0, 0.5), x2 = c(-0.5, 0, 0.5))
  post_default <- test_vineclass(test_data, result, c("c", "c"))

  result_shifted <- result
  result_shifted$class_priors <- c(0.95, 0.05)
  post_shifted <- test_vineclass(test_data, result_shifted, c("c", "c"))

  expect_false(isTRUE(all.equal(post_default, post_shifted)))
})
