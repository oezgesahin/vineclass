---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->



# vineclass: Vine Copula-Based Classifiers

<!-- badges: start -->
<!-- badges: end -->

vineclass provides tools for vine copula-based classification for a mixed of continuous-ordinal variables. It splits data into learning (training) and test sets, fits vine copulas per class on the training subset, and computes classification summaries (including posterior probabilities) on the test set or new data.

## Installation

You can install the development version from [GitHub](https://github.com/oezgesahin) with:

``` r
# install.packages("remotes")
remotes::install_github("oezgesahin/vineclass")
```
## Package overview

Below is an overview of some functions and features. 

* ```vineclass()```: main interface for vine copula classification for a mixed of continuous-ordinal variables. Splits data, transforms it to u-scale nonparametrically, fits class-specific vine models, and returns posterior probabilities.

* ```test_vineclass()```: calculates posterior probabilities on new test data using fitted vine models from \code{vineclass()}.

## Fitting with continuous variables

```r
# Simulate data for two classes
set.seed(123)
n1 <- 500   # number of observations for class 1
n2 <- 500   # number of observations for class 2
p <- 4      # number of variables

data1 <- matrix(rnorm(n1 * p), nrow = n1, ncol = p)
data2 <- matrix(rnorm(n2 * p, mean = 1), nrow = n2, ncol = p)

# Create data frames and add class labels
df1 <- data.frame(data1)
df1$class <- 1
df2 <- data.frame(data2)
df2$class <- 2
data <- rbind(df1, df2)
data <- data[sample(nrow(data)), ]

# Define variable types: all variables are continuous
# Note: var_types should have length equal to number of predictors (columns minus the class label)
var_types <- rep("c", p)

# Run vine-based classification using all data for training, where pair copulas are parametric.
result <- vineclass(
  data        = data,
  learn_ratio = 1,
  var_types   = var_types,
  copfam      = "parametric"
)

# Posterior probabilities for learning data
print(head(result$post_prob_learn, 5))
#>             [,1]       [,2]
#> [1,] 0.183429543 0.81657046
#> [2,] 0.255243586 0.74475641
#> [3,] 0.903166318 0.09683368
#> [4,] 0.005886429 0.99411357
#> [5,] 0.156375672 0.84362433
```

## Assess classification accuracy

```r
predicted_train <- max.col(result$post_prob_learn)
true_train <- result$train_data[, ncol(result$train_data)]
accuracy <- mean(predicted_train == true_train)
cat("Classification accuracy in training data:", accuracy, "\n")
#> Classification accuracy in training data: 0.838
```

## Testing the classifier

```r
# Assume external data for testing the classifier
test_data <- matrix(rnorm(100),ncol = 4)
var_types_test <- rep("c", p)
test_post_prob <- test_vineclass(test_data, result, var_types_test)

# Posterior probabilities for test data
head(test_post_prob)
#>           [,1]       [,2]
#> [1,] 0.9660377 0.03396233
#> [2,] 0.9046306 0.09536936
#> [3,] 0.5074505 0.49254948
#> [4,] 0.9770398 0.02296016
#> [5,] 0.7177078 0.28229222
#> [6,] 0.2938357 0.70616429
```


## Fitting with a mixed of continuous-ordinal variables

```r
set.seed(42)
# For Class 1:
# - Variable 1, 3, 4: Continuous from N(0, 1)
# - Variable 2: Discrete with 3 levels, probabilities 0.3, 0.4, 0.3
data1 <- data.frame(
  x1 = rnorm(n1, mean = 0, sd = 1),
  x2 = sample(1:3, n1, replace = TRUE, prob = c(0.3, 0.4, 0.3)),
  x3 = rnorm(n1, mean = 0, sd = 1),
  x4 = rnorm(n1, mean = 0, sd = 1)
)
data1$class <- 1

# For Class 2:
# - Variable 1, 3, 4: Continuous from N(1, 1)
# - Variable 2: Discrete with 3 levels, probabilities 0.2, 0.5, 0.3
data2 <- data.frame(
  x1 = rnorm(n2, mean = 1, sd = 1),
  x2 = sample(1:3, n2, replace = TRUE, prob = c(0.2, 0.5, 0.3)),
  x3 = rnorm(n2, mean = 1, sd = 1),
  x4 = rnorm(n2, mean = 1, sd = 1)
)
data2$class <- 2
data <- rbind(data1, data2)
data <- data[sample(nrow(data)), ]


# Define variable types:
# There are 4 variables:
# - x1, x3, x4 are continuous ("c")
# - x2 is discrete ("d")
var_types <- c("c", "d", "c", "c")

# Run vine-based classification with 75% of data for training and 25% for testing.
result_mix <- vineclass(
  data        = data,
  learn_ratio = 0.75,
  var_types   = var_types,
  copfam      = "parametric"
)

# Posterior probabilities for learning data
print(head(result_mix$post_prob_learn, 5))
#>            [,1]      [,2]
#> [1,] 0.22714657 0.7728534
#> [2,] 0.36879169 0.6312083
#> [3,] 0.22938433 0.7706157
#> [4,] 0.70136309 0.2986369
#> [5,] 0.03549832 0.9645017
```

## Contact

Please contact O.Sahin@tudelft.nl if you have any questions.

## References

Sahin, {\"O}., \&  Joe, H. (2024), Vine Copula-Based Classifiers with Applications. Journal of Classification, 1-29. [article](https://link.springer.com/article/10.1007/s00357-024-09494-y)

Brechmann, E.C., \&  Joe, H. (2015), Truncation of vine copulas using fit indices. Journal of Multivariate Analysis, 138, 19-33. [article](https://doi.org/10.1016/j.jmva.2015.02.012)
