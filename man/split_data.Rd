% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/split_data.R
\name{split_data}
\alias{split_data}
\title{Data splitting}
\usage{
split_data(data, learn_rat, valid_rat = NA)
}
\arguments{
\item{data}{A data frame or matrix to be split.}

\item{learn_rat}{A numeric value specifying the proportion of data to allocate to the learning (training) set.}

\item{valid_rat}{Optional numeric value specifying the proportion of data to allocate to the validation set.
If not provided, no validation set will be created.}
}
\value{
A list containing:
\itemize{
\item \code{learning}: The learning (training) set.
\item \code{valid}: The validation set (if \code{valid_rat} is provided).
\item \code{test}: The test set.
}
}
\description{
Splits the input dataset into learning (training), validation (optional),
and test sets based on the specified proportions.
}
\examples{
set.seed(1)
data <- matrix(rnorm(900), ncol = 3)

# Split into training (60\%), validation (20\%), and test (20\%)
split <- split_data(data, learn_rat = 0.6, valid_rat = 0.2)
train <- split$learning
valid <- split$valid
test <- split$test

# Split into training (70\%) and test (30\%) without validation
split <- split_data(data, learn_rat = 0.7)
train <- split$train
test <- split$test
}
