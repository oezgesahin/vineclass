#' Resolve class prior probabilities
#' @noRd
resolve_class_priors <- function(labels = NULL, total_cl, class_priors = NULL) {
  if (is.null(class_priors)) {
    if (is.null(labels)) {
      stop("Either 'labels' or 'class_priors' must be provided.")
    }
    class_counts <- tabulate(as.integer(labels), nbins = total_cl)
    class_priors <- class_counts / sum(class_counts)
  }

  if (!is.numeric(class_priors) || length(class_priors) != total_cl) {
    stop("'class_priors' must be a numeric vector with one entry per class.")
  }
  if (any(!is.finite(class_priors)) || any(class_priors < 0)) {
    stop("'class_priors' must contain finite, non-negative values.")
  }

  prior_sum <- sum(class_priors)
  if (prior_sum <= 0) {
    stop("'class_priors' must sum to a positive value.")
  }

  class_priors / prior_sum
}

#' Convert class-conditional likelihoods to posterior probabilities
#' @noRd
posterior_from_likelihoods <- function(likelihoods, class_priors) {
  weighted_likelihoods <- sweep(as.matrix(likelihoods), 2, class_priors, FUN = "*")
  likelihood_sums <- rowSums(weighted_likelihoods)
  likelihood_sums[is.na(likelihood_sums) | likelihood_sums <= 0] <- 1e-5
  weighted_likelihoods / likelihood_sums
}
