
findleaf <- function(conditioned, conditioning, iprint = FALSE) {
  given <- unique(conditioning)
  uniq <- sort(unique(conditioned))
  tem <- table(conditioned)
  i0 <- (tem == 1)
  leafs <- setdiff(uniq[i0], given)
  if (iprint) cat("findleaf: ", leafs, "\n")
  leafs
}
