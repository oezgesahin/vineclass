isPosDef <- function(amat) {
  tt <- try(chol(amat), silent = TRUE)
  "matrix" %in% class(tt)
}

norm_score <- function(data) {
  n_col <- ncol(data)
  n_row <- nrow(data)
  emp_scrs <- matrix(0, n_row, n_col)
  qn <- stats::qnorm(((1:n_row) - 0.5) / n_row)
  for (j in seq_len(n_col)) {
    temp <- rank(data[, j])
    z_scr <- qn[temp]
    emp_scrs[, j] <- z_scr
  }
  emp_scrs
}
