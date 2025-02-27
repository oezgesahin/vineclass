
edges2array <- function(d, edgemat, pcv) {
  A <- matrix(0, d, d)
  pc <- matrix(NA, d, d)
  mat <- edgemat # This matrix will have rows deleted as columns of A are filled.

  for (j in d:2) {
    k <- which(mat[, 3] == (j - 1))  # Find the row corresponding to the tree level.
    w <- mat[k, 2]
    A[j, j] <- w
    A[j - 1, j] <- mat[k, 1] # Assign edge variables.
    pc[j - 1, j] <- pcv[k]
    if (j == 2) {
      A[1, 1] <- A[1, 2]
    } else {
      for (ii in (j - 2):1) {
        itree <- which(mat[, 3] == ii)
        for (jj in itree) {
          cs <- mat[jj, ]
          pctem <- pcv[jj]
          if (cs[1] == w) {
            A[ii, j] <- cs[2]
            break
          } else if (cs[2] == w) {
            A[ii, j] <- cs[1]
            break
          }
        }
        pc[ii, j] <- pctem
        mat[jj, ] <- 0
      }
    }
    mat[k, ] <- 0
  }
  list(VineA = A, pcmat = pc)
}
