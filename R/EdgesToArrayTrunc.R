
edges2arrayTrunc <- function(d, mtree, edgemat, pcv, iprint = FALSE) {
  A <- matrix(0, d, d)
  pc <- matrix(NA, d, d)
  mat <- edgemat

  ntot <- 0
  lastcol <- d

  while (ntot < d - mtree - 1) {
    ind <- which(mat[, 3] == mtree)
    conditioned <- mat[ind, 1:2]
    conditioning <- mat[ind, 4:ncol(mat)]
    leafs <- findleaf(conditioned, conditioning, iprint)
    leafs <- rev(leafs) # Reverse order for leaf with the largest index.
    nleaf <- length(leafs)
    ntot <- ntot + nleaf
    for (i in seq_len(nleaf)) {
      icol <- lastcol - i + 1
      if (icol == (mtree + 1)) break
      w <- leafs[i]
      A[icol, icol] <- w
      for (ii in seq(mtree, 1, by = -1)) {
        itree <- which(mat[, 3] == ii)
        for (jj in itree) {
          cs <- mat[jj, ]
          if (cs[1] == w) {
            A[ii, icol] <- cs[2]
            break
          } else if (cs[2] == w) {
            A[ii, icol] <- cs[1]
            break
          }
        }
        pc[ii, icol] <- pcv[jj]
        mat[jj, ] <- 0
      }
    }
    lastcol <- icol - 1
  }

  # Fill remaining columns as usual.
  for (j in (mtree + 1):2) {
    k <- which(mat[, 3] == (j - 1))
    w <- mat[k, 2]
    A[j, j] <- w
    A[j - 1, j] <- mat[k, 1]
    pc[j - 1, j] <- pcv[k]
    if (j == 2) {
      A[1, 1] <- A[1, 2]
    } else {
      for (ii in (j - 2):1) {
        itree <- which(mat[, 3] == ii)
        for (jj in itree) {
          cs <- mat[jj, ]
          if (cs[1] == w) {
            A[ii, j] <- cs[2]
            break
          } else if (cs[2] == w) {
            A[ii, j] <- cs[1]
            break
          }
        }
        pc[ii, j] <- pcv[jj]
        mat[jj, ] <- 0
      }
    }
    mat[k, ] <- 0
  }
  list(VineA = A, pcmat = pc)
}
