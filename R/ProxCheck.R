
proxCheck <- function(avec, bvec) {
  m <- length(avec)
  mx <- max(c(avec, bvec))
  avec <- c(avec, mx + 1)
  bvec <- c(bvec, mx + 2)
  icond1 <- 0
  icond2 <- 0
  idifv <- rep(0, 2 * m + 1)
  i1 <- 1
  i2 <- 1
  ncom <- 0
  ndif <- 0
  icomv <- rep(0, m)
  a <- avec[i1]
  b <- bvec[i2]

  while (i1 <= m & i2 <= m) {
    if (a == b) {
      ncom <- ncom + 1
      icomv[ncom] <- a
      i1 <- i1 + 1
      i2 <- i2 + 1
      a <- avec[i1]
      b <- bvec[i2]
    } else if (a < b) {
      ndif <- ndif + 1
      idifv[ndif] <- a
      i1 <- i1 + 1
      a <- avec[i1]
    } else {
      ndif <- ndif + 1
      idifv[ndif] <- b
      i2 <- i2 + 1
      b <- bvec[i2]
    }
  }

  if (i2 > m & i1 <= m) {
    for (i in i1:m) {
      ndif <- ndif + 1
      idifv[ndif] <- avec[i]
    }
  }
  if (i1 > m & i2 <= m) {
    for (i in i2:m) {
      ndif <- ndif + 1
      idifv[ndif] <- bvec[i]
    }
  }
  if (ndif == 2) {
    icond1 <- idifv[1]
    icond2 <- idifv[2]
    iok <- 1
  } else {
    iok <- 0
  }
  list(iok = iok, icomv = icomv, icond1 = icond1, icond2 = icond2)
}
