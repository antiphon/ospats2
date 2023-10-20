#' Naive Cum-Root-F stratification
#'
#' @export
cumrootf <- function(y, nstrata, nbins = 1000) {
  m <- min(y)
  M <- max(y)
  m <- m - (M-m)*1e-5
  M <- M + (M-m)*1e-5

  bins <- seq(m, M, l = nbins)
  ix <- findInterval(y, bins)
  f  <- tabulate(ix, nbins)
  csf  <- cumsum(sqrt(f))
  mx <- csf[nbins] / nstrata
  bw <- mx * (0:nstrata)
  fi <- findInterval(csf, bw)
  ib <- which( diff( fi ) != 0 )[-nstrata]
  str <- findInterval(y, bins[c(1, ib, nbins) ])
  str
}
