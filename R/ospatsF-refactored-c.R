#' Legacy Ospats Algorithm
#'
#' Refactored Ospats::ospatsF and using RCpp
#'
#' @param x data, with variables [x y pred var]
#' @param covmodel_range Range parameter for assumed exponential correlation model
#' @param nstrata number of starta to consider
#' @param niter_outer Number of independent runs of the optimisation algorithm
#' @param niter Number of iterations per one run of the optimisation algorithm
#' @param verbose Print runtime diagnostics?
#' @param rsquared The R^2 in the paper (default: 1)
#' @param temperature Annealing factor, will accept slightly bad moves with prob ~exp(-abs(delta)/temperature)
#' @param coolingrate Change temperature each interation by this factor. Should be at most 1.
#' @param Cov Optional, overrides the covariance matrix calculation using exp-correlation. No checks with data variances.
#'
#' @details The algorithm in [ospatsF] with new inputs and outputs. Only change is that
#' the covariance matrix can be given to avoid the exponential covariance assumption and the "mean" covariance formula
#' used in the original scripts.
#'
#' See the vignette \code{vignette("legacy",)}
#'
#'
#' @useDynLib ospats2
#' @export

ospatsF_refc <- function(x,
                          covmodel_range,
                          nstrata,
                          niter = 100,
                          niter_outer = 100,
                          verbose = 0,
                          # Some parameters for annealing
                          temperature = 1,
                          coolingrate = .95,
                          rsquared = 1,
                          Cov
) {

  n <- nrow(x)
  xy <- cbind(x[,"x", drop = TRUE], x[,"y", drop = TRUE])
  z  <- x[,"pred", drop = TRUE]
  ## Starting configuration
  k1 <- kmeans(xy,
               centers = nstrata,
               nstart = 10,
               iter.max = 500,
               algorithm = "MacQueen")
  strat0 <- k1$cluster - 1
  #
  ## The spatially correlated deviations
  z2  <- outer(z, z, "-")^2
  v   <- x[,"var", drop = TRUE]
  sv2 <- outer(v, v, "+")
  lag <- dist(xy) |> as.matrix()
  if(missing(Cov))
    Cov <- 0.5 * sv2 * exp(-3*lag/covmodel_range) ## Check formula? This does not make sense.
  D2  <- z2/rsquared + sv2 - 2 * Cov # dev matrix
  #
  # starting cost
  OA2_init <- sapply(split(1:n, strat0), \(i)
            sum( D2[i,i] ))  / 2

  res <- ospats_ref_c(D2, strat0, OA2_init, niter, niter_outer, temperature, coolingrate, verbose )
  #return(res)
  OA         <- sqrt( res[["OA2_best"]] )
  strat_best <- res[["strat_best"]] + 1 # base 0 adjustment
  # done
  list(stratification = strat_best,
             Obar = sum(OA)/n,
               OA = OA)
}
