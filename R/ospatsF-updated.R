#' Modified Ospats Algorithm
#'
#' Ospats-algorithm with some modifications.
#'
#' @param x data, with variables [x y pred var]
#' @param covmodel_range Range parameter for assumed exponential correlation model
#' @param nstrata number of starta to consider
#' @param niter_outer Number of independent runs of the optimisation algorithm
#' @param niter Number of iterations per one run of the optimisation algorithm
#' @param verbose Print runtime diagnostics? (0,1,2)
#' @param rsquared The R^2 in the paper (default: 1)
#' @param temperature Annealing factor, will accept slightly bad moves with prob ~exp(-abs(delta)/temperature)
#' @param coolingrate Change temperature each interation by this factor. Should be at most 1.
#' @param Cov Optional, overrides the covariance matrix calculation using exp-correlation. No checks with data variances.
#'
#' @details Changes:
#'
#' 1. the covariance matrix can be given to avoid the exponential covariance assumption and the "mean" covariance formula
#' used in the original scripts.
#'
#' 2. The initial stratification is based on cum-root-f, not k-means which is unreliable.
#'
#'
#' @export

ospats2 <- function(x,
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
  cat1 <- if(verbose>0) message else \(x) NULL
  cat2 <- if(verbose>1) message else \(x) NULL

  n <- nrow(x)
  xy <- cbind(x[,"x", drop = TRUE], x[,"y", drop = TRUE])
  z  <- x[,"pred", drop = TRUE]
  ## Starting configuration
  strat0     <- cumrootf(z, nstrata)
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
  # run
  res <- ospats2_c(D2,
                      strat0-1, # go base 0
                      OA2_init,
                      niter, niter_outer,
                      temperature, coolingrate, verbose )
  #return(res)
  OA         <- sqrt( res[["OA2_best"]] )
  strat_best <- res[["strat_best"]] + 1 # base 0 adjustment

  # done
  list(stratification = strat_best,
       Obar = sum(OA)/n,
       OA = OA)
}
