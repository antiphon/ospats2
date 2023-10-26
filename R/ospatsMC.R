#' Ospats Algorithm Using Monte Carlo
#'
#' Ospats-algorithm with Monte Carlo samples.
#'
#' @param Z A matrix of n x M simulated values of the response
#' @param nstrata number of starta to consider
#' @param niter_outer Number of independent runs of the optimisation algorithm
#' @param niter Number of iterations per one run of the optimisation algorithm
#' @param verbose Print runtime diagnostics? (0,1,2)
#' @param rsquared The R^2 in the paper (default: 1)
#' @param temperature Annealing factor, will accept slightly bad moves with prob ~exp(-abs(delta)/temperature)
#' @param coolingrate Change temperature each interation by this factor. Should be at most 1.
#'
#' @details Estimate the mean of the target using Monte Carlo sample of size M.
#'
#' The initial stratification is based on cum-root-f, not k-means which is unreliable.
#'
#' @useDynLib ospats2
#' @import Rcpp
#' @export

ospatsMC <- function(Z,
                    nstrata,
                    niter = 100,
                    niter_outer = 100,
                    verbose = 0,
                    # Some parameters for annealing
                    temperature = 1,
                    coolingrate = .95,
                    rsquared = 1,
                    v = 1
                    ) {
  if(!is.matrix(Z)) stop("Z needs to be a matrix, one row of simulations per observation.")
  ## Starting configuration
  z_mean     <- rowMeans(Z)
  strat0     <- cumrootf(z_mean, nstrata)
  #
  # Initial penalty terms
  Sh_init <- sapply(split(1:nrow(Z),
                           strat0),
                     \(i)  mean(apply(Z[i,, drop=FALSE], 2, sd)) * sqrt(length(i) * (length(i)-1))  )

  ## The spatially correlated deviations
  # run
  #
  calcr <- if(v == 1) ospats2_MC_c else ospats2_MC_v2_c
  res <- calcr(Z,
                      strat0-1, # go base 0
                      Sh_init,
                      niter,
                      niter_outer,
                      temperature,
                      coolingrate,
                      verbose )
  #return(res)
  Sh         <- res[["Sh_best"]]
  strat_best <- res[["strat_best"]] + 1 # base 0 adjustment

  # done
  list(stratification = strat_best,
       Obar = sum(Sh)/nrow(Z),
       Sh = Sh)
}
