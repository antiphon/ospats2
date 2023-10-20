#' Legacy Ospats Algorithm
#'
#' Refactored Ospats::ospatsF
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
#' @details The algorithm in [ospatsF] with new inputs and outputs. Only change is that
#' the covariance matrix can be given to avoid the exponential covariance assumption and the "mean" covariance formula
#' used in the original scripts.
#'
#' See the vignette \code{vignette("legacy",)}
#'
#'
#'
#' @export

ospatsF_ref <- function(x,
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
  k1 <- kmeans(xy,
               centers = nstrata,
               nstart = 10,
               iter.max = 500,
               algorithm = "MacQueen")
  strat0 <- k1$cluster
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
  # Target
  Obj <- sum(sqrt(OA2_init))
  Obarbest <- Inf
  strat_best <- rep(NA, n)
  OA2_best <- NULL
  ##### Main looping structure
  #
  # outer loop to rerun many times. From same initials
  for(outer_it in 1:niter_outer) {
    OA2    <- OA2_init
    strat  <- strat0
    # Loop for a kind of sgd, one run starting from initial values and
    # each iteration checking if a unit would be better elsewhere.
    # Note: No switching, so no control for strata size > 0.
    for(it in 1:niter) {
      transfers <- 0
      # Random order pass over all units
      u <- sample(n)
      for(i in u) {
        # would unit i be better in some other strata?
        si        <- strat[i]
        s_not_i   <- (1:nstrata)[-si]
        OA2_no_i  <- max(0, OA2[si] - sum( D2[i, strat == si] )) # numerical precision might throw negative
        OA2_add_i <- sapply(s_not_i, \(sj) OA2[sj] + sum(D2[i, strat == sj]))
        Odelta    <- sqrt(OA2_no_i) - sqrt(OA2[si]) + sqrt(OA2_add_i) - sqrt(OA2[-si])
        # Come up with change probabilities
        Odelta[Odelta < (-Obj * 1e-10)] <- 0 # TR: greedy
        pr <- exp(-abs(Odelta) / temperature) # annealing-type
        # TR: The original alg flips the coin in order over strata, so smaller preferred.
        if(any(is.na(pr))) browser()
        for(j in seq_along(pr)) {
          if( runif(1) < pr[j]  ) {  # accept change
            transfers <- transfers + 1
            strat[i]  <- s_not_i[j]
            OA2[si]   <- OA2_no_i
            OA2[s_not_i[j]] <- OA2_add_i[j]
            Obj       <- sum(sqrt(OA2))
            # to get out of the for(j) loop
            break
          }
        }
      } # eo going over each unit in random order and moving them.
      # lower temperature for less fuzzy jump probabitilites
      temperature <- temperature * coolingrate
      # exit condition: nothing changed. TR: Why exit here?
      cat2(sprintf("    tras'rs [%i]", transfers))
      if(transfers == 0) break
    }  # single run ends
    # Did we get a better one this round?
    Obj <- sum(sqrt(OA2))
    Obarfinal <- Obj/n
    if(Obarfinal < Obarbest) {
      Obarbest   <- Obarfinal
      strat_best <- strat
      OA2_best    <- OA2
    }
    cat1(sprintf("Run %4i: Objective function [%7.4f]",
                 outer_it, Obarfinal))
  } # outer iterations
  # done
  list(stratification = strat_best,
             Obar = Obarbest, OA = sqrt(OA2_best))
}
