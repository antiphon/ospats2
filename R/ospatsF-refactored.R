#' Legacy Ospats Algorithm
#'
#' Refactored Ospats::ospatsF
#'
#' @param x data, with named columns (x, y, pred, var)
#' @param nstrata number of strata
#' @param covmodel_range Exponential covariance model range parameter
#' @param niter number of iterations of the allocation shuffle
#' @param
#'
#'
#' @export

ospatsF_ref <- function(x,
                        covmodel_range,
                        nstrata,
                        niter = 100,
                        niter_outer = 100,
                        verbose = FALSE,
                        # Some parameters for annealing
                        temperature = 1,
                        coolingrate = .95
) {
  cat2 <- if(verbose) message else \(x) NULL

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
  S   <- 0.5 * sv2 * exp(-3*lag/covmodel_range)  ## Check formula?
  D2  <- z2 + sv2 - 2 * S # dev matrix
  #
  # starting cost
  OA2_init <- sapply(split(1:n, strat0), \(i)
            sum( D2[i,i] ))  / 2
  # Target
  Obj <- sum(sqrt(OA2_init))
  Obarbest <- Inf
  strat_best <- rep(NA, n)
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
        OA2_no_i  <- OA2[si] - sum( D2[i, strat == si] )
        OA2_add_i <- sapply(s_not_i, \(sj) OA2[sj] + sum(D2[i, strat == sj]))
        Odelta    <- sqrt(OA2_no_i) - sqrt(OA2[si]) + sqrt(OA2_add_i) - sqrt(OA2[-si])
        # Come up with change probabilities
        Odelta[Odelta < (-Obj * 1e-10)] <- 0 # TR: greedy
        pr <- exp(-abs(Odelta) / temperature) # annealing-type
        # TR: The original alg flips the coin in order over strata, so smaller preferred.
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
    }
    cat2(sprintf("Run %4i: Objective function [%7.4f]",
                 outer_it, Obarfinal))
  } # outer iterations
  # done
  data.frame(stratFrame = strat_best)
}
