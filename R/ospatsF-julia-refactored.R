#' Ospats Julia Version in R
#'
#' Translate the Julia version of Ospats to R
#' @param dat data, [x y pred var]
#'
#' @details The original Julia code:
#'
#' https://github.com/jjdegruijter/ospats/ , file "ospatsmr"
#'
#' The original code uses global variables; This version mimics the input naming of the original ospatsF.
#'
#' This code only computes the stratification, and not the Neyman allocation and sampling in "Section 6" of the Julia code,
#' nor does the code include "Section 7." which writes the results into a file.
#'
#' This "_ref" is a refactored version, similar to ospatsF_ref, notably the input names change.
#'
#' @export

ospatsF_julia_ref <- function(
                          x,
                          covmodel_range,
                          nstrata,
                          rsquared,
                          niter,
                          niter_outer,
                          verbose = FALSE
                        )
{
  cat2 <- if(verbose) message else \(x) NULL

  n  <- nrow(x)
  xy <- cbind(x[,"x", drop = TRUE], x[,"y", drop = TRUE])
  z  <- x[,"pred", drop = TRUE]
  v  <- x[,"var", drop = TRUE]
  #
  ## The spatially correlated deviations
  z2  <- outer(z, z, "-")^2
  sv2 <- outer(v, v, "+")
  lag <- dist(xy) |> as.matrix()
  Q   <- exp(-3*lag/covmodel_range)  ## Check formula? This does not make sense.
  D2  <- z2/rsquared + sv2 * (1 - Q) # dev matrix

  TOTd2  <- sum(D2) / 2
  ObarH1 <- sqrt(TOTd2) / n
  #
  ObarS  <- NULL
  StratS <- NULL
  cbOjb  <- NULL
  cbObjS <- NULL

  # Repeated runs of the algorithm?
  for( run in  1:niter_outer ) {
    ########## Section 2. Initial stratification
    strat <- sample(1:nstrata, n, replace = TRUE)
    #### Section 3: Contributions from Strat to Objective Function
    OA2   <- tapply(1:n, strat, \(i) sum( D2[i,i] ))  / 2
    #### Section 4: Transfer Grid points
    for(cycle in 1:niter) {
      transfers <- 0
      u         <- sample(n)
      for(i in u) {
        # would unit i be better in some other strata?
        si        <- strat[i]
        s_not_i   <- (1:nstrata)[-si]
        OA2_no_i  <- OA2[si] - sum( D2[i, strat == si] )
        OA2_add_i <- sapply(s_not_i, \(sj) OA2[sj] + sum(D2[i, strat == sj]))
        Odelta    <- sqrt(OA2_no_i) - sqrt(OA2[si]) + sqrt(OA2_add_i) - sqrt(OA2[-si])
        # Julia version is simply greedy
        for(j in seq_along(Odelta)) { # note the ordered check
          if( Odelta[j]  < 0 ) {  # accept change
            transfers <- transfers + 1
            strat[i]  <- s_not_i[j]
            OA2[si]   <- OA2_no_i
            OA2[s_not_i[j]] <- OA2_add_i[j]
            # to get out of the for(j) loop
            break
          }
        }
      } # eo going over each unit in random order and moving them.
      cat2(sprintf("    tras'rs [%i]", transfers))
      if(transfers == 0) break # Why?
    }  # cycle
    O <- sum(sqrt(OA2))
    ObarFinal <- O / n
    ##      save results from run
    ObarS <- c(ObarS, ObarFinal)
    StratS <- cbind(StratS, strat)
    cbObj <- sqrt(OA2)
    cbObjS <- cbind(cbObjS, cbObj)
    cat2(sprintf("Run %4i: Objective function [%7.4f]",
                 run, Obarfinal))
  }
  # End of runs. What follows is summary.
  #  browser()
  # Some weird indexing issues:
  best <- order(ObarS)
  obar_best  <- ObarS[best[1]]
  strat_best <- StratS[,best[1]]
  cbObj_best <- cbObjS[, best[1]]

  list(stratification = strat_best, Obar = obar_best, OA = cbObj_best)
}






