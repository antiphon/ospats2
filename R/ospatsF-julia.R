#' Ospats Julia Version in R
#'
#' Translate the Julia version of Ospats to R
#' @param data data matrix with columns [x y pred var] in that specific order
#' @param dRange Range parameter of the exponential correlation model
#' @param dRSquare R^2 parameter for adjusting for "regression-to-the-mean" (see the paper)
#' @param dMaxrun Number of independent (apart from same initial states) repeated optimisation runs
#' @param nCycles Number of iterations per one run of the optimisation algorithm
#' @param dStrata Number of strata to have
#' @param verbose Print runtime diagnostics
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
#' @export

ospatsF_julia <- function(
                          data,
                          dRange,
                          dRSquare,
                          dMaxrun,
                          nCycles,
                          dStrata,
                          verbose = FALSE
                        )
{
  x       <- data[,1] #X
  y       <- data[,2] #Y
  z_pred  <- data[,3] # prediction
  s2      <- data[,4] # prediction variance
  N      <- nrow(data)

  # name
  H      <- dStrata
  maxcycle <- nCycles
  #
  Z2   <- outer(z_pred, z_pred, "-")^2 / dRSquare
  Lags <- dist(cbind(x,y)) |> as.matrix()
  SS   <- outer(s2, s2, "+")
  Q    <- exp(-3.0 * Lags / dRange)
  d2 <- Z2 + SS *( 1 - Q )

  TOTd2  <- sum(d2) / 2
  ObarH1 <- sqrt(TOTd2) / N
  #
  ObarS  <- NULL #array(0, maxrun2)
  StratS <- NULL#array(NA, c(N, maxrun2))
  cbOjb  <- NULL #array(NA, H)
  cbObjS <- NULL #array(NA, c(H, maxrun2))

  cat2 <- cat
  println <- if(verbose) function(...) cat2(..., "\n") else \(...) NULL

  # Repeated runs of the algorithm?
  for( run in  1:dMaxrun ) {
    TotTransf <- 0
    println("** Start run nr ", run)

    ########## Section 2. Initial stratification
    missing <- N - H * floor(N/H)
    A <- 1:H
    B <- c(A, A)
    nrepeat <- N/H - 2
    for(rep in 1:nrepeat) B <- c(B, A)
    if(missing>0) B <- c(B, 1:missing)
    strat0 <- array(NA, N)
    w <- sample(N)
    strat0[w] <- B
    # Initial random stratification done. (TR could have been done soooo much simpler)
    # strat0 <- sample(1:H, N, replace = TRUE)

    #### Section 3: Contributions from Strat to Objective Function
    Sd2 <- array(1, H)
    for(strat in 1:H) {
      Sd2[strat] <- 0
      for(i in 1:(N-1)) {
        if(strat0[i] == strat) {
          for(j in (i+1):N) {
            if(strat0[j] == strat) {
              Sd2[strat] <- Sd2[strat] + d2[i, j]
            }
          }
        }
      }
    }
    #Sd2Init <- Sd2
    cbObj <- sqrt(Sd2)
    O <- sum(cbObj)
    ObarInit <- O/N
    println("ObarInit =", ObarInit)

    #### Section 4: Transfer Grid points
    println("\n ------ transferring grid points ------ \n")

    stratcy <- strat0
    TotTransf <- 0
    TotCycle  <- 0

    for(cycle in 1:maxcycle) {
      transfers <- 0
      u <- sample(N)
      for(t in u) {
        Delta  <- 0
        change <- 0
        A <- stratcy[t]
        ij <- which( stratcy == A )
        dA <- sum( d2[t, ij] )
        sumd2tinA <- dA
        Sd2Amint <- Sd2[A] - sumd2tinA
        cbObjA <- sqrt(abs(Sd2Amint))

        for(stratnr in 1:H) {
          Delta <- 0 # red
          sumd2plus <- 0 # red
          if (stratnr != A) {
            B <- stratnr
            ij <- which(stratcy == B)
            dB <- sum(d2[t, ij])
            sumd2plus <- dB
            cbObjB <- sqrt(abs(Sd2[B] + sumd2plus))
            Delta <- cbObjA + cbObjB - cbObj[A] - cbObj[B]
            if(Delta < 0 * 1e-10) {
              change <- 1
              transfers <- transfers + 1
              stratcy[t] <- B
              Sd2[A] <- Sd2[A] - sumd2tinA
              Sd2[B] <- Sd2[B] + sumd2plus
              cbObj  <- sqrt(abs(Sd2))
              O <- sum(cbObj)
              Delta <- 0 # red
            }
          }
          if(change) break
        } # strantr
      } # t in u

      #println("cycle", cycle, "    transfers = ", transfers)
      TotTransf  = TotTransf + transfers
      if(transfers == 0) break
      TotCycle = cycle
    }  # cycle
    println("Total number of transfers = ", TotTransf)
    println("Number of iteration cycles = ", TotCycle)

    O <- sum(cbObj)
    ObarFinal <- O / N
    #browser()

    ##      save results from runs
    ObarS <- cbind(ObarS, ObarFinal)
    StratS <- cbind(StratS, stratcy)
    cbObjS <- cbind(cbObjS, cbObj)
  }
  # End of runs. What follows is summary.
  #  browser()
  # Some weird indexing issues:
  best <- order(ObarS)
  strat_best <- StratS[,best[1]]
  cbObj_best <- cbObjS[, best[1]]
  obar_best  <- ObarS[best[1]]
#
#   # Section 5. Sample sizes and relative standard error
#   Nh <- table(strat_best)
#   meanzpred <- mean(z_pred)
#   n_pred <- ( 100 * obar_best /(meanzpred * RSEmax) )^2
#   n_pred <- round(n_pred)
#
#   sum_ahOh <- sum(Nh * cbObj_best)
#
#   # Neyman allocation
#   nh <- n_pred * Nh * cbObj_best / sum_ahOh
#   nh <- round(nh)
#
#   # Relative standard error
#   RSE <- 100 * obar_best / ( meanzpred * sqrt(nmax) )
#
#   # Section 6.
#   n_tot <- sum(nh) |> round()
#
#   points <- NULL
#
#   for(h in 1:H) {
#     k <- which(strat_best == h)
#     stratsize <- length(k)
#     w <- sample(stratsize)
#     f <-  nh[h] # what what
#     k_rand <- k[w]
#     points_h <- k_rand[1:f]
#     points <- c(points, points_h)
#   }
#   xs <- x[points] # GLOBAL x and y !
#   ys <- y[points]
#   strata <- strat_best[points]
#   sampnr <- 1:n_tot
#   browser()
#   strs <- cbind(sampnr, strata, points, xs, ys)
#
#   # pretty much done
#
#   # Section 7.

#  write.csv(file = "Stratification_r", strat_best)
#  write.csv(file = "Sample_r", strs)

  # println("FINAL RESULTS: ")
  # println("Value of Obar (O/N) without stratification : ", ObarH1)
  # println("Values of Obar : ", ObarS)
  # println("Lowest value of Obar : ", obar_best)
  # println("Sizes of grid-strata : ", Nh)
  # println("Contributions to O from strata (cbObj) : ", cbObj_best)
  # println("Mean of z-predictions : ", meanzpred)
  # println("Predicted sample size needed to attain RSEmax : ", n_pred)
  # println("Neyman sample allocation : ", nh)
  # println("Predicted RSE for nmax : ", RSE)
  #
  # println("--------------------- END FUNCTION OSPATSMR ---")

  list(stratification = strat_best, Obar0 = ObarH1,
       Obar = obar_best,
       OA   = cbObj_best
  )
}






