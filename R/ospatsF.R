#' Ospats Function from Ospats-package
#'
#' Copy-paste of the original code (with small edits)
#'
#' @param data data matrix with columns [x y pred var] in that specific order
#' @param dRange Range parameter of the exponential correlation model
#' @param dRSquare R^2 parameter for adjusting for "regression-to-the-mean" (see the paper)
#' @param dMaxrun Number of independent (apart from same initial states) repeated optimisation runs
#' @param nCycles Number of iterations per one run of the optimisation algorithm
#' @param dStrata Number of strata to have
#' @param ClusterStart  external for dStart == 3 (Saby Input)
#' @param dStart choose between kMeans (0) or CumrootSquare(1) or external (3)
#' @param initialTemperature simulated annealing parameter
#' @param coolingRate  simulated annealing parameter
#' @param debug Useful during development for troubleshooting issues
#' @param verbose Print runtime diagnostics
#'
#' @details These are the refactoring changes to the original code to minimize the footprint:
#'
#' 1. The call \code{fastPdist2(Ar = as.matrix(z), Br = as.matrix(z)))^2} is now  \code{outer(z, z, "-")^2}
#' 1. The call \code{fastVsum(Ar = as.matrix(s2))} is now \code{outer(s2, s2, "+")}
#' 1. The call \code{fastPdist2(Ar = xy, Br = xy)} is now \code{as.matrix(dist(xy)) }
#'
#' Other than that the code has not been touched. Only this documentation has been added.
#'
#' @section Original documentation:
#' ## opstats algorithm:
#' As described in de Gruijter et al. (2015) {Optimizing Stratification and Allocation for
#' Design-Based Estimation of Spatial Means Using Predictions with Error}
#'
#' ## What does it do:
#' Given a map (for instance a digital soil map) with prediction variance this function will derive a
#' spatial stratifcation somewhere on the spectrum between and including compact geographical stratification
#' (high uncertainty) and clusters based on the quantiles of the data (low uncertainty).
#'
#' @export
ospatsF<- function(data = NULL, # input data (data frame). 4 columns [ X, Y, Pred, Var]
                   dRange = NULL, # spatial structure of prediction variance
                   nCycles = NULL, # Number of allocation iterations
                   dStart = NULL, # choose between kMeans (0) or CumrootSquare(1) or external (3)
                   ClusterStart = NULL , # external for dStart == 3 (Saby Input)
                   dMaxrun = NULL, # Number of runs the algorithm will go through to find optimal allocation
                   dRSquare = NULL, # Used for compensation of leveling
                   dStrata = NULL, # Number of strata
                   initialTemperature = NULL, #simulated annealing parameter
                   coolingRate = NULL,  # simulated annealing parameter
                   debug = FALSE, # Useful during development for troubleshooting issues
                   verbose = TRUE){  # Prints messages during the running of function


  ##################Embedded Function##############################################################
  #For clustering the predictions based on quantiles of the data
  cumsqfDel<- function(z, nclass, ns){
    # calculate the distribution based on nclass
    z2<- cbind(seq(1,length(z),by=1), z)
    d<- z2[order(z),2]
    io<- z2[order(z),1]
    nd<- length(d)
    dmin<- min(d)
    dmax<- max(d)
    step<- (dmax-dmin)/nclass
    h<- t(as.matrix(seq(dmin,dmax,by=step)))

    fx<- matrix(NA, nrow=nclass, ncol=2)
    fy<- matrix(NA, nrow=nclass, ncol=1)
    for (i in 1:nclass){
      ll<- h[1,i]
      ul<- h[1,i+1]
      ij=which(d>=ll & d<ul)
      if(i==nclass) {ij<- which(d>=ll & d<=ul)}
      nij<- length(ij)
      fx[i,1]<- ll
      fx[i,2]<- ul
      fy[i,]<- nij}

    sqf<- sqrt(fy)
    cumsqf<- cumsum(sqf)

    #set out boundaries based on ns (number of strata)
    cmax<- max(cumsqf)
    cmin<- min(cumsqf)
    cstep=(cmax-cmin) /ns
    cbound<- t(as.matrix(seq(cmin,cmax,by=cstep)))

    # find out corresponding class
    db<- matrix(NA, nrow=1, ncol= length(2:ns))
    for (i in 2:ns){
      cm<- abs(cumsqf-cbound[1,i])
      db[1,i-1]<- which(cm==min(cm))}

    idat<- matrix(0, nrow= nd, ncol= 1)
    xm<- matrix(0, nrow= ns, ncol=5)

    for (i in 1:ns){
      if(i==1) {ll<- 1} else {ll<- db[1,i-1]}
      if(i==ns) {ul<- nclass} else {ul<- db[1,i]}

      fll<- fx[ll,1]
      ful<- fx[ul,1]
      if(i==ns){ful<- fx[ul,2]}

      ij<- which(d>=fll & d<ful)
      if(i==ns) {ij<- which(d>=fll & d<=ful)}
      idat[ij,1]<- i
      dij<- d[ij]
      nij<- length(dij)

      xm[i,1]<- nij
      xm[i,2]<-mean(dij)
      xm[i,3]<-var(dij)
      xm[i,4]<-min(dij)
      xm[i,5]<-max(dij)}

    strat<- matrix(1,nrow=nd, ncol=1)
    strat[io,1]<- idat
    retval<- list(strat, xm)
    return(retval)}
  ###################################################################################################################


  # Define structure for storing time series of criterion
  params<- c(dRange,nCycles, dStart,dMaxrun, dRSquare,dStrata,initialTemperature,coolingRate  )
  Eall<-NULL
  # set initial temperature
  Temp <- initialTemperature

  x<- data[,1] #X
  y<- data[,2] #Y
  z<- data[,3] # prediction
  s2<-data[,4] # prediction variance
  s2s<- s2
  n = length(x)
  d20<- 0 # storing distance in case of re-run

  print("-----OSPATS INITIALISATION------")

  # generate initial solution
  print("Initial Solution")
  if (dStart == 0)  #use kmeans of coordinates
  {k1<- kmeans(x= cbind(x,y), centers= dStrata, nstart=10, iter.max = 500, algorithm = "MacQueen")
  strat0<- as.matrix(k1$cluster)}
  if (dStart == 1)  #use Cum-sqrt-f
  {nclass<-  100
  strat0<- cumsqfDel(z= z, nclass= nclass, ns= dStrata)[[1]]}
  if (dStart == 3)  #use external solution (saby)
  { strat0<- as.matrix(ClusterStart) }


  #                            INITIATION
  # Vectorised calculation of distance matrix, without loop
  print("distance matrix")
  if(d20==0){
    # calculate distance matrix
    xy = cbind(x,y)
    #Lag<- fastPdist2(Ar = xy, Br = xy)
    Lag <- dist(xy) |> as.matrix()

    # calculate prediction (z) difference (vectorized)
    #Z2<- (fastPdist2(Ar = as.matrix(z), Br = as.matrix(z)))^2 # TR: Mod 1
    Z2  <- outer(z, z, "-")^2

    Z2<-  Z2/dRSquare        #compensation for leveling


    # calculate variance (s2)
    #S2<- fastVsum(Ar = as.matrix(s2)) # TR: Mod 2
    S2 <- outer(s2, s2, "+")

    # calculate matrix of maximum of covariance
    print("matrix of maximum covariance")
    Cov_max <-  0.5 *S2
    Cov<-  Cov_max * exp(-3*Lag/dRange)
    d2<-  Z2 + S2 -2*Cov
    #d2 = S2 -2*Cov; % emulate equal z-pred's
    #d2 = Lag;    %use only Lag
    #d2 = Z2;     % use only Z (true or pred)
    d20<- d2} else {d2<- d20}

  rm(Lag, Z2, S2,Cov_max, Cov)
  gc()

  #format compact
  #MinZ2 = min(min(Z2)), MeanZ2 = mean(mean(Z2)), MaxZ2 = max(max(Z2));
  #MinS2 = min(min(S2)), MeanS2 = mean(mean(S2)), MaxS2 = max(max(S2));
  #MinCov = min(min(Cov)), MeanCov = mean(mean(Cov)), MaxCov = max(max(Cov));

  TOTd2<-  sum(colSums(d2))/2
  ObjNoStr<-  sqrt(TOTd2)
  ObarH1<-  ObjNoStr/n
  # Calculate sums of d2's within strata (Sd2) and contributions of strata to Obj (cbObj).
  print("D2 sum matrix")
  d2[which((is.na(d2)))]<- 0 ## Change NaNs to 0 (02/12/16)
  Sd2<- matrix(0, nrow=1, ncol=dStrata)
  for (strat in 1:dStrata){
    Sd2[1,strat]<-  0
    for (i in 1:(n-1)){
      if (strat0[i,1] == strat)
      {for (j in (i+1):n){
        if (strat0[j,1] == strat)
        {Sd2[1,strat]<-  Sd2[1,strat]+d2[i,j]}}}}}
  gc()
  if (debug==TRUE){print(Sd2)}
  Sd2_init = Sd2;

  cbObj<- sqrt(Sd2)
  Obj<- sum(cbObj)
  ObarInit<- Obj/n


  #                          START RUNS
  ObarBest<-  1000        #arbitrary large number to start with
  StratBest<- matrix(0, nrow=n, ncol=1)
  d1<- data

  for (run in 1:dMaxrun){
    if (debug==TRUE){print("_________________________________________________________________________________________________________RUN::")}
    TotTransf<- 0
    Sd2<- Sd2_init;
    #                             ITERATIVE RE-ALLOCATION
    if (verbose==TRUE){print(paste("----------------------------------------Run:", run, sep=" "))}
    stratcy<-  strat0       # stores stratum number for each point
    change<-  1
    for (cycle in 1:nCycles){ # loop through cycles


      if (debug==TRUE){print("________________________________________________________________________________________________________CYCLES:::")}
      transfers<-  0

      u<- t(as.matrix(sample.int(n, size = n, replace = FALSE, prob = NULL))) # put grid points in random order

      for (gp in 1:ncol(u)){

        if (debug==TRUE){print("________________________________________________________GP:::")}

        if (debug==TRUE){print(paste("le gp est",gp)) }

        samp <- u[1,gp]
        Delta <-  0
        change <-  0         # indicator of transfer
        A <-  stratcy[samp]
        # remove t from A
        ij <-  which(stratcy == A)
        dA <- sum(d2[samp,ij])-d2[samp,samp]
        sumd2tinA <-  dA
        Sd2Amint <-  Sd2[1,A] - sumd2tinA

        if (Sd2Amint<0) break() # test
        cbObjA <- sqrt(Sd2Amint)

        for (stratnr in 1:dStrata){ #% idem between t and the points in different strata


          if (debug==TRUE) {print("___________________________________________________STRATNR:::")}

          Delta<-  0
          sumd2plus<-  0
          if (stratnr != A){


            if (debug==TRUE){print("stratnr != A")}
            # Add to B
            B <-  stratnr
            ij <-  which(stratcy == B)
            dB <- sum(d2[samp,ij])
            sumd2plus<- dB

            cbObjB<- sqrt(Sd2[1,B] + sumd2plus)
            Delta<- cbObjA + cbObjB - cbObj[1,A] - cbObj[1,B]


            if (Delta < -Obj*1e-10){
              # always accept improvements
              if (debug==TRUE){print("Delta < -Obj*1e-10")}
              if (debug==TRUE){print("Delta < -Obj*1e-10")}
              if (debug==TRUE) print(paste("le Delta est " , Delta))
              pr <- 1
            } else {
              pr <- exp(-abs(Delta) / Temp ) # use Boltzmann to judge if deteriorations should be accepted
              if (debug==TRUE){print("!!!!Delta > -Obj*1e-10 !!!!")}
              if (debug==TRUE) print(paste("le Delta est " , Delta))
              if (debug==TRUE) print(paste("la temp?rature est " , Temp))
              if (debug==TRUE) print(paste("le crit?re de Boltzman " , pr))
            } # end calculate prob

            #browser()
            pChange <- runif(n = 1) # draw uniform deviate

            if (pChange < pr) { # test  proposal

              change<- 1
              transfers <- transfers + 1
              stratcy[samp,1]<- B    # update stratification
              Sd2[1,A] <-  Sd2[1,A]- sumd2tinA
              Sd2[1,B] <- Sd2[1,B] + sumd2plus
              cbObj <-  sqrt(Sd2)
              Obj <-  sum(cbObj)
              Eall <- rbind( Eall , Obj )
              Delta <- 0
            }


          } # end if srtat!=A

          #if(cycle == 1) browser()
          if (change == 1){break}

        } # End of For loop per stratum


      } # End of For loop per gp

      ###___________________________

      if(verbose==TRUE){print(paste(paste("Cycle Number=", cycle, sep=" "), paste("Transfers=", transfers, sep=" "), sep= "  :::---:::  "))}
      if(verbose==TRUE){print(paste('Temperature = ', Temp, sep=" "))}


      TotTransf<-  TotTransf + transfers
      # lower temperature
      Temp <- coolingRate * Temp

      if (transfers == 0) {break}
    } # End of cycles---------------------------------------------


    if(verbose==TRUE){print("-----END TRANSFER------")


      print(paste('Number of cycles= ', cycle, sep=" "))
      print(paste('Total number of transfers= ', TotTransf, sep=" "))}
    Obj<-  sum(cbObj)
    ObarFinal<- Obj/n
    print(paste('Objective function= ', ObarFinal, sep=" "))

    # Update if last ObarFinal is smaller than the smallest previous ObarFinal
    if (ObarFinal < ObarBest){
      ObarBest<- ObarFinal
      StratBest[,1]<- stratcy
    }
  }
  d1<- cbind(data,StratBest)
  names(d1)[ncol(d1)]<- "ospats_strata"
  retval<- list(objective = ObarBest, stratFrame = d1 , annealOut= Eall, covMat= d2, sumSq= Sd2, ospatsParams= params)
  return(retval)
}



