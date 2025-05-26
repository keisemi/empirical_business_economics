VSigmaGeneration <-
  function(CCP1, CCP2, EstimatedTransition, EVrandom, UNIrandom, InitialState,
           NumSimMarkets, NumSimulations, NumSimPeriods, beta) {
    # # for sanity check
    # CCP1 <- CCP1UpdatedMat[,2]
    # CCP2 <- CCP2UpdatedMat[,2]
    # EstimatedTransition <- TransitionMat

    # preparation
    # here is the core part of the trick to reduce computational burden as in BBL;
    # these will be the coefficients of the structural parameters (in each market and simulation)
    # that express VSigma of both firms:
    # 1st component is the base profit for firm i(=1,2), 2-5th component is the
    # firm-common parameters, and the sixth component is the profit from random
    W1 <- array(0, dim = c(6, NumSimMarkets, NumSimulations))
    W2 <- array(0, dim = c(6, NumSimMarkets, NumSimulations))

    for (mrkt in 1:NumSimMarkets) {
      # mrkt <- 1

      # threshold values should not be varied by simulations
      # these are 8*1 vectors, one for each firm
      ThresholdValue1 <- log(CCP1) - log(1 - CCP1)
      ThresholdValue2 <- log(CCP2) - log(1 - CCP2)

      for (sim in 1:NumSimulations) {
        # sim <- 1

        # extracting the data
        State <- InitialState[mrkt, 2]

        for (t in 1:NumSimPeriods) {
          # cat("period", t, "state", State)

          DiffEpsilon1inc <- -EVrandom[mrkt, t, 1, sim, State, 2] + EVrandom[mrkt, t, 1, sim, State, 3]
          DiffEpsilon1dec <- -EVrandom[mrkt, t, 1, sim, State, 2] + EVrandom[mrkt, t, 1, sim, State, 1]

          DiffEpsilon2inc <- -EVrandom[mrkt, t, 2, sim, State, 2] + EVrandom[mrkt, t, 2, sim, State, 3]
          DiffEpsilon2dec <- -EVrandom[mrkt, t, 2, sim, State, 2] + EVrandom[mrkt, t, 2, sim, State, 1]

          DiffFirm1inc <- ThresholdValue1[State] - DiffEpsilon1inc
          DiffFirm1dec <- ThresholdValue1[State] - DiffEpsilon1dec

          DiffFirm2inc <- ThresholdValue2[State] - DiffEpsilon2inc
          DiffFirm2dec <- ThresholdValue2[State] - DiffEpsilon2dec

          # determine firm 1's action, i.e., entry/exit decision for next period
          # a1 = 0 means stay, otherwise (a1 =1 ) means increasing or decreasing
          # e1 is the profit from random component
          if (State %in% c(1, 2, 5, 6)) {
            # n1 = 0 under this condition
            a1 <- 0 * (DiffFirm1inc > 0) + 1 * (DiffFirm1inc < 0)
            e1 <- (1 - a1) * EVrandom[mrkt, t, 1, sim, State, 2] + a1 * EVrandom[mrkt, t, 1, sim, State, 3]
          } else {
            # n1 =1
            a1 <- 0 * (DiffFirm1dec > 0) + 1 * (DiffFirm1dec < 0)
            e1 <- (1 - a1) * EVrandom[mrkt, t, 1, sim, State, 2] + a1 * EVrandom[mrkt, t, 1, sim, State, 1]
          }

          # firm 2's action
          if (State %in% c(1, 3, 5, 7)) {
            # n2 =0
            a2 <- 0 * (DiffFirm2inc > 0) + 1 * (DiffFirm2inc < 0)
            e2 <- (1 - a2) * EVrandom[mrkt, t, 2, sim, State, 2] + a2 * EVrandom[mrkt, t, 2, sim, State, 3]
          } else {
            # n1 =1
            a2 <- 0 * (DiffFirm2dec > 0) + 1 * (DiffFirm2dec < 0)
            e2 <- (1 - a2) * EVrandom[mrkt, t, 2, sim, State, 2] + a2 * EVrandom[mrkt, t, 2, sim, State, 1]
          }

          # Finally, decide the next state based on the current
          # state, a1, a2, and transition matrix

          if (State == 1) {
            W1Seed <- matrix(c(0, 0, 0, 0, a1, e1))
            W2Seed <- matrix(c(0, 0, 0, 0, a2, e2))

            # if(a1==0 & a2==0){
            #   NextState <- 1
            # } else if(a1==0 & a2==1){
            #   NextState <- 2
            # } else if(a1==1 & a2==0){
            #   NextState <- 3
            # } else{
            #   NextState <- 4
            # }
            NextState <- 1 + 2 * a1 + a2

            # 一様乱数がP(G|G)の値を超えたら、景気がTransitionするとみなす。
            # （つまり, 来季にBへ移動する。）
            # TransitionSeedが1になる＝来季がBになる。
            TransitionSeed <- (UNIrandom[mrkt, t, sim] > EstimatedTransition[1, 1])
            # 1から4はG、5から8はBに対応している。
            NextState <- 4 * TransitionSeed + NextState
          } else if (State == 2) {
            W1Seed <- matrix(c(0, 0, 0, 0, a1, e1))
            W2Seed <- matrix(c(1, 0, 1, a2, 0, e2))

            # if(a1==0 & a2==0){
            #   NextState <- 2
            # } else if(a1==0 & a2==1){
            #   NextState <- 1
            # } else if(a1==1 & a2==0){
            #   NextState <- 4
            # } else{
            #   NextState <- 3
            # }

            NextState <- 2 + 2 * a1 - a2
            TransitionSeed <- (UNIrandom[mrkt, t, sim] > EstimatedTransition[1, 1])
            NextState <- 4 * TransitionSeed + NextState
          } else if (State == 3) {
            W1Seed <- matrix(c(1, 0, 1, a1, 0, e1))
            W2Seed <- matrix(c(0, 0, 0, 0, a2, e2))

            # if(a1==0 & a2==0){
            #   NextState <- 3
            # } else if(a1==0 & a2==1){
            #   NextState <- 4
            # } else if(a1==1 & a2==0){
            #   NextState <- 1
            # } else{
            #   NextState <- 2
            # }

            NextState <- 3 - 2 * a1 + a2
            TransitionSeed <- (UNIrandom[mrkt, t, sim] > EstimatedTransition[1, 1])
            NextState <- 4 * TransitionSeed + NextState
          } else if (State == 4) {
            W1Seed <- matrix(c(1, 1, 1, a1, 0, e1))
            W2Seed <- matrix(c(1, 1, 1, a2, 0, e2))

            # if(a1==0 & a2==0){
            #   NextState <- 4
            # } else if(a1==0 & a2==1){
            #   NextState <- 3
            # } else if(a1==1 & a2==0){
            #   NextState <- 2
            # } else{
            #   NextState <- 1
            # }

            NextState <- 4 - 2 * a1 - a2
            TransitionSeed <- (UNIrandom[mrkt, t, sim] > EstimatedTransition[1, 1])
            NextState <- 4 * TransitionSeed + NextState
          } else if (State == 5) {
            W1Seed <- matrix(c(0, 0, 0, 0, a1, e1))
            W2Seed <- matrix(c(0, 0, 0, 0, a2, e2))

            # if(a1==0 & a2==0){
            #   NextState <- 1
            # } else if(a1==0 & a2==1){
            #   NextState <- 2
            # } else if(a1==1 & a2==0){
            #   NextState <- 3
            # } else{
            #   NextState <- 4
            # }

            NextState <- 1 + 2 * a1 + a2

            # 一様乱数がP(B|B)の値よりも低かったら、景気がTransitionしない。
            # つまり、来季はのまま。
            # よって、TransitionSeedが1になる＝来季がBになる。
            TransitionSeed <- (UNIrandom[mrkt, t, sim] < EstimatedTransition[2, 2])
            NextState <- 4 * TransitionSeed + NextState
          } else if (State == 6) {
            W1Seed <- matrix(c(0, 0, 0, 0, a1, e1))
            W2Seed <- matrix(c(1, 0, 0, a2, 0, e2))

            # if(a1==0 & a2==0){
            #   NextState <- 2
            # } else if(a1==0 & a2==1){
            #   NextState <- 1
            # } else if(a1==1 & a2==0){
            #   NextState <- 4
            # } else{
            #   NextState <- 3
            # }
            NextState <- 2 + 2 * a1 - a2
            TransitionSeed <- (UNIrandom[mrkt, t, sim] < EstimatedTransition[2, 2])
            NextState <- 4 * TransitionSeed + NextState
          } else if (State == 7) {
            W1Seed <- matrix(c(1, 0, 0, a1, 0, e1))
            W2Seed <- matrix(c(0, 0, 0, 0, a2, e2))

            # if(a1==0 & a2==0){
            #   NextState <- 3
            # } else if(a1==0 & a2==1){
            #   NextState <- 4
            # } else if(a1==1 & a2==0){
            #   NextState <- 1
            # } else{
            #   NextState <- 2
            # }

            NextState <- 3 - 2 * a1 + a2
            TransitionSeed <- (UNIrandom[mrkt, t, sim] < EstimatedTransition[2, 2])
            NextState <- 4 * TransitionSeed + NextState
          } else {
            W1Seed <- matrix(c(1, 1, 0, a1, 0, e1))
            W2Seed <- matrix(c(1, 1, 0, a2, 0, e2))

            # if(a1==0 & a2==0){
            #   NextState <- 4
            # } else if(a1==0 & a2==1){
            #   NextState <- 3
            # } else if(a1==1 & a2==0){
            #   NextState <- 2
            # } else{
            #   NextState <- 1
            # }

            NextState <- 4 - 2 * a1 - a2
            TransitionSeed <- (UNIrandom[mrkt, t, sim] < EstimatedTransition[2, 2])
            NextState <- 4 * TransitionSeed + NextState
          }

          W1[, mrkt, sim] <- W1[, mrkt, sim] + beta^(t - 1) * W1Seed
          W2[, mrkt, sim] <- W2[, mrkt, sim] + beta^(t - 1) * W2Seed

          State <- NextState
        }
      }
    }

    W1out <- apply(W1, c(1, 2), mean)
    W2out <- apply(W2, c(1, 2), mean)

    output <- list(W1out, W2out)
    output
  }
