Bootstrap_BBL <- function(FakeData, beta, 
                          EVrandom, UNIrandom, InitialState, 
                          NumSimMarkets, NumSimulations, NumSimPeriods,
                          NumPerturbations,noise_CCP1, noise_CCP2){
  
  # step 1: Estimating CCP and transition Probabilities----
  EstimatedCCP1 <- matrix(0, 8, 1)
  EstimatedCCP2 <- matrix(0, 8, 1)
  
  for (s in 1:8) {
    SubData <- FakeData[FakeData[,3]==s,]
    Subseta1iszero <- SubData[SubData[,7]==0, ]
    Subseta2iszero <- SubData[SubData[,8]==0, ]
    EstimatedCCP1[s] <- nrow(Subseta1iszero)/nrow(SubData)
    EstimatedCCP2[s] <- nrow(Subseta2iszero)/nrow(SubData)
  }
  
  EstimatedTransition <- matrix(0, 2, 2)
  
  DataLag1 <- matrix(c(0, FakeData[1:nrow(FakeData)-1, 4])) # １期前のデータ
  SubData <- cbind(FakeData, DataLag1)
  SubData <- SubData[SubData[,2]!=1,] # t=1からt=2、t=2から・・・という遷移
  for (z in 1:2) {
    SubDataZ <- SubData[SubData[,4]==z,]
    SubDataZnext <- SubDataZ[SubDataZ[,9]==z,]
    EstimatedTransition[z,z] <- nrow(SubDataZnext)/nrow(SubDataZ)
    EstimatedTransition[z,3-z] <- 1 - EstimatedTransition[z,z] 
  }
  
  # Step 2: Forward simulation using estimated CCP ----
    
  # Forward simulationのパラメタ設定
  NumSimFirms    <- 2

  ## Step 2-1: Generating V for any sigma----
  # Given EstimatedCCPs (and parameter values), calculate a
  # discounted sum of profit for each firm and basis functions
  
  # True CCP and Transition
  Wstar <- 
    VSigmaGeneration(EstimatedCCP1, EstimatedCCP2, EstimatedTransition, EVrandom,
                     UNIrandom, InitialState, NumSimMarkets, NumSimulations,
                     NumSimPeriods, beta)
  W1star <- Wstar[[1]]
  W2star <- Wstar[[2]]
  
  ## Step 2-2: Forward simulation for perturbed CCP ----
  
  PerturbedCCP1 <- rep(EstimatedCCP1, NumPerturbations) + noise_CCP1
  PerturbedCCP2 <- rep(EstimatedCCP2, NumPerturbations) + noise_CCP2

  # To make CCP be inside of [0,1] with some buffer
  for (i in 1:8){
    for (j in 1:NumPerturbations){
      PerturbedCCP1[i,j] <- max(PerturbedCCP1[i,j], 0.001)
      PerturbedCCP1[i,j] <- min(PerturbedCCP1[i,j], 0.999)
      
      PerturbedCCP2[i,j] <- max(PerturbedCCP2[i,j], 0.001)
      PerturbedCCP2[i,j] <- min(PerturbedCCP2[i,j], 0.999)
    }
  }
  
  W1_all <- array(rep(0,6*NumSimMarkets*NumPerturbations),dim = c(6,NumSimMarkets,NumPerturbations))
  W2_all <- array(rep(0,6*NumSimMarkets*NumPerturbations),dim = c(6,NumSimMarkets,NumPerturbations))
  
  # tic()
  for (per in 1:NumPerturbations){
    W1_p  <- VSigmaGeneration(PerturbedCCP1[,per], EstimatedCCP2, EstimatedTransition,
                              EVrandom, UNIrandom, InitialState, NumSimMarkets, NumSimulations, NumSimPeriods, beta)
    W1_all[ , ,per] <- W1_p[[1]]
    
    W2_p  <- VSigmaGeneration(EstimatedCCP1, PerturbedCCP2[ ,per], EstimatedTransition,
                              EVrandom, UNIrandom, InitialState, NumSimMarkets, NumSimulations, NumSimPeriods, beta)
    W2_all[ , ,per] <- W2_p[[2]]
  }
  # toc()
  
  # Step 3: Mimimizing the objective function----
  
  # Set initial
  initial <- c(0.3, 0.2, -0.27, 0.45, -2.1)
  
  # 推定：Nonlinear least squaresとして推定する。
  obj_forward_BBL_NLS <- 
    function(x){BBLobjective_NLS(c(x[1], x[3:4], 0, x[5], x[2],x[3:4], 0, x[5]),NumPerturbations
                                 , W1star, W2star, W1_all, W2_all)}
  
  result <-  optim(par = initial, fn = obj_forward_BBL_NLS,
                   method = "L-BFGS-B",
                   lower = c(0, 0, -Inf, 0, -Inf),
                   upper = c(Inf, Inf, 0, Inf, 0),
                   control = list(maxit = 1e5,
                                  factr = 1e-10))
  
  output <- list(result,EstimatedCCP1,EstimatedCCP2,EstimatedTransition)
}