Bootstrap_PS_forward <- function(FakeData, beta,
                                 EVrandom, UNIrandom, InitialState, 
                                 NumSimMarkets, NumSimulations, NumSimPeriods){

  # conduct forward simulation from generated data (& beta)
  # output: the result of optimization, estimated CCPs, and estimated transition probability
  
  # Step 1: Estimating CCP and transition probabilities
  EstimatedCCP1 <- matrix(rep(0,8),ncol = 1)
  EstimatedCCP2 <-  matrix(rep(0,8),ncol = 1)
  
  for( s in 1:8){
    SubData <- FakeData[FakeData[ ,3] == s, ]
    Subseta1iszero <- SubData[SubData[ ,7]==0, ]
    Subseta2iszero <- SubData[SubData[ ,8]==0, ]
    EstimatedCCP1[s,1] <- dim(Subseta1iszero)[1]/dim(SubData)[1]
    EstimatedCCP2[s,1] <- dim(Subseta2iszero)[1]/dim(SubData)[1]
  }
    
  EstimatedTransition <- matrix(rep(0,4),ncol = 2)
  DataL1 <- matrix( c(0,FakeData[1:dim(FakeData)[1] - 1 , 4 ]) ,ncol = 1)
  SubData <- cbind(FakeData,DataL1)
  SubData <- SubData[SubData[ ,2] !=1, ]
    
  for (z in 1:2){
    SubDataZ     <- SubData[SubData[ ,4] == z, ]
    SubDataZnext <- SubDataZ[SubDataZ[ ,9] == z, ]
    EstimatedTransition[z,z] <- dim(SubDataZnext)[1]/dim(SubDataZ)[1]
  }
    
  EstimatedTransition[1,2] <- 1-EstimatedTransition[1,1]
  EstimatedTransition[2,1] <- 1-EstimatedTransition[2,2]
    
  # Step 2: P-SD estimator with forward simulation
  NumSimFirms    <- 2
  
  #NumSimPeriods  <- 100
  #NumSimulations <- 1000
  #InitialState <- cbind(seq(1,8,1),seq(1,8,1))
  #NumSimMarkets <- 8
  
  # Step 2-1: Generating epsilons from extreme value distributions for
  # choices and epsilons from uniform distribution for state
  # transition
  
  # EVrandom <- 
  #   array(
  #     -log(-log(runif(NumSimMarkets * NumSimPeriods * NumSimFirms * NumSimulations * 8 * 3))),
  #     dim = c(NumSimMarkets, NumSimPeriods, NumSimFirms, NumSimulations, 8, 3))
  # 
  # UNIrandom <- array(runif(NumSimMarkets * NumSimPeriods * NumSimulations),
  #                    dim=c(NumSimMarkets,NumSimPeriods,NumSimulations))
  
  # Step 2-2: Generating V for any sigma
  # Given EstimatedCCPs (and parameter values), calculate a
  # discounted sum of profit for each firm and basis functions
  
  # True CCP and Transition
  output <- VSigmaGeneration(EstimatedCCP1, EstimatedCCP2, EstimatedTransition,
    EVrandom, UNIrandom, InitialState, NumSimMarkets, NumSimulations, NumSimPeriods, beta)
  W1star <- output[[1]]
  W2star <- output[[2]]
  
  obj_forward_PSD <- 
    function(params){Estimation_forward_PSD(c(params[1], params[3:4], 0, params[5], 
                                              params[2], params[3:4], 0, params[5]),
                                            W1star, W2star, EstimatedTransition, 
                                            EstimatedCCP1, EstimatedCCP2, beta)}
  
  initial <- c(0.3, 0.2, -0.27, 0.45, -2.1)
                           
  result <-  optim(initial,obj_forward_PSD)
  output <- list(result,EstimatedCCP1,EstimatedCCP2,EstimatedTransition)
  
  output

}
