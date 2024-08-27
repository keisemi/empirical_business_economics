CCP_to_Value_to_prediction <- function(theta, CCP1, CCP2, TransitionMat, beta){
  
  # CCPから事前の価値をinvertし, その事前の価値からCCPを予測
  # 疑似データ作成時の手続きの一部
  # input: 
  #   theta: 10 x 1 parameter vector
  #   CCP1/2: 8 x 1 CCP essence vector of firm 1/2
  #   Transition Mat: 2 x 2 exogenous state transition matrix
  # output: 8 x 2 updated CCP
  
  # step1: setup
  # step2: CCPを所与としたときの, ValueのInversion
  # step3: ValueとパラメタをもとにCCPを予測
  
  # step1: setup----
  
  # CCP1Adjuster and CCP2Adjuster are matrices where available
  # choices are denoted by 1 and 0 otherwise
  CCP1Adjuster <- matrix(c(0, 1, 1, 
                           0, 1, 1, 
                           1, 1, 0,
                           1, 1, 0), ncol = 3, byrow = TRUE) # 4 x 3
  CCP1Adjuster <- rbind(CCP1Adjuster, CCP1Adjuster) # 8 x 3
  CCP2Adjuster <- matrix(c(0, 1, 1, 
                           1, 1, 0,
                           0, 1, 1,
                           1, 1, 0), ncol = 3, byrow = TRUE)
  CCP2Adjuster <- rbind(CCP2Adjuster, CCP2Adjuster)
  
  # profit pi1 and pi2
  pi1 <- pi1gen(theta) * CCP1Adjuster # Hadamard product
  pi2 <- pi2gen(theta) * CCP2Adjuster
  
  # transform CCP vectors into matrices
  CCP1Mat <- CCP1Transform(CCP1)
  CCP2Mat <- CCP2Transform(CCP2)
  
  # step2: CCPを所与としたときの, ValueのInversion----
  # F^P^sigmaを求める。これは遷移行列とCCP1とCCP2から作られる３つの行列の、
  # 要素ごとの積（アダマール積）
  fPsigma <- fP(TransitionMat, CCP1, CCP2)
  
  # パラメターの下で、pi1Psigma と pi2Psigma を計算する
  pi1Psigma <- pi1PsigmaGen(pi1, CCP2Mat)
  pi2Psigma <- pi2PsigmaGen(pi2, CCP1Mat)
  
  # e^P_i(a_i,s)を計算する
  # R automatically adjusts the dims of matrices:
  eP1 <- eulergamma - CCP1LogTransform(CCP1)
  eP2 <- eulergamma - CCP2LogTransform(CCP2)
  # lengthy alternative: 
  # eP1 <- eulergamma * matrix(rep(1, 8*3), ncol = 3) - CCP1LogTransform(CCP1)
  # eP2 <- eulergamma * matrix(rep(1, 8*3), ncol = 3) - CCP2LogTransform(CCP2)
  
  # 以上のパーツから事前の価値関数を求める
  ExanteV1 <- as.matrix(matlib::inv(diag(8) - beta*fPsigma) %*% 
                          ((CCP1Mat * (pi1Psigma+eP1)) %*% matrix(c(1,1,1))))
  ExanteV2 <- matlib::inv(diag(8) - beta*fPsigma) %*% 
    ((CCP2Mat * (pi2Psigma+eP2)) %*% matrix(c(1,1,1)))
  
  # Step3: 事前の価値関数を基に、CCPを計算する----
  
  # 各企業の行動が与えられたときの遷移行列を求める
  # リスト形式で１枚目から順にa1=-1, 0, 1のときの遷移行列
  fP_a1 <- fP_a1given(TransitionMat, CCP2)
  fP_a2 <- fP_a2given(TransitionMat, CCP1)
  
  # 以下では、求められた事前の価値関数をもとにCCPを更新
  NewSigmaSeed1  <- (pi1Psigma + beta*c(fP_a1[[1]]%*% ExanteV1, fP_a1[[2]]%*%ExanteV1 , fP_a1[[3]]%*%ExanteV1 ))*CCP1Adjuster
  NewSigmaDeno1  <- apply(exp(NewSigmaSeed1), 1, sum) - rep(1,8)
  NewSigmaDeno1M <- rep(NewSigmaDeno1,3)
  NewSigma1      <- exp(NewSigmaSeed1)/NewSigmaDeno1M
  CCP1UpdatedMat <- NewSigma1*CCP1Adjuster
  CCP1Updated    <- CCP1UpdatedMat[,2]
  
  NewSigmaSeed2  <-  (pi2Psigma + beta*c(fP_a2[[1]]%*% ExanteV2, fP_a2[[2]]%*%ExanteV2 , fP_a2[[3]]%*%ExanteV2 ))*CCP2Adjuster
  NewSigmaDeno2  <- apply(exp(NewSigmaSeed2), 1, sum) - rep(1,8)
  NewSigmaDeno2M <- rep(NewSigmaDeno2,3)
  NewSigma2      <- exp(NewSigmaSeed2)/NewSigmaDeno2M
  CCP2UpdatedMat <- NewSigma2*CCP2Adjuster
  CCP2Updated    <- CCP2UpdatedMat[,2]
  
  output <- cbind(CCP1Updated, CCP2Updated)
  
  output
  
}