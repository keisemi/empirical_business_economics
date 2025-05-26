Estimation_forward_PSD <- function(param, W1star, W2star, TransitionMat, CCP1, CCP2, beta) {
  NumMaxChoices <- 3

  # Forward simuからValueを作る。
  param1 <- c(param[1:5], 1)
  param2 <- c(param[6:10], 1)

  ExanteV1 <- t(W1star) %*% param1
  ExanteV2 <- t(W2star) %*% param2

  # CCP1Adjuster and CCP2Adjuster are matrices where available
  # choices are denoted by 1 and 0 otherwise
  CCP1Adjuster <- matrix(c(
    0, 1, 1,
    0, 1, 1,
    1, 1, 0,
    1, 1, 0
  ), ncol = 3, byrow = TRUE) # 4 x 3
  CCP1Adjuster <- rbind(CCP1Adjuster, CCP1Adjuster) # 8 x 3
  CCP2Adjuster <- matrix(c(
    0, 1, 1,
    1, 1, 0,
    0, 1, 1,
    1, 1, 0
  ), ncol = 3, byrow = TRUE)
  CCP2Adjuster <- rbind(CCP2Adjuster, CCP2Adjuster)


  # Profit pi1 and pi2
  pi1 <- pi1gen(param) * CCP1Adjuster
  pi2 <- pi2gen(param) * CCP2Adjuster


  # Transforming CCP vectors into matrices
  CCP1Mat <- CCP1Transform(CCP1)
  CCP2Mat <- CCP2Transform(CCP2)


  # Given parameter values, calculate pi1Psigma and pi2Psigma
  pi1Psigma <- pi1PsigmaGen(pi1, CCP2Mat)
  pi2Psigma <- pi2PsigmaGen(pi2, CCP1Mat)



  fP_a1 <- fP_a1given(TransitionMat, CCP2)
  fP_a2 <- fP_a2given(TransitionMat, CCP1)

  # Given ExanteVF1 and ExanteVF2, calculate CCPs again
  NewSigmaSeed1 <- (pi1Psigma + beta * c(fP_a1[[1]] %*% ExanteV1, fP_a1[[2]] %*% ExanteV1, fP_a1[[3]] %*% ExanteV1)) * CCP1Adjuster
  NewSigmaDeno1 <- apply(exp(NewSigmaSeed1), 1, sum) - rep(1, 8)
  NewSigmaDeno1M <- rep(NewSigmaDeno1, 3)
  NewSigma1 <- exp(NewSigmaSeed1) / NewSigmaDeno1M
  CCP1UpdatedMat <- NewSigma1 * CCP1Adjuster
  CCP1Updated <- CCP1UpdatedMat[, 2]

  NewSigmaSeed2 <- (pi2Psigma + beta * c(fP_a2[[1]] %*% ExanteV2, fP_a2[[2]] %*% ExanteV2, fP_a2[[3]] %*% ExanteV2)) * CCP2Adjuster
  NewSigmaDeno2 <- apply(exp(NewSigmaSeed2), 1, sum) - rep(1, 8)
  NewSigmaDeno2M <- rep(NewSigmaDeno2, 3)
  NewSigma2 <- exp(NewSigmaSeed2) / NewSigmaDeno2M
  CCP2UpdatedMat <- NewSigma2 * CCP2Adjuster
  CCP2Updated <- CCP2UpdatedMat[, 2]

  # Output
  CCP_updated <- matrix(c(CCP1Updated, CCP2Updated), ncol = 2)


  # Value function

  obj <- sum((CCP_updated[, 1] - CCP1)^2) + sum((CCP_updated[, 2] - CCP2)^2)
}
