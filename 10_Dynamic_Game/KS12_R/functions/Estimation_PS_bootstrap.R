Estimation_PS_bootstrap <- function(FakeData, beta) { # `{` added

  # P-SDの推定方法による、ブートストラップを行う。
  # FakeDataは、以下の推定に用いるデータ。ブートストラップにおいては、元のデータからリサンプリングしたデータになっている。

  # Step 1: Estimating CCP and transition probabilities
  EstimatedCCP1 <- matrix(rep(0, 8), ncol = 1)
  EstimatedCCP2 <- matrix(rep(0, 8), ncol = 1)

  for (s in 1:8) {
    SubData <- FakeData[FakeData[, 3] == s, ]
    Subseta1iszero <- SubData[SubData[, 7] == 0, ]
    Subseta2iszero <- SubData[SubData[, 8] == 0, ]
    EstimatedCCP1[s, 1] <- dim(Subseta1iszero)[1] / dim(SubData)[1]
    EstimatedCCP2[s, 1] <- dim(Subseta2iszero)[1] / dim(SubData)[1]
  }

  EstimatedTransition <- matrix(rep(0, 4), ncol = 2)

  DataL1 <- matrix(c(0, FakeData[1:dim(FakeData)[1] - 1, 4]), ncol = 1)
  SubData <- cbind(FakeData, DataL1)
  SubData <- SubData[SubData[, 2] != 1, ]

  for (z in 1:2) {
    SubDataZ <- SubData[SubData[, 4] == z, ]
    SubDataZnext <- SubDataZ[SubDataZ[, 9] == z, ]
    EstimatedTransition[z, z] <- dim(SubDataZnext)[1] / dim(SubDataZ)[1]
  }

  EstimatedTransition[1, 2] <- 1 - EstimatedTransition[1, 1]
  EstimatedTransition[2, 1] <- 1 - EstimatedTransition[2, 2]


  # Step 2: P-SD estimator
  obj <- function(x) {
    obj_fun(
      c(x[1], x[3:4], 0, x[5], x[2], x[3:4], 0, x[5]),
      EstimatedCCP1, EstimatedCCP2, EstimatedTransition,
      beta
    )
  }
  initial <- matrix(c(0.3, 0.2, -0.27, 0.45, -2.1))
  result <- optim(initial, obj)


  # Output
  output <- list(
    cbind(EstimatedCCP1, EstimatedCCP2),
    EstimatedTransition,
    result$par,
    result$value
  )

  output
}
