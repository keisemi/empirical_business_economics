## 5. P-SD + Forward SimulationにおけるBootstrap

# マーケット単位でリサンプリングを行う。マーケットは500個

# Bootstrapのための乱数を固定
set.seed(2023)

# Bootstrap のリサンプリング回数
numBootSample <- 100

# 各Bootstrap sampleで用いるマーケットのインデックスを乱数から発生させる。
# 市場が500個であるため、１から５００の整数について、重複を許して５００個ドローする。
bootindex <- matrix(sample(x = 1:500, size = 500 * numBootSample, replace = TRUE), ncol = numBootSample)

# Matlabの再現: 乱数 読み込み
if (option_matlab == TRUE) {
  fa <- readMat("data_from_matlab/random_number_matlab_PSD.mat")
  bootindex <- fa$bootindex
}

# 結果を保存するための行列
bootresult_transition <- matrix(rep(0, 2 * numBootSample), ncol = numBootSample)
bootresult_CCP1 <- matrix(rep(0, 8 * numBootSample), ncol = numBootSample)
bootresult_CCP2 <- matrix(rep(0, 8 * numBootSample), ncol = numBootSample)
bootresult_payoff <- matrix(rep(0, 5 * numBootSample), ncol = numBootSample)


# Bootstrap を実行するループ
# option_parallelで並列計算の有無を決める。

if (option_parallel == FALSE) {
  cat("Bootstrap: ")
  for (b in 1:numBootSample) {
    # print(append("Bootstrap:", as.character(b)))
    cat(b, "")
    # Bootstrap sampleを構築する。
    bootsample <- matrix(rep(0, 500 * 50 * 8), ncol = 8)
    for (m in 1:500) {
      temp <- FakeData[FakeData[, 1] == bootindex[m, b], ]
      bootsample[(1 + 50 * (m - 1)):(50 * m), ] <- temp
    }

    # 構築したBootstrap sampleを用いて推定を行う。
    output <- Bootstrap_PS_forward(
      bootsample, beta,
      EVrandom, UNIrandom, InitialState,
      NumSimMarkets, NumSimulations, NumSimPeriods
    )

    # 推定結果を保存する。
    bootresult_payoff[, b] <- output_param <- output[[1]]$par
    bootresult_CCP1[, b] <- output[[2]]
    bootresult_CCP2[, b] <- output[[3]]
    bootresult_transition[, b] <- c(output[[4]][1, 1], output[[4]][2, 2])
  }
} else if (option_parallel == TRUE) {
  # コア数の設定。
  cores <- 10
  cl <- makeCluster(cores, outfile = "")
  registerDoParallel(cl)

  tic()
  foreach(b = 1:numBootSample, .verbose = TRUE) %dopar% {
    # print(append("Bootstrap:", as.character(b)))
    cat(b, "")
    # Bootstrap sampleを構築する。
    bootsample <- matrix(rep(0, 500 * 50 * 8), ncol = 8)
    for (m in 1:500) {
      temp <- FakeData[FakeData[, 1] == bootindex[m, b], ]
      bootsample[(1 + 50 * (m - 1)):(50 * m), ] <- temp
    }

    # 構築したBootstrap sampleを用いて推定を行う。
    result <- Bootstrap_PS_forward(
      bootsample, beta,
      EVrandom, UNIrandom, InitialState,
      NumSimMarkets, NumSimulations, NumSimPeriods
    )

    # 推定結果を保存する。
    output <- list()
    output$payoff <- result[[1]]$par
    output$CCP1 <- result[[2]]
    output$CCP2 <- result[[3]]
    output$transition <- c(result[[4]][1, 1], result[[4]][2, 2])

    # 推定結果はRDSファイルとして保存している。
    saveRDS(output, file = paste("result/result_forwardPSD_boot_", b, ".RDS", sep = ""))
  }
  toc()

  stopCluster(cl)

  # 結果のSummary
  for (b in 1:numBootSample) {
    output <- readRDS(paste("result/result_forwardPSD_boot_", b, ".RDS", sep = ""))
    bootresult_payoff[, b] <- output$payoff
    bootresult_CCP1[, b] <- output$CCP1
    bootresult_CCP2[, b] <- output$CCP2
    bootresult_transition[, b] <- output$transition
  }
}




# Payoff parameter
true <- c(0.3, 0.2, -0.27, 0.45, -2.1)
print("Payoff parameter: True, Normalized true, Estimated, SE ")
normalized_param <- c(0.3 - (1 - beta) / beta * (-0.15), 0.2 - (1 - beta) / beta * (-0.15), -0.27, 0.45, -2.1 + (-0.15))
print(matrix(c(true, normalized_param, opt, apply(bootresult_payoff, 1, sd)), nrow = 5))
