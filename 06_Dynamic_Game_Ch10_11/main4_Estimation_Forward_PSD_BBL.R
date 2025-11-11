# Forward Simulationを用いたP-SD推定およびBBLの不等式推定量

# 本ファイルの流れは以下のようになっています
# 1. 下準備
# 2. 疑似データの生成 (第10回の内容に相当)
# 3. 推定 Step 1: CCPとTransition
# 4. 推定 Step 2-1: 推定したCCPの元でのForward Simulation
# 5. Forward Simulationを用いたPesendorer and Schmidt-Dengler (以下P-SD)
# 6. BBLの不等式推定量
# 7. BBLの不等式推定量におけるBootstrap

# 1. 下準備----

# データをクリア
rm(list = ls())

# メモリを開放
gc()

# ==== Prevent Mac from sleeping during long R runs ====

# Start caffeinate (prevents sleep while this R session runs)
# The "-dims" flags mean:
# -d : prevent display sleep
# -i : prevent idle sleep
# -m : prevent disk sleep
# -s : prevent system sleep
system("caffeinate -dims &")


# 紙面結果の再現のためのオプション
# option_matlab: FALSEの場合はRで乱数を生成＆データセット生成する。
# TRUEの場合はMatlabコード(書籍原稿を再現する分析)と同じ乱数＆データセットを用いる。
# なおTRUEにする場合、`data_from_matlab`フォルダにおいて以下のファイルを入れること。
# (サポートサイトからダウンロード可能)
# random_number_matlab_PSD.mat
# random_number_matlab_BBL.mat
# FakeData_Matlab.csv
# 留意点：紙面におけるForward simulationを用いたP-SDの推定値＆標準誤差、
# 及びBBL(不等式ベース)の推定値の結果を再現するようになっている。
# BBL(不等式ベース)の標準誤差の再現はできないので留意されたい。
# 詳細はGithubのReadmeを参照されたい。
isUseMatlabData <- TRUE

# Bootstrapにおける並列計算のオプション
# option_parallel: FALSEの場合は並列計算なし。TRUEの場合は並列計算を行う。
# 並列計算は（１）P-SD with forward simulationのBootstrap、（２）BBLのBootstrap で行う。
# なお、利用するコア数などは該当箇所で設定すること。(コードでは並列計算コア数を１０としている。)
option_parallel <- TRUE
option_numcores <- 10

# 必要なパッケージの読み込み
library(matlib) # 逆行列の作成に用いる
library(tictoc) # 実行時間を測るのに用いる
library(dplyr) # dplyr::if_else() に用いたがifelse()に替えることも可, see VsigmaGeneration
library(R.matlab) # matlab の.matファイルの読み込みに用いる。Matlabで用いたものと同じ乱数を利用するために必要。
library(foreach) # 並列計算のため
library(doParallel) # 並列計算のため
library(here) # ファイル読み込みのため

# 必要な関数の読み込み
functionlt <- list.files("06_Dynamic_Game_Ch10_11/functions_R",
  pattern = "*.R$", full.names = TRUE,
  ignore.case = TRUE
)
# sapply(functionlt[-20], source)
sapply(functionlt, source)

# 1. 下準備のコード。均衡計算を行う。
source(here("06_Dynamic_Game_Ch10_11/sub_1_prepare.R"))

# 2. 疑似データの作成 ----
# 疑似データ作成のコード。第11回と同じ。
source(here("06_Dynamic_Game_Ch10_11/sub_2_DGP.R"))

# Matlabの再現: Fakedata 読み込み
if (isUseMatlabData == 1) {
  FakeData <- read.csv(here("06_Dynamic_Game_Ch10_11/data_from_matlab/FakeData_Matlab.csv"), header = FALSE) %>% as.matrix()
}

# 3. 推定 Step 1: CCPとTransition----
# CCPと状態変数zの遷移確率を推定する。第11回と同じだが再掲。
EstimatedCCP1 <- matrix(0, 8, 1)
EstimatedCCP2 <- matrix(0, 8, 1)

for (s in 1:8) {
  SubData <- FakeData[FakeData[, 3] == s, ]
  Subseta1iszero <- SubData[SubData[, 7] == 0, ]
  Subseta2iszero <- SubData[SubData[, 8] == 0, ]
  EstimatedCCP1[s] <- nrow(Subseta1iszero) / nrow(SubData)
  EstimatedCCP2[s] <- nrow(Subseta2iszero) / nrow(SubData)
}

EstimatedTransition <- matrix(0, 2, 2)

DataLag1 <- matrix(c(0, FakeData[1:nrow(FakeData) - 1, 4])) # １期前のデータ
SubData <- cbind(FakeData, DataLag1)
SubData <- SubData[SubData[, 2] != 1, ] # t=1からt=2、t=2から・・・という遷移
for (z in 1:2) {
  SubDataZ <- SubData[SubData[, 4] == z, ]
  SubDataZnext <- SubDataZ[SubDataZ[, 9] == z, ]
  EstimatedTransition[z, z] <- nrow(SubDataZnext) / nrow(SubDataZ)
  EstimatedTransition[z, 3 - z] <- 1 - EstimatedTransition[z, z]
}

# 確認
print(cbind(EstimatedCCP1, EstimatedCCP2)) # 推定されたCCP
print(EstimatedTransition) # 推定されたZの遷移確率行列

# 4. 推定 Step 2-1: 推定したCCPの元でのForward Simulation----

# Note: State index and its corresponding values
# 1: G(1) 0 0
# 2: G(1) 0 1
# 3: G(1) 1 0
# 4: G(1) 1 1
# 5: B(2) 0 0
# 6: B(2) 0 1
# 7: B(2) 1 0
# 8: B(2) 1 1

# Forward simulationのパラメタ設定

NumSimPeriods <- 100
NumSimFirms <- 2
NumSimulations <- 1000

#  Forward simulationにおけるInitialのStateを設定。
# 今回は状態変数が取りうる値が８通りなので、その８通りをInitialとする。
InitialState <- array(1:8, dim = c(8, 2))
NumSimMarkets <- 8

# Forward simulationで用いるランダムショックを準備する。

# 乱数の発生が必要なため、乱数発生のためのシードを設定する
set.seed(2023)

# Idiosyncratic shock (ロジットショック)のドローを準備する。
# F\left(\varepsilon_{n j}\right)=e^{-e^{-\varepsilon_{n j}}}
# F(x) = exp(-exp(-x)) を使う。
# -log(-log(F(x)) ) = x となるので、F(x)の部分を0-1一様分布から引き、そのInversionを取る。
# 逆関数法により目的の確率分布を得られる
EVrandom <-
  array(
    -log(-log(runif(NumSimMarkets * NumSimPeriods * NumSimFirms * NumSimulations * 8 * 3))),
    dim = c(NumSimMarkets, NumSimPeriods, NumSimFirms, NumSimulations, 8, 3)
  )

# 景気のtransitionのシミュレーションに用いるドローを用意する。
UNIrandom <- array(runif(NumSimMarkets * NumSimPeriods * NumSimulations),
  dim = c(NumSimMarkets, NumSimPeriods, NumSimulations)
)

# Matlabの再現: 乱数 読み込み
if (isUseMatlabData == 1) {
  fa <- readMat(here("06_Dynamic_Game_Ch10_11/data_from_matlab/random_number_matlab_PSD.mat") )
  EVrandom <- fa$EVrandom
  UNIrandom <- fa$UNIrandom
}

# Sanity check
# DGPにおけるExante value functionと、Forward SimulationしたValueが一致するかをチェックする。

# DGPにおける CCP と Transition の元でForward Simulationする。
tic()
(Wstar <-
  VSigmaGeneration(
    CCP1UpdatedMat[, 2], CCP2UpdatedMat[, 2], TransitionMat, EVrandom,
    UNIrandom, InitialState, NumSimMarkets, NumSimulations,
    NumSimPeriods, beta
  ))
toc()
W1star <- Wstar[[1]]
W2star <- Wstar[[2]]

param1 <- matrix(c(TrueParameterValues[1:5], 1))
param2 <- matrix(c(TrueParameterValues[6:10], 1))

fa1 <- cbind(ExanteV1, t(W1star) %*% param1, ExanteV1Updated - t(W1star) %*% param1)
fa2 <- cbind(ExanteV2, t(W2star) %*% param2, ExanteV2Updated - t(W2star) %*% param2)
print("True value, Simulated value, and their difference")
fa1 # True, Simulated, and Diff
fa2 # True, Simulated, and Diff

# Check value after considering normalization a la Agg-Suzuki
Normalized_TrueParam <-
  matrix(c(
    Parameters[1] - ((1 - beta) / beta) * Parameters[5], # 企業1のベース利潤
    Parameters[3], # 企業2の店舗数の企業1への影響 ;顧客収奪効果
    Parameters[4], # 景気が良い時の企業1の追加的利潤
    0, # 企業１の退出のためのコスト; 0に標準化
    Parameters[6] + Parameters[5], # 企業1の参入のためのコスト
    Parameters[2] - ((1 - beta) / beta) * Parameters[5], # 企業2のベース利潤
    Parameters[3], # 企業1の店舗数の企業2への影響 ;顧客収奪効果
    Parameters[4], # 景気が良い時の企業2の追加的利潤
    0, # 企業2の退出のためのコスト; 0に標準化
    Parameters[6] + Parameters[5] # 企業2の参入のためのコスト
  ))
normparam1 <- matrix(c(Normalized_TrueParam[1:5], 1))
normparam2 <- matrix(c(Normalized_TrueParam[6:10], 1))

fafa1 <- cbind(ExanteV1, t(W1star) %*% param1, t(W1star) %*% normparam1)
fafa2 <- cbind(ExanteV2, t(W2star) %*% param2, t(W2star) %*% normparam2)
print("True value, simulated value (original parameter), simulated value (noramlized parameter)")
fafa1 # True, Simulated, and Simulated after normatlization
fafa2 # True, Simulated, and Simulated after normalization

print("Difference b.w. two simulated values (with and without normalization")
t(W1star) %*% param1 - t(W1star) %*% normparam1 # diff b/w bfr & aft normalization
t(W2star) %*% param2 - t(W2star) %*% normparam2 # diff b/w bfr & aft normalization

# 著者メモ：各企業が店舗を開いている(n=1)のStateにおいて、
# -0.1875 * 0.8 (discount factor) = - 0.15
# という差が生じている。
# これは、Closing costをゼロに基準化していることに起因する。

# 推定したCCPのもとで、Forward Simulationを行い、Value functionの基底を出す
Wstar <-
  VSigmaGeneration(
    EstimatedCCP1, EstimatedCCP2, EstimatedTransition, EVrandom,
    UNIrandom, InitialState, NumSimMarkets, NumSimulations,
    NumSimPeriods, beta
  )
W1star <- Wstar[[1]]
W2star <- Wstar[[2]]

# 5. P-SD with Forward Simulation ----

# Forward simulation したValueを用いて、P-SD流に推定を行う。

# まず目的関数の定義
obj_forward_PSD <-
  function(x) {
    Estimation_forward_PSD(
      c(x[1], x[3:4], 0, x[5], x[2], x[3:4], 0, x[5]),
      W1star, W2star, EstimatedTransition,
      EstimatedCCP1, EstimatedCCP2, beta
    )
  }

# 初期値はTrueにしておく。
initial <- c(0.3375, 0.2375, -0.27, 0.45, -2.25)

# 推定
res <- optim(par = initial, fn = obj_forward_PSD)
(opt_forwardPSD <- res[[1]])

# 以下のスクリプトを実行する。Bootstrapするので多少時間かかる。
# 並列計算推奨
tic()
source(here("06_Dynamic_Game_Ch10_11/sub_5_Bootstrap_PSD_forward.R"))
toc()

# 6. Estimation by BBL Inequality----

# 計算速度の観点から、以下ではSimulationの期間を30期間とする。
NumSimPeriods <- 30

# 乱数の発生が必要なため、乱数発生のためのシードを設定する
set.seed(2023)

# Idiosyncratic shock (ロジットショック)のドローを準備する。
EVrandom <-
  array(
    -log(-log(runif(NumSimMarkets * NumSimPeriods * NumSimFirms * NumSimulations * 8 * 3))),
    dim = c(NumSimMarkets, NumSimPeriods, NumSimFirms, NumSimulations, 8, 3)
  )

# 景気のtransitionのシミュレーションに用いるドローを用意する。
UNIrandom <- array(runif(NumSimMarkets * NumSimPeriods * NumSimulations),
  dim = c(NumSimMarkets, NumSimPeriods, NumSimulations)
)

# Matlabの再現: 乱数 読み込み
if (isUseMatlabData == 1) {
  # PerturbedCCP / UNIrandom / EVrandomを呼び出す
  fa <- readMat(here("06_Dynamic_Game_Ch10_11/data_from_matlab/random_number_matlab_BBL.mat") )
  UNIrandom <- fa$UNIrandom
  EVrandom <- fa$EVrandom
}

# 推定したCCPのもとで、Forward Simulationを行い、Value functionの基底を出す
Wstar <-
  VSigmaGeneration(
    EstimatedCCP1, EstimatedCCP2, EstimatedTransition, EVrandom,
    UNIrandom, InitialState, NumSimMarkets, NumSimulations,
    NumSimPeriods, beta
  )
W1star <- Wstar[[1]]
W2star <- Wstar[[2]]

# CCPのPertubationを行う。
NumPerturbations <- 200
PerturbedCCP1 <-
  rep(EstimatedCCP1, NumPerturbations) +
  matrix(rnorm(8 * NumPerturbations, mean = 0, sd = .1), nrow = 8)
PerturbedCCP2 <-
  rep(EstimatedCCP2, NumPerturbations) +
  matrix(rnorm(8 * NumPerturbations, mean = 0, sd = .1), nrow = 8)

# To make CCP be inside of [0,1] with some buffer
for (i in 1:8) {
  for (j in 1:NumPerturbations) {
    PerturbedCCP1[i, j] <- max(PerturbedCCP1[i, j], 0.001)
    PerturbedCCP1[i, j] <- min(PerturbedCCP1[i, j], 0.999)

    PerturbedCCP2[i, j] <- max(PerturbedCCP2[i, j], 0.001)
    PerturbedCCP2[i, j] <- min(PerturbedCCP2[i, j], 0.999)
  }
}

# Matlabの再現: 乱数 読み込み
if (isUseMatlabData == 1) {
  # PerturbedCCP を呼び出す。
  PerturbedCCP1 <- fa$PerturbedCCP1
  PerturbedCCP2 <- fa$PerturbedCCP2
}

W1_all <- array(rep(0, 6 * NumSimMarkets * NumPerturbations), dim = c(6, NumSimMarkets, NumPerturbations))
W2_all <- array(rep(0, 6 * NumSimMarkets * NumPerturbations), dim = c(6, NumSimMarkets, NumPerturbations))

# PertubationしたCCPを用いてForward simulationを行う
tic()
if (option_parallel == FALSE){
  for (per in 1:NumPerturbations) {
    cat(per, "")
    W1_p <- VSigmaGeneration(
      PerturbedCCP1[, per], EstimatedCCP2, EstimatedTransition,
      EVrandom, UNIrandom, InitialState, NumSimMarkets, NumSimulations, NumSimPeriods, beta
    )
    W1_all[, , per] <- W1_p[[1]]
    
    W2_p <- VSigmaGeneration(
      EstimatedCCP1, PerturbedCCP2[, per], EstimatedTransition,
      EVrandom, UNIrandom, InitialState, NumSimMarkets, NumSimulations, NumSimPeriods, beta
    )
    W2_all[, , per] <- W2_p[[2]]
  }
  
} else if (option_parallel == TRUE){
  
  # コア数の設定
  cores <- option_numcores
  cl <- makeCluster(cores, outfile = "")
  registerDoParallel(cl)
  
  # foreachで結果をリストとして返す
  results <- foreach( per = 1:NumPerturbations, .verbose = TRUE) %dopar% {
    
    cat(per, "")
    W1_p <- VSigmaGeneration(
      PerturbedCCP1[, per], EstimatedCCP2, EstimatedTransition,
      EVrandom, UNIrandom, InitialState, NumSimMarkets, NumSimulations, NumSimPeriods, beta
    )
    
    W2_p <- VSigmaGeneration(
      EstimatedCCP1, PerturbedCCP2[, per], EstimatedTransition,
      EVrandom, UNIrandom, InitialState, NumSimMarkets, NumSimulations, NumSimPeriods, beta
    )
    
    list(W1 = W1_p[[1]], W2 = W2_p[[2]])
    
    
  }
  stopCluster(cl)
  
  # foreachの結果（リスト）をまとめて配列に変換
  W1_all <- array(0, dim = c(6, NumSimMarkets, NumPerturbations))
  W2_all <- array(0, dim = c(6, NumSimMarkets, NumPerturbations))
  
  for (per in 1:NumPerturbations) {
    W1_all[, , per] <- results[[per]]$W1
    W2_all[, , per] <- results[[per]]$W2
  }
  
}
toc()

# Step 2-4: Minimizing the objective function

# Set initial
initial <- c(0.3, 0.2, -0.27, 0.45, -2.1)

# 目的関数の定義

# 補足：関数BBLobjective_NLSと関数BBLobjectiveは同じ目的関数の値を返す。
# ただし、BBLobjective_NLSの方が計算が早い。
# なお、matlabコードでは、関数BBLobjective_NLSは、二乗和を取る前の値をベクトルとして返すものであり、
# これを非線形最小二乗法の最適化アルゴリズムを用いて推定を行っていた。
# Rコードではこのステップを省いて、BBLobjective_NLSも通常の関数（値を返すもの）として定義し、
# optim関数で最適化を行っている。
obj_forward_BBL_NLS <-
  function(x) {
    BBLobjective_NLS(
      c(x[1], x[3:4], 0, x[5], x[2], x[3:4], 0, x[5]),
      NumPerturbations, W1star, W2star, W1_all, W2_all
    )
  }

opt <- optim(
  par = initial, fn = obj_forward_BBL_NLS,
  control = list(
    maxit = 1e4,
    abstol = 1e-10,
    reltol = 1e-10
  )
)

# 参考：同じ結果が得られる。
obj_forward_BBL <-
  function(x) {
    BBLobjective(
      c(x[1], x[3:4], 0, x[5], x[2], x[3:4], 0, x[5]),
      NumPerturbations, W1star, W2star, W1_all, W2_all
    )
  }
(opt <- optim(
  par = initial, fn = obj_forward_BBL,
  control = list(
    maxit = 1e4,
    abstol = 1e-10,
    reltol = 1e-10
  )
))

# 推定結果を保存する
saveRDS(opt, file = here("06_Dynamic_Game_Ch10_11/output/result_BBL_pointestimate.RDS") )

# 7. BBL (inequality estimator)によるBootstrap ----

# Pertubationの回数
NumPerturbations <- 200

# Bootstrapのための乱数を固定
set.seed(2023)

# Bootstrap のリサンプリング回数
numBootSample <- 100

# 各Bootstrap sampleで用いるマーケットのインデックスを乱数から発生させる。
# 市場が500個であるため、１から５００の整数について、重複を許して５００個ドローする。
bootindex <- array(sample(1:500, 500 * numBootSample, replace = TRUE), dim = c(500, numBootSample))

# Matlabの再現: 乱数 読み込み
if (isUseMatlabData == TRUE) {
  fa_PSD <- readMat(here("06_Dynamic_Game_Ch10_11/data_from_matlab/random_number_matlab_PSD.mat"))
  bootindex <- fa_PSD$bootindex
  rm(fa_PSD)
}

# 結果を保存するための行列
bootresult_payoff <- matrix(rep(0, 5 * numBootSample), ncol = numBootSample)

# 乱数の発生が必要なため、乱数発生のためのシードを設定する
set.seed(2023)
noise_CCP1 <- matrix(rnorm(8 * NumPerturbations, mean = 0, sd = .1), nrow = 8)
noise_CCP2 <- matrix(rnorm(8 * NumPerturbations, mean = 0, sd = .1), nrow = 8)

# Bootstrap を実行するループ

# 注意：Bootstrapには非常に長い時間がかかるため注意。
# なお、Bootstrapの結果については、output/result_BBL_Bootstrap/result_BBL_boot_XX.rds ファイルに格納されている。

if (option_parallel == FALSE) {
  # 並列計算なし

  for (b in 1:numBootSample) {
    # Bootstrap sampleを構築する。
    bootsample <- matrix(rep(0, 500 * 50 * 8), ncol = 8)
    for (m in 1:500) {
      temp <- FakeData[FakeData[, 1] == bootindex[m, b], ]
      bootsample[(1 + 50 * (m - 1)):(50 * m), ] <- temp
    }
    tic()

    # 構築したBootstrap sampleを用いて推定を行う。

    # 留意点：Matlabコードの方では、Bootstrap_BBL関数の内部で
    # EVrandom, UNIrandom、CCPのPertubationなどの乱数を発生させている。
    # 一方、Rコードの方では、必要となる乱数は関数の外で発生させており、引数として渡している。
    output <- Bootstrap_BBL(
      bootsample, beta,
      EVrandom, UNIrandom, InitialState,
      NumSimMarkets, NumSimulations, NumSimPeriods,
      NumPerturbations, noise_CCP1, noise_CCP2
    )

    # 推定結果を保存する。
    output_param <- output[[1]]$par
    bootresult_payoff[, b] <- output_param
    toc()
    cat(b, "")
  }
} else if (option_parallel == TRUE) {
  # コア数の設定
  cores <- option_numcores
  cl <- makeCluster(cores, outfile = "")
  registerDoParallel(cl)

  tic()
  foreach(b = 1:numBootSample, .verbose = TRUE) %dopar% {
    
    cat(b, "")
    # Bootstrap sampleを構築する。
    bootsample <- matrix(rep(0, 500 * 50 * 8), ncol = 8)
    for (m in 1:500) {
      temp <- FakeData[FakeData[, 1] == bootindex[m, b], ]
      bootsample[(1 + 50 * (m - 1)):(50 * m), ] <- temp
    }

    # 構築したBootstrap sampleを用いて推定を行う。

    # 留意点1：Matlabコードの方では、Bootstrap_BBL.m関数の内部で
    # EVrandom, UNIrandom、CCPのPertubationなどの乱数を発生させている。
    # 一方、Rコードの方では、必要となる乱数は関数の外で発生させており、引数として渡している。
    # 留意点2: Matlabコードで得られるBootstrapの結果を再現するには、
    # Bootstrap_BBL.m関数内部で生成される乱数を取得(外部に保存し、Rで利用)する必要がある。
    # その作業が若干困難なところがあるため、
    # 今回はRコードにおけるBBLのBoostrapにおいて、Matlabコードの結果を再現することを断念した。
    result <- Bootstrap_BBL(
      bootsample, beta,
      EVrandom, UNIrandom, InitialState,
      NumSimMarkets, NumSimulations, NumSimPeriods,
      NumPerturbations, noise_CCP1, noise_CCP2
    )

    # 推定結果を保存する。
    output <- list()
    output$param <- result[[1]]$par

    # 推定結果はRDSファイルとして保存している。
    saveRDS(output, file = paste("06_Dynamic_Game_Ch10_11/output/result_BBL_Bootstrap/result_BBL_boot_", b, ".RDS", sep = ""))
  }
  toc()

  for (b in 1:numBootSample) {
    output <- readRDS(paste("06_Dynamic_Game_Ch10_11/output/result_BBL_Bootstrap/result_BBL_boot_", b, ".RDS", sep = ""))
    bootresult_payoff[, b] <- output$param
  }
}

# 推定結果のまとめ
# output/result_BBL_boot_XXX.RDS と 
# output/result_BBL_Bootstrap/result_BBL_pointestimate.RDS に
# ブートストラップの結果と点推定値がそれぞれ保存されている。
true <- c(0.3, 0.2, -0.27, 0.45, -2.1)
print("Payoff parameter: True, Normalized true, Estimated, SE ")
normalized_param <- c(0.3 - (1 - beta) / beta * (-0.15), 0.2 - (1 - beta) / beta * (-0.15), -0.27, 0.45, -2.1 + (-0.15))
estmat_BBL <- matrix(c(true, normalized_param, opt$par, apply(bootresult_payoff, 1, sd)), nrow = 5)
# Save in csv
write.csv(estmat_BBL, file = here("06_Dynamic_Game_Ch10_11/output/Tab11_5_forward_BBL.csv"))




