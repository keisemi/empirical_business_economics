# Aguirregabiria and Mira (2007) の推定方法

# 1. 下準備----

# データをクリア
rm(list = ls())

# メモリを開放
gc()

# 必要なパッケージの読み込み

# pacman パッケージを使って読み込みを行う。要インストール
# if (!require("pacman")) install.packages("pacman")

pacman::p_load(
  # 必要となるパッケージ名
  matlib, # 逆行列の作成に用いる
  R.matlab, # matlab の.matファイルの読み込みに用いる。Matlabで用いたものと同じ乱数を利用するために必要。
  pracma, # matlabと同じ原理での最適化に用いる
  here
)

# 必要な関数の読み込み
functionlt <- list.files("06_Dynamic_Game_Ch10_11/functions_R/",
  pattern = "*.R$", full.names = TRUE,
  ignore.case = TRUE
)
sapply(functionlt, source)

# 2. 均衡計算を行う
# "main1_Computation_Equilibrium.R"の前半部分が"sub_1_prepare.R"に入っている。
source("06_Dynamic_Game_Ch10_11/sub_1_prepare.R", echo = TRUE)

# 3. 疑似データの生成----
# "main1_Computation_Equilibrium.R"の後半部分が"sub_2_DGP.R"に入っている。
source("06_Dynamic_Game_Ch10_11/sub_2_DGP.R", echo = TRUE)

# 4. Aguirregabiria and Mira (2007)の方法によるパラメターの推定----

# 紙面における点推定値を同じものを得るためには、
# Matlabで生成したFakedataを利用する必要がある
# なおこの場合でも、Bootstrapにより計算した標準誤差については乱数の影響で若干の差が出ることには留意されたい

# isUseMatlabData = 0のときは上 "06_Dynamic_Game_Ch10_11/sub_2_DGP.R" で生成したFakedataを利用する
# isUseMatlabData = 1のときは、Matlabで生成したFakedata（FakeData_Matlab.csvファイル）を利用する
isUseMatlabData <- 1

if (isUseMatlabData == 1) {
  # Matlab で作ったFake dataをロードする
  dt <- read.csv(here("06_Dynamic_Game_Ch10_11/data_from_matlab/FakeData_Matlab.csv"), header = FALSE)
  FakeData <- as.matrix(dt)
} 

# Step 1: データから推定される遷移確率行列やCCPの導出
# データから推定されるCCP
EstimatedCCP1 <- matrix(0, 8, 1)
EstimatedCCP2 <- matrix(0, 8, 1)


for (s in 1:8) {
  SubData <- FakeData[FakeData[, 3] == s, ]
  Subseta1iszero <- SubData[SubData[, 7] == 0, ]
  Subseta2iszero <- SubData[SubData[, 8] == 0, ]
  EstimatedCCP1[s] <- nrow(Subseta1iszero) / nrow(SubData)
  EstimatedCCP2[s] <- nrow(Subseta2iszero) / nrow(SubData)
}

# データから推定される遷移確率行列
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

# Step 2: 初期値の設定
# パラメターの初期値を設定する
# この初期値は、初期値のPに対して最適化された後は使わない
# InitialParameters1は真のパラメターに近い初期値
InitialParameters1 <- matrix(c(
  0.3, # 企業1のベース利潤
  0.3, # 企業2のベース利潤
  -0.25, # 顧客収奪効果
  0.4, # 景気が良い時の追加的利潤
  -0.12, # 退出のためのコスト
  -2.2 # 参入のためのコスト
))
# InitialParameters2は真のパラメターから遠い初期値
InitialParameters2 <- matrix(c(
  0.1, # 企業1のベース利潤
  0.1, # 企業2のベース利潤
  -0.1, # 顧客収奪効果
  0.1, # 景気が良い時の追加的利潤
  -0.1, # 退出のためのコスト
  -1.0 # 参入のためのコスト
))
# ここではInitialParameters1を使う
InitialParameters <- InitialParameters1

InitialParameterValues <- matrix(c(
  InitialParameters[1], # 企業1のベース利潤
  InitialParameters[3], # ライバルの店舗数が企業1の利潤に与える影響
  InitialParameters[4], # 景気が良い時の企業1への追加的利潤
  InitialParameters[5], # 企業1の退出のためのコスト
  InitialParameters[6], # 企業1の出店のためのコスト
  InitialParameters[2], # 企業2のベース利潤
  InitialParameters[3], # ライバルの店舗数が企業2の利潤に与える影響
  InitialParameters[4], # 景気が良い時の企業2への追加的利潤
  InitialParameters[5], # 企業2の退出のためのコスト
  InitialParameters[6] # 企業2の出店のためのコスト
))

# CCPの初期値．以後はループの中で更新される
# FakeDataから推定されるCCPそのものを使う
ccp1 <- EstimatedCCP1
ccp2 <- EstimatedCCP2

# Step 3: Aguirregabiria and Mira (2007)の方法によるパラメターの推定

# CCPの更新回数
i <- 1

# データから観察された行動
Actions <- FakeData[, c(7:8, 3)]

# パラメターの初期値
initial <- c(InitialParameters[1:4], InitialParameters[6])

# 推定
while (i > 0) {
  # 時間がかかりすぎる場合に強制的に停止する
  if (i == 10000) {
    print("CCP did not converge")
    break
  }
  # 何回目のループか表示
  paste0("第", i, "回目のループ") |> print()

  # 目的関数objの定義
  # 対数尤度関数obj_likがパラメタ、CCP、遷移確率、割引因子に依存するのに対し、
  # 目的関数objはパラメタのみに依存する
  # CCPと遷移確率は前回のループ（または初期値）で得られたものを与える
  obj <- function(x) {
    obj_lik(
      c(x[1], x[3:4], 0, x[5], x[2], x[3:4], 0, x[5]),
      ccp1, ccp2,
      EstimatedTransition, beta, Actions
    )
  }

  # CCPを所与としたパラメターの最適化
  # 疑似尤度関数を最大化するようなパラメターを求める
  sol <- optim(initial, obj, control = list(fnscale = -1))

  # 尤度関数を最大化するようなパラメターのもとでCCPを計算する
  param <- c(sol$par[1], sol$par[3:4], 0, sol$par[5], sol$par[2], sol$par[3:4], 0, sol$par[5])
  newCCP <- CCP_to_Value_to_prediction(param, ccp1, ccp2, EstimatedTransition, beta)

  # 結果を表示
  print("尤度")
  print(sol$value)
  print("パラメター")
  print(param)
  print("CCP")
  print(cbind(ccp1, ccp2))

  # CCPが収束していないか確認する
  print("収束の判定")
  if (check_convergence(newCCP, cbind(ccp1, ccp2), tol = 1e-6) ) break

  # CCPが収束していない場合、新しいCCPでもとのCCPを更新する
  ccp1 <- newCCP[, 1]
  ccp2 <- newCCP[, 2]

  # iを更新
  i <- i + 1
}


# 推定されたパラメター
print("推定されたパラメター")
print(sol$par)
print("真のパラメター")
print(Parameters)

print("AM2007のアルゴリズムの収束地点におけるCCP")
print(newCCP)
print("真のCCP")
print(cbind(CCP1Updated, CCP2Updated))

# 5. Bootstrapによる標準誤差の計算 ----
# マーケット単位でリサンプリングを行う。マーケットは500個


# Bootstrapのための乱数を固定
set.seed(2023)

# Bootstrap のリサンプリング回数
numBootSample <- 100

# 各Bootstrap sampleで用いるマーケットのインデックスを乱数から発生させる
# 市場がNumSimMarkets個であるため、1からNumSimMarketsの整数について、重複を許してNumSimMarkets個ドローする
bootindex <- matrix(sample(1:NumSimMarkets, NumSimMarkets * numBootSample, replace = TRUE),
  ncol = numBootSample
)

# Matlabの再現: Bootstrap sampleで用いるマーケットのインデックス
if (isUseMatlabData == 1) {
  fa <- readMat(here("06_Dynamic_Game_Ch10_11/data_from_matlab/random_number_matlab_PSD.mat"))
  bootindex <- fa$bootindex
  numBootSample <- 100
  rm(fa)
}


# 結果を保存するための行列
bootresult_transition <- matrix(rep(0, 2 * numBootSample), ncol = numBootSample)
bootresult_CCP1 <- matrix(rep(0, 8 * numBootSample), ncol = numBootSample)
bootresult_CCP2 <- matrix(rep(0, 8 * numBootSample), ncol = numBootSample)
bootresult_payoff <- matrix(rep(0, 5 * numBootSample), ncol = numBootSample)

# Bootstrap を実行するループ
# mの繰り返し中にあるbootsampleの行を指定する部分の括弧は消してはいけない
# Estimation_AM_bootstrapが収束しない場合はスキップするので、numBootSampleよりも少ない数の結果が得られる可能性がある
b <- 1
while (b <= numBootSample) {
  print(append("Bootstrap:", as.character(b)))

  # Bootstrap sampleを構築する。
  bootsample <- matrix(rep(0, NumSimPeriods * NumSimMarkets * 8), ncol = 8)
  for (m in 1:NumSimMarkets) {
    temp <- FakeData[FakeData[, 1] == bootindex[m, b], ]
    bootsample[(1 + NumSimPeriods * (m - 1)):(NumSimPeriods * m), ] <- matrix(as.matrix(temp), ncol = 8)
  }

  # 構築したBootstrap sampleを用いて推定を行う
  # 関数Estimation_AM_bootstrapは4.のStep 1からStep3を実行する関数となっている
  output <- Estimation_AM_bootstrap(bootsample, beta)

  if (!is.null(output)) {
    # 解が得られた場合の処理
    # 推定結果を保存する
    bootresult_CCP1[, b] <- output[[1]][, 1]
    bootresult_CCP2[, b] <- output[[1]][, 2]
    bootresult_transition[, b] <- diag(output[[2]])
    bootresult_payoff[, b] <- output[[3]]
  } else {
    # 解が得られなかった場合
    # NAを保存する
    bootresult_CCP1[, b] <- NA
    bootresult_CCP2[, b] <- NA
    bootresult_transition[, b] <- NA
    bootresult_payoff[, b] <- NA
  }
  b <- b + 1 # b+1回目のブートストラップに進む
}

# NA以外の結果のみを先頭列（左端）から100列取り出す
na_columns <- colSums(is.na(bootresult_CCP1)) > 0
na_column_indices <- which(na_columns)
if (length(na_column_indices) > 30) {
  # NAが多すぎる（> 30）場合はエラーを出力する
  print("Error: too many NA columns")
  print(na_column_indices)
} else {
  # 100個以上のbootstrap sampleにおいて解が得られた場合
  print(na_column_indices)
  bootresult_CCP1 <- bootresult_CCP1[, !na_columns]
  bootresult_CCP2 <- bootresult_CCP2[, !na_columns]
  bootresult_transition <- bootresult_transition[, !na_columns]
  bootresult_payoff <- bootresult_payoff[, !na_columns]

  # 100個ちょうどになるように調整
  bootresult_CCP1 <- bootresult_CCP1[, 1:100]
  bootresult_CCP2 <- bootresult_CCP2[, 1:100]
  bootresult_transition <- bootresult_transition[, 1:100]
  bootresult_payoff <- bootresult_payoff[, 1:100]
}

# 結果のSummary
# CCP player 1
Summary_CCP1 <- matrix(c(CCP1UpdatedMat[, 2], EstimatedCCP1, apply(bootresult_CCP1, 1, sd, na.rm = TRUE)), ncol = 3)
colnames(Summary_CCP1) <- c("CCP for firm 1: True", "Estimated", "SE")
print(Summary_CCP1)

# CCP player 2
Summary_CCP2 <- matrix(c(CCP2UpdatedMat[, 2], EstimatedCCP2, apply(bootresult_CCP2, 1, sd, na.rm = TRUE)), ncol = 3)
colnames(Summary_CCP2) <- c("CCP for firm 2: True", "Estimated", "SE")
print(Summary_CCP2)

# Transition CCP
Summary_Transition <- matrix(c(c(0.7, 0.6), diag(EstimatedTransition), apply(bootresult_transition, 1, sd, na.rm = TRUE)), ncol = 3)
colnames(Summary_Transition) <- c("Transition Probability (GG and BB): True", "Estimated", "SE")
print(Summary_Transition)

# Payoff parameter
true <- c(Parameters[1:4], Parameters[6])
normalized_param <- c(Parameters[1] - (1 - beta) / beta * (Parameters[5]), Parameters[2] - (1 - beta) / beta * (Parameters[5]), Parameters[3], Parameters[4], Parameters[6] + (Parameters[5]))
Summary_Payoff_param <- matrix(c(true, normalized_param, sol$par, apply(bootresult_payoff, 1, sd, na.rm = TRUE)), ncol = 4)
colnames(Summary_Payoff_param) <- c("Payoff parameter: True", "Normalized true", "Estimated", "SE")
print(Summary_Payoff_param)

# Save results
write.csv(Summary_CCP1, file = here("06_Dynamic_Game_Ch10_11/output/Tab11_3_CCP_firm1.csv"))
write.csv(Summary_CCP2, file = here("06_Dynamic_Game_Ch10_11/output/Tab11_3_CCP_firm2.csv"))
write.csv(Summary_Transition, file = here("06_Dynamic_Game_Ch10_11/output/Tab11_2_Transition.csv"))
write.csv(Summary_Payoff_param, file = here("06_Dynamic_Game_Ch10_11/output/Tab11_4_AM2007.csv"))
