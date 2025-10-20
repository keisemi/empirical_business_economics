# Pesendorfer-SchmidDengler によるパラメタ推定

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
  pracma # matlabと同じ原理での最適化に用いる
)

# 必要な関数の読み込み
functionlt <- list.files("06_Dynamic_Game_Ch10_11/functions_R/", 
                         pattern="*.R$", full.names=TRUE, 
                         ignore.case=TRUE)
sapply(functionlt, source)


# 1. 下準備のコード。均衡計算を行う。
source(here("06_Dynamic_Game_Ch10_11/sub_1_prepare.R"))

# 2. 疑似データの作成 ----
# 疑似データ作成のコード。第11回と同じ。
source(here("06_Dynamic_Game_Ch10_11/sub_2_DGP.R"))

# 2. パラメタの推定----

# 紙面における点推定値を同じものを得るためには、
# Matlabで生成したFakedataを利用する必要がある。
# なおこの場合でも、Bootstraにより計算した標準誤差については乱数の影響で若干の差が出ることには留意されたい。

# isUseMatlabData = 0のときは上で生成したFakedataを利用する。
# isUseMatlabData = 1のときは、Matlabで生成したFakedata（FakeData_Matlab.csvファイル）を利用する。
isUseMatlabData <- 1

if (isUseMatlabData == 1){
  # Matlab で作ったFake dataをロードする。
  dt <- read.csv(here("06_Dynamic_Game_Ch10_11/data_from_matlab/FakeData_Matlab.csv"), header = FALSE)
  FakeData <- as.matrix(dt)
}

# 3. 推定 Step 1: CCPとTransition----
## Step 1: CCPと状態変数zの遷移確率を推定する----

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

# 確認
print("推定されたCCP")
print(cbind(EstimatedCCP1, EstimatedCCP2)) # 推定されたCCP

print("推定されたzの遷移確率行列")
print(EstimatedTransition) # 推定されたZの遷移確率行列

## Step 2: Pesendorfer and Schmidt-Dengler  の方法で推定を行う----

# ここで用いる関数が CCP_to_Value_to_prediction
# 引数：パラメタ、CCP、遷移確率、割引因子ベータ
# 内容：CCPから、価値関数をインバージョンし、その上で選択確率を計算する
# アウトプット：選択確率
# 疑似データ生成での手順の一部

### Step 2-1: Sanity check 1----
# DGPにおけるCCP、遷移確率、そして真のパラメタを入れて出てくる予測が、DGPのCCPと一致しているか確認する。
# 真のCCP: CCP1UpdatedMat, CCP2UpdatedMat
# Estimated CCP: EstimatedCCP1, EstimatedCCP2

output <-
  CCP_to_Value_to_prediction(TrueParameterValues, CCP1UpdatedMat[,2], 
                             CCP2UpdatedMat[,2], TransitionMat, beta)
diff <- output - cbind(CCP1UpdatedMat[,2], CCP2UpdatedMat[,2])
print("Difference: CCP in DGP and predicted CCP")
print(diff) # データ生成過程の真のCCPと予測されたCCPの差分

### Step 2-2: Sanity check 2----
# Agguiregabiria-SuzukiのNormalizationをしたパラメタの元で、同じ予測が出てくるかをチェックする。

Normalized_TrueParam <-
  matrix(c(Parameters[1] - ((1-beta)/beta)*Parameters[5], # 企業1のベース利潤
           Parameters[3], # 企業2の店舗数の企業1への影響 ;顧客収奪効果
           Parameters[4], # 景気が良い時の企業1の追加的利潤
           0, # 企業１の退出のためのコスト; 0に標準化
           Parameters[6] + Parameters[5], # 企業1の参入のためのコスト
           Parameters[2] - ((1-beta)/beta)*Parameters[5], # 企業2のベース利潤
           Parameters[3], # 企業1の店舗数の企業2への影響 ;顧客収奪効果
           Parameters[4], # 景気が良い時の企業2の追加的利潤
           0, # 企業2の退出のためのコスト; 0に標準化
           Parameters[6] + Parameters[5] # 企業2の参入のためのコスト
  ))

output_normalized <-
  CCP_to_Value_to_prediction(Normalized_TrueParam, CCP1UpdatedMat[,2], 
                             CCP2UpdatedMat[,2], TransitionMat, beta)
diff2 <- output_normalized - output
print("Difference: predicted CCP in true parameter and normalized parameter")
print(diff2) # 真のパラメタの元で予測されたCCPと標準化されたパラメタの差分

### Step 2-3: Estimate by least squares as in P-SD----

# obj_fun.Rで定義している関数obj_fun() は、パラメタ、CCP、遷移確率、割引因子を引数とする。
# 最適化のルーチンに入れるために、objという推定するパラメタのみに関する関数を定義する。
# 或いは、5行1列のパラメタをインプットとするobj_fun()も定義できる。

obj <- function(x){obj_fun(c(x[1], x[3:4], 0, x[5], x[2],x[3:4], 0, x[5]), 
                           EstimatedCCP1, EstimatedCCP2, EstimatedTransition,
                           beta)}

# 初期値の設定

# まず、初期値はTrue (標準化を加味する前のもの)とする
initial <- c(Parameters[1:4], Parameters[6])

# その上で、Trueに揺らぎ（0.6倍から1.2倍）を加えた5x10行列を作成する
mat_initial2 <- matrix(rep(initial, 10) * runif(50, min=0.6, max=1.2), nrow=5)

# ゆらぎを加えたものを初期値として推定を行う
result <- matrix(rep(0, 60), ncol = 10)

for (i in 1:10) {
  sol <- optim(mat_initial2[,i], obj)
  result[,i] <- matrix(c(sol$par, sol$value))
}

result_pick <- result[1:5, which.min(result[6,])]
cat(result_pick) # theta1, theta2, rival, z, entry costに対する推定結果

# matlab と同じ原理で最小化問題を解いた場合
# pracma::fminsearch() 
result_fmin = matrix(rep(0, 60), ncol = 10)
for (i in 1:10) {
  sol_fmin <- pracma::fminsearch(obj, mat_initial2[ ,i])
  result_fmin[,i] <- matrix(c(sol_fmin$xmin, sol_fmin$fmin))
}

result_fmin_pick <- result_fmin[1:5, which.min(result_fmin[6,])]

print("Estimation result:  (theta1, theta2, rival, z, entry cost)")
print(result_fmin_pick) # theta1, theta2, rival, z, entry costに対する推定結果

normalized_trueparam <- matrix(c(Normalized_TrueParam[1],
                                 Normalized_TrueParam[6], 
                                 Normalized_TrueParam[2],
                                 Normalized_TrueParam[3],
                                 Normalized_TrueParam[5]
))

print("Normalized true parameter following Agguiregabiria and Suzuki")
print(as.vector(normalized_trueparam)) # Aguirregabiria and Suzukiに従って標準化された真のパラメタ

## 5. Bootstrapによる標準誤差の計算 ----

# マーケット単位でリサンプリングを行う。マーケットは500個

# Bootstrapのための乱数を固定
set.seed(2023)

# Bootstrap のリサンプリング回数
numBootSample <- 100

# 各Bootstrap sampleで用いるマーケットのインデックスを乱数から発生させる。
# 市場がNumSimMarkets個であるため、1からNumSimMarketsの整数について、重複を許してNumSimMarkets個ドローする。
bootindex <- matrix(sample(1:NumSimMarkets, NumSimMarkets*numBootSample, replace = TRUE), 
                    ncol = numBootSample)

# Matlabの再現: Bootstrap sampleで用いるマーケットのインデックス
if (isUseMatlabData == 1) {
  fa <- readMat(here("06_Dynamic_Game_Ch10_11/data_from_matlab/random_number_matlab_PSD.mat"))
  bootindex <- fa$bootindex
  numBootSample <- 100
  rm(fa)
}

# 結果を保存するための行列
bootresult_transition <- matrix(rep(0,2*numBootSample) ,ncol = numBootSample)
bootresult_CCP1 <- matrix(rep(0,8*numBootSample) ,ncol = numBootSample)
bootresult_CCP2 <- matrix(rep(0,8*numBootSample) ,ncol = numBootSample)
bootresult_payoff <- matrix(rep(0,5*numBootSample) ,ncol = numBootSample)

# Bootstrap を実行するループ
# mの繰り返し中にあるbootsampleの行を指定する部分の括弧は消してはいけない。
for (b in 1:numBootSample) {
  print(append("Bootstrap:", as.character(b)))
  
  # Bootstrap sampleを構築する。
  bootsample <- matrix(rep(0,NumSimPeriods*NumSimMarkets*8) ,ncol = 8)
  for (m in 1:NumSimMarkets){
    temp <-FakeData[ FakeData[,1] == bootindex[m,b] , ]
    bootsample[(1+NumSimPeriods*(m-1)):(NumSimPeriods*m) , ] <- matrix(as.matrix(temp),ncol = 8)
  }
  
  
  # 構築したBootstrap sampleを用いて推定を行う。
  # 関数Estimation_PS_bootstrapは 上のStep 1とStep2を実行する関数となっている。
  output <- Estimation_PS_bootstrap(bootsample, beta)
  
  # 推定結果を保存する。
  
  bootresult_CCP1[,b] <- output[[1]][,1]
  bootresult_CCP2[,b] <- output[[1]][,2]
  bootresult_transition[,b] <- diag(output[[2]])
  bootresult_payoff[,b] <- output[[3]]
}

# 結果のSummary

# CCP player 1
Summary_CCP1<- matrix(c(CCP1UpdatedMat[ ,2], EstimatedCCP1,  apply(bootresult_CCP1,1,sd) ), ncol = 3)
colnames(Summary_CCP1) <- c("CCP for firm 1: True", "Estimated", "SE")
print(Summary_CCP1)

#CCP player 2
Summary_CCP2 <- matrix(c(CCP2UpdatedMat[ ,2], EstimatedCCP2,  apply(bootresult_CCP2,1,sd) ), ncol = 3)
colnames(Summary_CCP2) <- c("CCP for firm 2: True", "Estimated", "SE") 
print(Summary_CCP2)

#Transition Prob
Summary_Transition <- matrix( c(c(0.7, 0.6), diag(EstimatedTransition) , apply(bootresult_transition,1,sd)) , ncol = 3)
colnames(Summary_Transition) <- c("Transition Probability (GG and BB): True", "Estimated", "SE") 
print(Summary_Transition)

#Payoff parameter
true <- initial
normalized_param <- c( 0.3 - (1-beta)/beta*(-0.15), 0.2 - (1-beta)/beta*(-0.15),-0.27, 0.45, -2.1 + (-0.15))
Summary_Payoff_param <- matrix(c( true, normalized_param , result_pick, apply(bootresult_payoff,1,sd)) ,ncol = 4)
colnames(Summary_Payoff_param) <- c("Payoff parameter: True", "Normalized true", "Estimated", "SE")
print(Summary_Payoff_param)

# Save estimates of P-SD (1st stage estimates are the same as AM2007)
write.csv(Summary_Payoff_param, file = here("06_Dynamic_Game_Ch10_11/output/Tab11_4_PSD.csv"))
