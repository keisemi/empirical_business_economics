
library(tidyverse)
library(optimx)
library(here)

source(here("04_Static_Game_MRI_Ch08/function.R"))

data <- read_csv(file = here("04_Static_Game_MRI_Ch08/data/MRIData.csv"), locale = locale(encoding = "UTF-8"))
data <- data.frame(data)

# データの準備 ----

# 病院レベルのデータを市区町村レベルの集計データに変換
listCode <- unique(data$CityCode)

NumMRI <- numeric(length(listCode))
Pop <- numeric(length(listCode))
Menseki <- numeric(length(listCode))
PopDen <- numeric(length(listCode))
Income <- numeric(length(listCode))

for (i in 1:length(listCode)){
  subdata <- data[data$CityCode==listCode[i],]
  NumMRI[i] <- sum(subdata$MRIOwnDum)
  Pop[i] <- as.integer(unique(subdata$Population))
  Menseki[i] <- as.integer(unique(subdata$Menseki)[1])
  PopDen[i] <- unique(subdata$PopDensity)
  Income[i] <- as.integer(unique(subdata$TaxableIncome))
}

dataset <- data.frame(Code = listCode, 
                      NumMRI = NumMRI, 
                      Pop = Pop,
                      Menseki = Menseki,
                      PopDen = PopDen, 
                      Income = Income)


# 変数に欠損がある市区町村は落とす
dataset <- na.omit(dataset)

# 人口を百万で除する
dataset$Pop <- dataset$Pop/1000000

# 最終的に推定で用いるサンプルサイズ（市区町村数）の確認
M <- nrow(dataset)
print(M)

# Bresnahan and Reiss (1991b)モデルの推定 ----

N_max <- 6

# MRIを導入している病院数がN_maxより大きい市区町村のMRI導入済み病院数を、N_maxに置換する
dataset$NumMRI[dataset$NumMRI > N_max] <- N_max


# パラメターの初期値を全て1に設定
initial <- rep(1, N_max+1)

# 最適化
res <- optimx(par = initial, 
              fn = obj, 
              dataset = dataset, 
              N_max = N_max, 
              hessian=T, 
              method='L-BFGS-B', 
              lower=0,  
              control = list(fnscale=-1))


# ヘシアンの計算
Hessian <- gHgen(par = as.numeric(res[0:N_max+1]), obj, dataset = dataset, N_max = N_max) 
# 標準誤差の確認
se <- sqrt(diag(solve(-Hessian$Hn)/M)) 

# 推定値と標準誤差を表示
# 上からalpha_0, alpha_1, ...で、最終行がgamma_0の推定値・標準誤差
estimates <- as.numeric(res[0:N_max+1])
cbind(estimates, se)

# 推定値に基づくエクササイズ ----

# 推定値からalphaとgammaの値を定義する
alpha <- estimates[1:N_max]
gamma <- estimates[N_max+1]

# Entry Threshold の行列を定義する
EntryThreshold <- matrix(0, N_max, 2) 

# Entry Threshold(S_1とs_1)を計算して、行列に代入する
deno <- alpha[1]
EntryThreshold[1,] <-  c(as.integer(gamma/deno*10^6), as.integer(gamma/deno*10^6))

# Entry Threshold(S_n, s_n, n>1)を計算して、行列に代入する
for (i in 2:N_max){
  deno <- deno-alpha[i]
  EntryThreshold[i,] <- c(as.integer(gamma/deno*10^6), as.integer(gamma/deno*10^6/i))
}

# 行列を表示する
print(EntryThreshold)


