
rm(list = ls())

library(tidyverse)
library(optimx)
library(here)

source(here("04_Static_Game_MRI_Ch08/function.R"))

data <- read_csv(file = here("04_Static_Game_MRI_Ch08/data/MRIData_BR1991.csv"), locale = locale(encoding = "UTF-8"))
data <- data.frame(data)

# データの準備 ----

# 病院レベルのデータを市区町村レベルの集計データに変換
listCode <- unique(data$CityCode)

NumHospital <- numeric(length(listCode))
NumMRI <- numeric(length(listCode))
Pop <- numeric(length(listCode))
Menseki <- numeric(length(listCode))
PopDen <- numeric(length(listCode))
Income <- numeric(length(listCode))


for (i in 1:length(listCode)){
  subdata <- data[data$CityCode==listCode[i],]
  NumHospital[i] <- length(subdata$MRIOwnDum)
  NumMRI[i] <- sum(subdata$MRIOwnDum)
  Pop[i] <- as.integer(unique(subdata$Population))
  Menseki[i] <- as.integer(unique(subdata$Menseki)[1])
  PopDen[i] <- unique(subdata$PopDensity)
  Income[i] <- as.integer(unique(subdata$TaxableIncome))
}

dataset <- data.frame(Code = listCode, 
                      NumHospital = NumHospital,
                      NumMRI = NumMRI, 
                      Pop = Pop,
                      Menseki = Menseki,
                      PopDen = PopDen, 
                      Income = Income)

# Tab 8.3: 病院数をMRI保有病院数のテーブル (変数欠損を落とす前のサンプル)
# なお、ここでは10以上は10としてまとめる。
dataset %>% 
  select(NumHospital, NumMRI) %>%
  mutate(NumHospital = ifelse(NumHospital > 10 , 10, NumHospital), 
         NumMRI = ifelse(NumMRI > 10, 10, NumMRI)) -> dt_table

# テーブルを出力
tbl <- table(dt_table$NumHospital, dt_table$NumMRI)
print(tbl)

# CSVとして保存
write.csv(tbl, file = here("04_Static_Game_MRI_Ch08/output/Tab8_3_hospital_mri_table.csv"), row.names = TRUE)

# Bresnahan and Reiss (1991b)モデルの推定 ----

# 変数に欠損がある市区町村は落とす
dataset <- na.omit(dataset)

# 人口を百万で除する
dataset$Pop <- dataset$Pop/1000000

# 最終的に推定で用いるサンプルサイズ（市区町村数）の確認
M <- nrow(dataset)
print(M)

# N_maxを6から10まで試す。
N_max_list <- 6:10


# dataset をキープする。
dataset_org <- dataset

# 結果の保管場所
result_list_estimates = list()
result_list_thresholds = list()

# ループ: 
for (i in 6:10){ # N_max = 6,...,10
  
  # 表示
  print(paste("N_max:", i))
  
  # N_max を設定  
  N_max <- i 
  
  # あとでリストに使うName
  name <- paste0("N_max=", i)
  
  # datasetを初期化
  dataset <- dataset_org
  
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
  estresult <- cbind(estimates, se)
  print("Estimates")
  print(estresult)
  
  # 結果を保管
  result_list_estimates[[name]] <- estresult 
　  
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
  
  # 比率
  s <- EntryThreshold[,2]
  ratio <- numeric(length(s))
  for (j in 1:length(s)){
    
    if (j < length(s)){
      ratio[j] <- s[j+1]/s[j]
    } else {
      ratio[j] <- NA
    }
  }
  
  # Combineする
  EntryThreshold <- cbind(EntryThreshold, ratio) 
  colnames(EntryThreshold) <- c("S_N", "s_N=S_N /N", "s_{N+1}/s_N (ratio) ")

  # 行列を表示する
  print("Entry Threshold")
  print(EntryThreshold)
  
  # 結果を保管
  result_list_thresholds[[name]] <- EntryThreshold 
  
  
  
}

# 結果を保存
# 推定値
sink(here("04_Static_Game_MRI_Ch08/output/Tab8_4_BR1991_Estimates.txt")  )
print(result_list_estimates)
sink()

# Entry threshold
sink(here("04_Static_Game_MRI_Ch08/output/Tab8_5_Entry_Thresholds.txt")  )
print(result_list_thresholds)
sink()




