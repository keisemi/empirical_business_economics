# 反実仮想シミュレーション

# 1. 下準備----
  
# データをクリア
rm(list = ls())

# メモリを開放
gc()

# 必要なパッケージの読み込み
# if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  # 必要となるパッケージ名
  matlib # 逆行列の作成に用いる
)

# 必要な関数の読み込み
functionlt <- list.files("./functions", 
                         pattern="*.R$", full.names=TRUE, 
                         ignore.case=TRUE)
sapply(functionlt, source)

# 割引因子
beta <- 0.8

# 定数
eulergamma <- 0.5772

# 景気の遷移行列
TransitionMat <- matrix(c(0.7,0.3, 0.4, 0.6), nrow = 2, byrow = TRUE)

# 今回のモデルにおける、各企業の持っている全選択要素数
NumMaxChoices <- 3

Parameters <- matrix(c(
  0.3, # 企業1のベース利潤
  0.2, # 企業2のベース利潤
  -0.27, # 顧客収奪効果
  0.45, # 景気が良い時の追加的利潤
  -0.15, # 退出のためのコスト
  -2.10 # 参入のためのコスト
))
NewParameters <- matrix(c(
  0.3, # 企業1のベース利潤
  0.2, # 企業2のベース利潤
  -0.25, # 顧客収奪効果
  0.5, # 景気が良い時の追加的利潤
  -0.15, # 退出のためのコスト
  -2 # 参入のためのコスト
))

# パラメターの設定 (ベースライン)
BaselineParameterValues <- matrix(c(
  Parameters[1], # 企業1のベース利潤
  Parameters[3], # ライバルの店舗数が企業1の利潤に与える影響
  Parameters[4], # 景気が良い時の企業1への追加的利潤
  Parameters[5], # 企業1の退出のためのコスト
  Parameters[6], # 企業1の出店のためのコスト
  Parameters[2], # 企業2のベース利潤
  Parameters[3], # ライバルの店舗数が企業2の利潤に与える影響
  Parameters[4], # 景気が良い時の企業2への追加的利潤
  Parameters[5], # 企業2の退出のためのコスト
  Parameters[6]  # 企業2の出店のためのコスト
))

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

pi1 <- pi1gen(BaselineParameterValues) * CCP1Adjuster 
pi2 <- pi2gen(BaselineParameterValues) * CCP2Adjuster


# 2. ベースラインの均衡----

# MPEを解く
output <- f_MPE(TransitionMat, pi1, pi2, beta) #CCP1base, CCP2base, V1base, V2base

CCP1base <- output[[1]]
CCP2base <- output[[2]]
V1base <- output[[3]]
V2base <- output[[4]]
# dataをシミュレーションする
NumSimPeriods <- 15

# Transition matrix. Use fPsigma.R

fPsigma <- fP(TransitionMat, CCP1base[ ,2], CCP2base[ ,2] )

initial_state <- c(1, rep(0,7))

transitionpath <- matrix(rep(0,NumSimPeriods*8),ncol = 8)

for (tt in 1:NumSimPeriods){
  if (tt == 1){
    transitionpath[tt, ] <- t(initial_state)
  } else if (tt >= 2 ){
    transitionpath[tt, ] <- t( t(fPsigma) %*% t(matrix(transitionpath[tt-1, ],ncol = 8)))
  }
}

# 参入状況をカウントする。

n1 <- apply(transitionpath * matrix(rep(c(0,0,1,1,0,0,1,1),NumSimPeriods),ncol = 8,byrow = TRUE), 1,sum)
n2 <- apply(transitionpath * matrix(rep(c(0,1,0,1,0,1,0,1),NumSimPeriods),ncol = 8,byrow = TRUE), 1,sum)


## 4. 反実仮想シミュレーション: 差別化戦略

# シナリオ1
# 企業1のベース利潤が上昇し、顧客奪取効果がどちらも0になった場合
CounterfactualParameterValues1 <- matrix(c(
  0.5, # 企業1のベース利潤
  0, # ライバルの店舗数が企業1の利潤に与える影響
  Parameters[4], # 景気が良い時の企業1への追加的利潤
  Parameters[5], # 企業1の退出のためのコスト
  Parameters[6], # 企業1の出店のためのコスト
  Parameters[2], # 企業2のベース利潤
  0, # ライバルの店舗数が企業2の利潤に与える影響
  Parameters[4], # 景気が良い時の企業2への追加的利潤
  Parameters[5], # 企業2の退出のためのコスト
  Parameters[6]  # 企業2の出店のためのコスト
))


# Profitの式を作成
pi1_cf1 <- pi1gen(CounterfactualParameterValues1) * CCP1Adjuster
pi2_cf2 <- pi2gen(CounterfactualParameterValues1) * CCP2Adjuster

# MPEを解く
output_cf <- f_MPE(TransitionMat,  pi1_cf1, pi2_cf2, beta) 

CCP1cf <- output_cf[[1]]
CCP2cf <- output_cf[[2]]
V1cf <- output_cf[[3]]
V2cf <- output_cf[[4]]

# Transition matrix. Use fPsigma.R

fPsigma <- fP(TransitionMat, CCP1cf[ ,2], CCP2cf[ ,2])

initial_state <- c(1, rep(0,7))

transitionpathcf <- matrix(rep(0,NumSimPeriods*8),ncol = 8)

for (tt in 1:NumSimPeriods){
  if (tt == 1){
    transitionpathcf[tt, ] <- t(initial_state)
  } else if (tt >= 2 ){
    transitionpathcf[tt, ] <- t( t(fPsigma) %*% t(matrix(transitionpath[tt-1, ],ncol = 8)))
  }
}

# 参入状況をカウントする。
n1_cf <- apply(transitionpathcf * matrix(rep(c(0,0,1,1,0,0,1,1),NumSimPeriods),ncol = 8,byrow = TRUE), 1,sum)
n2_cf <- apply(transitionpathcf * matrix(rep(c(0,1,0,1,0,1,0,1),NumSimPeriods),ncol = 8,byrow = TRUE), 1,sum)
n1_cf_s1 <- n1_cf
n2_cf_s1 <- n2_cf

# par(mfrow=c(1,2))関数を使用して2つのプロットを並べる  
png("result/ProbEntry1Orginal.png")
par(mfrow=c(1,2))
plot(1:NumSimPeriods, n1, type="l",xlim = c(1, 15), ylim = c(0, 0.8)
     ,xaxp=c(1, 15, 14),xaxt="n" ,main="企業1の店舗存在確率", xlab="", ylab="")
lines(1:NumSimPeriods, n1_cf, type="b",lty = "dashed", pch = 18,xaxp=c(1, 15, 14),xaxt="n")
axis(1, at=1:15)

plot(1:NumSimPeriods, n2, type="l",xlim = c(1, 15), ylim = c(0, 0.8)
     ,xaxp=c(1, 15, 14),xaxt="n",main="企業2の店舗存在確率", xlab="", ylab="")
lines(1:NumSimPeriods, n2_cf, type="b",lty = "dashed", pch = 18,xaxp=c(1, 15, 14),xaxt="n")
axis(1, at=1:15)
dev.off()

# 利潤の計算
diff_value_1 <- V1cf - V1base;
matrix(c(V1base, V1cf, diff_value_1),ncol = 3)
diff_value_2 <- V2cf - V2base;
matrix(c(V2base, V2cf, diff_value_2),ncol = 3)
print("シナリオ1")
matrix(c(V1base, V2base, V1cf, diff_value_1, V2cf),ncol = 5)


# シナリオ2
# 企業1のベース利潤が上昇し、顧客奪取効果に企業間で異質性がある場合
CounterfactualParameterValues2 <- matrix(c(
  0.5, # 企業1のベース利潤
  -0.1, # ライバルの店舗数が企業1の利潤に与える影響
  Parameters[4], # 景気が良い時の企業1への追加的利潤
  Parameters[5], # 企業1の退出のためのコスト
  Parameters[6], # 企業1の出店のためのコスト
  Parameters[2], # 企業2のベース利潤
  -0.2, # ライバルの店舗数が企業2の利潤に与える影響
  Parameters[4], # 景気が良い時の企業2への追加的利潤
  Parameters[5], # 企業2の退出のためのコスト
  Parameters[6]  # 企業2の出店のためのコスト
))

# Profitの式を作成
pi1_cf1 <- pi1gen(CounterfactualParameterValues2) * CCP1Adjuster
pi2_cf2 <- pi2gen(CounterfactualParameterValues2) * CCP2Adjuster

# MPEを解く
output_cf <- f_MPE(TransitionMat,  pi1_cf1, pi2_cf2, beta) 

CCP1cf <- output_cf[[1]]
CCP2cf <- output_cf[[2]]
V1cf <- output_cf[[3]]
V2cf <- output_cf[[4]]

# Transition matrix. Use fPsigma.R

fPsigma <- fP(TransitionMat, CCP1cf[ ,2], CCP2cf[ ,2])

initial_state <- c(1, rep(0,7))

transitionpathcf <- matrix(rep(0,NumSimPeriods*8),ncol = 8)

for (tt in 1:NumSimPeriods){
  if (tt == 1){
    transitionpathcf[tt, ] <- t(initial_state)
  } else if (tt >= 2 ){
    transitionpathcf[tt, ] <- t( t(fPsigma) %*% t(matrix(transitionpath[tt-1, ],ncol = 8)))
  }
}

# 参入状況をカウントする。
n1_cf <- apply(transitionpathcf * matrix(rep(c(0,0,1,1,0,0,1,1),NumSimPeriods),ncol = 8,byrow = TRUE), 1,sum)
n2_cf <- apply(transitionpathcf * matrix(rep(c(0,1,0,1,0,1,0,1),NumSimPeriods),ncol = 8,byrow = TRUE), 1,sum)
n1_cf_s2 <- n1_cf
n2_cf_s2 <- n2_cf

# par(mfrow=c(1,2))関数を使用して2つのプロットを並べる  
png("result/ProbEntry2Original.png")
par(mfrow=c(1,2))
plot(1:NumSimPeriods, n1, type="l",xlim = c(1, 15), ylim = c(0, 0.8)
     ,xaxp=c(1, 15, 14),xaxt="n" ,main="企業1の店舗存在確率", xlab="", ylab="")
lines(1:NumSimPeriods, n1_cf, type="b",lty = "dashed", pch = 18,xaxp=c(1, 15, 14),xaxt="n")
axis(1, at=1:15)

plot(1:NumSimPeriods, n2, type="l",xlim = c(1, 15), ylim = c(0, 0.8)
     ,xaxp=c(1, 15, 14),xaxt="n",main="企業2の店舗存在確率", xlab="", ylab="")
lines(1:NumSimPeriods, n2_cf, type="b",lty = "dashed", pch = 18,xaxp=c(1, 15, 14),xaxt="n")
axis(1, at=1:15)
dev.off()

# 利潤の計算
diff_value_1 <- V1cf - V1base;
matrix(c(V1base, V1cf, diff_value_1),ncol = 3)
diff_value_2 <- V2cf - V2base;
matrix(c(V2base, V2cf, diff_value_2),ncol = 3)
print("シナリオ2")
matrix(c(V1base, V2base, V1cf, diff_value_1, V2cf),ncol = 5)

# 3プロット
png("result/ProbEntry3Plots.png")
par(mfrow=c(1,2))
plot(1:NumSimPeriods, n1, type="l",xlim = c(1, 15), ylim = c(0, 0.8)
     ,xaxp=c(1, 15, 14),xaxt="n" ,main="企業1の店舗存在確率", xlab="", ylab="")
lines(1:NumSimPeriods, n1_cf_s1, type="b",lty = "dashed", pch = 1,xaxp=c(1, 15, 14),xaxt="n")
lines(1:NumSimPeriods, n1_cf_s2, type="b",lty = "dashed", pch = 16,xaxp=c(1, 15, 14),xaxt="n")
axis(1, at=1:15)

plot(1:NumSimPeriods, n2, type="l",xlim = c(1, 15), ylim = c(0, 0.8)
     ,xaxp=c(1, 15, 14),xaxt="n",main="企業2の店舗存在確率", xlab="", ylab="")
lines(1:NumSimPeriods, n2_cf_s1, type="b",lty = "dashed", pch = 1,xaxp=c(1, 15, 14),xaxt="n")
lines(1:NumSimPeriods, n2_cf_s2, type="b",lty = "dashed", pch = 16,xaxp=c(1, 15, 14),xaxt="n")
axis(1, at=1:15)
dev.off()