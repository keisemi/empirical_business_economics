# 1. はじめに----
# このスクリプトは経済セミナー第11回の資料を拡張し、
# Aguirregabiria and Mira (2007)のモデルに基づいた
# 経済セミナーのモデルの推定を行うものである
# したがって、疑似データの作成や基本的な関数の定義に関しては当該資料と同じものを用いている
# ただし、真のパラメータは新たに設定したものを用いている
#
# 本スクリプトの構成は以下のようになっている:
# 1. はじめに
# 2. 下準備（パラメータ以外は経済セミナー資料と同じ）
# 3. 疑似データの生成（経済セミナー資料と同じ）
# 4. Aguirregabiria and Mira (2007)の方法によるパラメターの推定
# 5. Bootstrapによる標準誤差の計算

# 2. 下準備----

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
  pracma # matlabと同じ原理での最適化に用いる
)

# 必要な関数の読み込み
functionlt <- list.files("06_Dynamic_Game_Ch10_11/functions_R/",
  pattern = "*.R$", full.names = TRUE,
  ignore.case = TRUE
)
sapply(functionlt, source)

# 割引因子
beta <- 0.8

# 定数
eulergamma <- 0.5772

# 景気の遷移行列
TransitionMat <- matrix(c(0.7, 0.3, 0.4, 0.6), nrow = 2, byrow = TRUE)
# デフォルトではbyrowがFALSEになっており、列方向に各要素が代入される

# パラメターの設定
# OldParametersは経済セミナーの資料において設定されたパラメター
Parameters <- matrix(c(
  0.3, # 企業1のベース利潤
  0.2, # 企業2のベース利潤
  -0.27, # 顧客収奪効果
  0.45, # 景気が良い時の追加的利潤
  -0.15, # 退出のためのコスト
  -2.10 # 参入のためのコスト
))
# Parametersは今回新たに設定したパラメター
NewParameters <- matrix(c(
  0.3, # 企業1のベース利潤
  0.2, # 企業2のベース利潤
  -0.25, # 顧客収奪効果
  0.5, # 景気が良い時の追加的利潤
  -0.15, # 退出のためのコスト
  -2 # 参入のためのコスト
))
# NewParametersを利用する場合
# Parameters <- NewParameters

# matrix(c())はncolやnrowを指定しない場合、n行1列の行列を作成する
# 設定したパラメターを一般化しやすいよう並べ替える
TrueParameterValues <- matrix(c(
  Parameters[1], # 企業1のベース利潤
  Parameters[3], # ライバルの店舗数が企業1の利潤に与える影響
  Parameters[4], # 景気が良い時の企業1への追加的利潤
  Parameters[5], # 企業1の退出のためのコスト
  Parameters[6], # 企業1の出店のためのコスト
  Parameters[2], # 企業2のベース利潤
  Parameters[3], # ライバルの店舗数が企業2の利潤に与える影響
  Parameters[4], # 景気が良い時の企業2への追加的利潤
  Parameters[5], # 企業2の退出のためのコスト
  Parameters[6] # 企業2の出店のためのコスト
))

# 後々の計算を簡略化するための行列として、CCP1AdjusterとCCP2Adjusterを
# 定義する。これらは、以下のように、各状態変数で選択することができる選択肢には
# 1が、選択不可能な選択肢には0が割り振られている行列である
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

# 真のパラメターの下で、pi1とpi2を計算する。各々は8*3の行列になっていて、縦が
# 状態変数（G00,G01,G10,...,B11）を意味し、横がa_i = -1, 0, 1 という行動を
# 意味し、それぞれの状態変数で行動を選んだ時の利潤が行列として与えられている。
# 一般的には、企業1の利潤関数はpi_1(a_1,a_2,s)のように、自身の行動a_1だけで
# なく、企業2の行動a_2にも依存するが、今回の特定化ではお互いの企業の利潤は、
# ライバル企業の行動の影響を受けない関数形になっていることに注意

pi1 <- pi1gen(TrueParameterValues) * CCP1Adjuster # アダマール積
pi2 <- pi2gen(TrueParameterValues) * CCP2Adjuster

# 3. 疑似データの生成----

## Step 1: CCPの初期値の設定（および設定したCCPをベクトルから行列へ変換する）----

# 以下ではCCPを全て0.5とする
CCP1 <- matrix(0.5, 8, 1)
CCP2 <- matrix(0.5, 8, 1)
# ただし、この行列は同じ確率である必要はなく、例えば以下のような行列でも構わない
# CCP1 <- matrix(c(.5, .55, .6, .65, .4, .6, .45, .55))
# CCP2 <- matrix(c(.5, .6, .55, .65, .4, .45, .6, .55))

# 以下の関数を用いて上で与えた初期値のベクトルを行列へと変換する
# 例えば、CCP1Matは
# CCP1Mat = [ 0   0.5 0.5
# 0   0.5 0.5
# :    :   :
# 0.5  0.5  0
# 0.5  0.5  0  ] のような行列になる
CCP1Mat <- CCP1Transform(CCP1)
CCP2Mat <- CCP2Transform(CCP2)

## Step 2: Step 1 で与えられたCCPの初期値を基に、事前の価値関数を計算する----

# 基本的には（４）式の計算を行いたいため、以下で必要なパーツを計算する

# まず、F^P^sigmaを求める。これは109ページの囲みで与えられているように、
# 遷移行列とCCP1とCCP2から作られる３つの行列の、要素ごとの積（アダマール
# 積）として計算できる
fPsigma <- fP(TransitionMat, CCP1, CCP2)

# パラメターの下で、pi1Psigma と pi2Psigma を計算する

pi1Psigma <- pi1PsigmaGen(pi1, CCP2Mat)
pi2Psigma <- pi2PsigmaGen(pi2, CCP1Mat)

# e^P_i(a_i,s)を計算する
# Rでは自動的に行列の次元が設定される
eP1 <- eulergamma - CCP1LogTransform(CCP1)
eP2 <- eulergamma - CCP2LogTransform(CCP2)
# 少々長くはなるが、以下のように書くことも可能
# eP1 <- eulergamma * matrix(rep(1, 8*3), ncol = 3) - CCP1LogTransform(CCP1)
# eP2 <- eulergamma * matrix(rep(1, 8*3), ncol = 3) - CCP2LogTransform(CCP2)

# 上で求めたものを（４）式に代入して、事前の価値関数を求める
ExanteV1 <- as.matrix(matlib::inv(diag(8) - beta * fPsigma) %*%
  ((CCP1Mat * (pi1Psigma + eP1)) %*% matrix(c(1, 1, 1))))
ExanteV2 <- matlib::inv(diag(8) - beta * fPsigma) %*%
  ((CCP2Mat * (pi2Psigma + eP2)) %*% matrix(c(1, 1, 1)))
# (pij + ePj) が何らかの値をとっていたとしても、CCP1MatやCCP2Matを
# かけることで、取ることのできない選択肢にはゼロが振られることになり、
# 分母の和は正しく計算できていることがわかる。例えば、
# CCP1Mat.*(pi1Psigma+eP1) を表示させてみると、実際にいくつかの
# 要素がゼロになっていることが確認できる

## Step 3: 事前の価値関数を基に、CCPを計算する ----

# 各企業の行動が与えられたときの遷移行列を求める。これは、本文では
# f^P*(s'|s,a_i)として表現されている行列である。fP_a1/2はリスト形式になっているが
# fP_a1/2[[1]] がa1/2=-1であるときの、
# fP_a1/2[[2]] がa1/2=0 であるときの、
# fP_a1/2[[3]] がa1/2=1 であるときの遷移行列に相当する

fP_a1 <- fP_a1given(TransitionMat, CCP2)
fP_a2 <- fP_a2given(TransitionMat, CCP1)

# 以下では、求められた事前の価値関数をもとにCCPを更新する。最終的には、
# 計算の簡便化のため更新されたCCPを行列の形（CCP1UpdatedMat）とベクトル
# の形（CCP1Updated）の両方で求めておく。この計算は一見複雑だが、需要
# 関数の推定の時と似た計算であるため、説明は省略する。最初のブロックは企
# 業1について、次のブロックは企業2について更新されたCCPを求めている

# fP_a1given()は出力されるものがlist形式であるに注意
# fP_a1はlistであるため、この様にかっこをつけなければならない
# []を一つにするとclassがlistのままで行列の積が計算できない
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


## Step 4: アップデートされたCCPを基に、再び事前の価値関数を計算する----
# ただし、 Step 2と同一の作業であるため、説明は省略する

fPsigma <- fP(TransitionMat, CCP1Updated, CCP2Updated)

pi1Psigma <- pi1PsigmaGen(pi1, CCP2UpdatedMat)
pi2Psigma <- pi2PsigmaGen(pi2, CCP1UpdatedMat)

eP1 <- eulergamma - CCP1LogTransform(CCP1Updated)
eP2 <- eulergamma - CCP2LogTransform(CCP2Updated)

ExanteV1Updated <- inv((diag(8) - beta * fPsigma)) %*% apply(CCP1UpdatedMat * (pi1Psigma + eP1), 1, sum)
ExanteV2Updated <- inv((diag(8) - beta * fPsigma)) %*% apply(CCP2UpdatedMat * (pi2Psigma + eP2), 1, sum)

## Step 5: Step 3 と Step 4 を事前の価値関数が一致するまで繰り返す----

# 前の事前の価値関数と更新された事前の価値関数の差を定義する
DiffExanteV <- apply((ExanteV1Updated - ExanteV1)^2 + (ExanteV2Updated - ExanteV2)^2, 2, sum)

# 上で定義した差が1.0e-12よりも小さくなるまで、繰り返す
while (DiffExanteV > 1.0e-12) {
  # 更新されたCCPと事前の価値関数を、CCPとExanteVとして置き換える
  CCP1 <- CCP1Updated
  CCP2 <- CCP2Updated

  ExanteV1 <- ExanteV1Updated
  ExanteV2 <- ExanteV2Updated

  # Step 3 を再度実行する（上と全く一緒のため、説明は省略する）
  fP_a1 <- fP_a1given(TransitionMat, CCP2)
  fP_a2 <- fP_a2given(TransitionMat, CCP1)

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

  # Step 4 を再度実行する（上と全く一緒のため、説明は省略する）
  fPsigma <- fP(TransitionMat, CCP1Updated, CCP2Updated)

  pi1Psigma <- pi1PsigmaGen(pi1, CCP2UpdatedMat)
  pi2Psigma <- pi2PsigmaGen(pi2, CCP1UpdatedMat)

  eP1 <- eulergamma - CCP1LogTransform(CCP1Updated)
  eP2 <- eulergamma - CCP2LogTransform(CCP2Updated)

  ExanteV1Updated <- inv((diag(8) - beta * fPsigma)) %*% apply(CCP1UpdatedMat * (pi1Psigma + eP1), 1, sum)
  ExanteV2Updated <- inv((diag(8) - beta * fPsigma)) %*% apply(CCP2UpdatedMat * (pi2Psigma + eP2), 1, sum)

  # Step 5 の冒頭部分の差分の計算を再度実行する
  DiffExanteV <- apply((ExanteV1Updated - ExanteV1)^2 + (ExanteV2Updated - ExanteV2)^2, 2, sum)
}

# 上記のwhileループを実行することで、均衡のCCPを求めることができる。
# ここで、確認のため、どのようなCCPになっているかを表示させる
print(matrix(c(CCP1UpdatedMat, CCP2UpdatedMat), ncol = 6)) # 均衡でのCCP

## Step 6: 均衡におけるCCPをもとに、疑似データをシミュレートする----

# まず、シミュレーションの設定を行う。
# 乱数の発生が必要なため、乱数発生のためのシードを設定する
set.seed(2023)

# 次に、何個の市場で何期のデータを生成するかを設定し、企業数は2とする
NumSimMarkets <- 500 # 市場の数
NumSimPeriods <- 50 # シミュレーションの期間
NumSimFirms <- 2 # 企業の数

# さらに、各市場の1期の状態変数を、一様分布に従ってランダムに生成する
InitialState <- sample(1:8, NumSimMarkets, replace = TRUE)

# 企業1と企業2がどのような行動をとるのか、及び、次期の状態変数を決める際
# に必要となる一様乱数を RandomNumbers という行列に格納する
RandomNumbers <- array(runif(NumSimMarkets * NumSimPeriods * (NumSimFirms + 1), min = 0, max = 1),
  dim = c(NumSimMarkets, NumSimPeriods, NumSimFirms + 1)
)

# 最後に疑似データを格納するために、FakeDataという行列を準備する
FakeData <- matrix(rep(0, NumSimMarkets * NumSimPeriods * 8), ncol = 8)
# 8列あるが、それぞれの列は以下のデータを格納している
# 1: Market ID
# 2: 時間t, t=1,...,NumSimPeriods
# 3: 状態変数 \in {1, 2, 3, 4, 5, 6, 7, 8}
# 4: 状態変数の要素である景気状態
# 5: 状態変数の要素である出店数 n_{1,t}
# 6: 状態変数の要素である出店数 n_{2,t}
# 7: 生成された企業1の行動 a_{1,t}
# 8: 生成された企業2の行動 a_{2,t}


for (m in 1:NumSimMarkets) {
  for (t in 1:NumSimPeriods) {
    # まず、市場のIDと時点tの情報をFakeDataの1列目と2列目に格納する
    FakeData[(m - 1) * NumSimPeriods + t, 1] <- m
    FakeData[(m - 1) * NumSimPeriods + t, 2] <- t

    # 1期目の状態変数は先の通り乱数で発生させているが、2期目以降は
    # 内生的に決まることから、ここでは1期目と2期目以降は別々のプロ
    # グラムとなっている
    if (t == 1) {
      # 先に決めたInitialState(m)から、該当する市場の初期状態を
      # FakeDataの3列目に格納する
      FakeData[(m - 1) * NumSimPeriods + 1, 3] <- InitialState[m]

      # 格納した初期状態に応じて、景気の状態がどうなっているかを
      # 判断して、FakeDataの4列目に格納する
      if (InitialState[m] >= 1 && InitialState[m] <= 4) {
        FakeData[(m - 1) * NumSimPeriods + 1, 4] <- 1
      } else {
        FakeData[(m - 1) * NumSimPeriods + 1, 4] <- 2
      }
    } else {
      # t=2,...,NumSimPeriodsの場合, まずは前期の状態変数を一行
      # 前を参照にして sprev (previous state)として持ってくる。
      # さらに、今期の状態変数は前期の各企業の行動にも依存するので、
      # 前期の各企業の行動をa1prev, a2prevとしてもってくる
      sprev <- FakeData[(m - 1) * NumSimPeriods + t - 1, 3]
      a1prev <- FakeData[(m - 1) * NumSimPeriods + t - 1, 7]
      a2prev <- FakeData[(m - 1) * NumSimPeriods + t - 1, 8]

      # 以下では、まず前期の状態変数から、景気の状態についての情報
      # 取り出し、一様乱数を用いて今期の景気の状態を決める。その上
      # で、今期の景気の状態を4行目に格納する
      if (sprev >= 1 && sprev <= 4) {
        if (RandomNumbers[m, t, 3] < TransitionMat[1, 1]) {
          FakeData[(m - 1) * NumSimPeriods + t, 4] <- 1
        } else {
          FakeData[(m - 1) * NumSimPeriods + t, 4] <- 2
        }
      } else {
        if (RandomNumbers[m, t, 3] < TransitionMat[2, 2]) {
          FakeData[(m - 1) * NumSimPeriods + t, 4] <- 2
        } else {
          FakeData[(m - 1) * NumSimPeriods + t, 4] <- 1
        }
      }

      # さらに、今期の景気の状態（snow）と、今期の各企業の出店数
      # （n1tとn2t）を基に、今期の状態変数が何かを判断して、
      # FakeDataの3列目に格納していく
      snow <- FakeData[(m - 1) * NumSimPeriods + t, 4]
      n1t <- FakeData[(m - 1) * NumSimPeriods + t - 1, 5] + a1prev
      n2t <- FakeData[(m - 1) * NumSimPeriods + t - 1, 6] + a2prev

      if (snow == 1 && n1t == 0 && n2t == 0) {
        FakeData[(m - 1) * NumSimPeriods + t, 3] <- 1
      } else if (snow == 1 && n1t == 0 && n2t == 1) {
        FakeData[(m - 1) * NumSimPeriods + t, 3] <- 2
      } else if (snow == 1 && n1t == 1 && n2t == 0) {
        FakeData[(m - 1) * NumSimPeriods + t, 3] <- 3
      } else if (snow == 1 && n1t == 1 && n2t == 1) {
        FakeData[(m - 1) * NumSimPeriods + t, 3] <- 4
      } else if (snow == 2 && n1t == 0 && n2t == 0) {
        FakeData[(m - 1) * NumSimPeriods + t, 3] <- 5
      } else if (snow == 2 && n1t == 0 && n2t == 1) {
        FakeData[(m - 1) * NumSimPeriods + t, 3] <- 6
      } else if (snow == 2 && n1t == 1 && n2t == 0) {
        FakeData[(m - 1) * NumSimPeriods + t, 3] <- 7
      } else if (snow == 2 && n1t == 1 && n2t == 1) {
        FakeData[(m - 1) * NumSimPeriods + t, 3] <- 8
      }
    }

    # 上記の手続きで決まった状態変数を、以下で用いやすくするために、
    # sと名づける
    s <- FakeData[(m - 1) * NumSimPeriods + t, 3]

    # 以下では、各状態変数だった時に、どのようなデータを格納していく
    # かを逐次決定していく。例としてs=1 (G00)の時について説明を行う
    # が、s=2以降については省略する

    if (s == 1) {
      # s=1の場合、状態変数は1であり、各企業は出店していないという
      # ことなので、[ 1, 0, 0 ]をFakeDataの４ー６行目に格納する
      FakeData[(m - 1) * NumSimPeriods + t, 4:6] <- c(1, 0, 0)

      # 企業1はCCP1UpdatedMat(s,2)という確率で、現状維持をして、
      # 1-CCP1UpdatedMat(s,2)で新規出店を行うので、もしある乱数
      # がCCP1UpdatedMat(s,2)より大きければa1=1をアサインし、
      # そうでなければa1=0を生成された疑似データとして格納する
      if (RandomNumbers[m, t, 1] > CCP1UpdatedMat[s, 2]) {
        FakeData[(m - 1) * NumSimPeriods + t, 7] <- 1
      }

      # 企業2も企業1と同様、引いてきた一様乱数が
      # CCP2UpdatedMat(s,2)よりも大きければa2=1を、小さければ
      # a2=0を生成された疑似データとして格納する
      if (RandomNumbers[m, t, 2] > CCP2UpdatedMat[s, 2]) {
        FakeData[(m - 1) * NumSimPeriods + t, 8] <- 1
      }
    } else if (s == 2) {
      FakeData[(m - 1) * NumSimPeriods + t, 4:6] <- c(1, 0, 1)

      if (RandomNumbers[m, t, 1] > CCP1UpdatedMat[s, 2]) {
        FakeData[(m - 1) * NumSimPeriods + t, 7] <- 1
      }

      if (RandomNumbers[m, t, 2] > CCP2UpdatedMat[s, 2]) {
        FakeData[(m - 1) * NumSimPeriods + t, 8] <- -1
      }
    } else if (s == 3) {
      FakeData[(m - 1) * NumSimPeriods + t, 4:6] <- c(1, 1, 0)

      if (RandomNumbers[m, t, 1] > CCP1UpdatedMat[s, 2]) {
        FakeData[(m - 1) * NumSimPeriods + t, 7] <- -1
      }

      if (RandomNumbers[m, t, 2] > CCP2UpdatedMat[s, 2]) {
        FakeData[(m - 1) * NumSimPeriods + t, 8] <- 1
      }
    } else if (s == 4) {
      FakeData[(m - 1) * NumSimPeriods + t, 4:6] <- c(1, 1, 1)

      if (RandomNumbers[m, t, 1] > CCP1UpdatedMat[s, 2]) {
        FakeData[(m - 1) * NumSimPeriods + t, 7] <- -1
      }

      if (RandomNumbers[m, t, 2] > CCP2UpdatedMat[s, 2]) {
        FakeData[(m - 1) * NumSimPeriods + t, 8] <- -1
      }
    } else if (s == 5) {
      FakeData[(m - 1) * NumSimPeriods + t, 4:6] <- c(2, 0, 0)

      if (RandomNumbers[m, t, 1] > CCP1UpdatedMat[s, 2]) {
        FakeData[(m - 1) * NumSimPeriods + t, 7] <- 1
      }

      if (RandomNumbers[m, t, 2] > CCP2UpdatedMat[s, 2]) {
        FakeData[(m - 1) * NumSimPeriods + t, 8] <- 1
      }
    } else if (s == 6) {
      FakeData[(m - 1) * NumSimPeriods + t, 4:6] <- c(2, 0, 1)

      if (RandomNumbers[m, t, 1] > CCP1UpdatedMat[s, 2]) {
        FakeData[(m - 1) * NumSimPeriods + t, 7] <- 1
      }

      if (RandomNumbers[m, t, 2] > CCP2UpdatedMat[s, 2]) {
        FakeData[(m - 1) * NumSimPeriods + t, 8] <- -1
      }
    } else if (s == 7) {
      FakeData[(m - 1) * NumSimPeriods + t, 4:6] <- c(2, 1, 0)

      if (RandomNumbers[m, t, 1] > CCP1UpdatedMat[s, 2]) {
        FakeData[(m - 1) * NumSimPeriods + t, 7] <- -1
      }

      if (RandomNumbers[m, t, 2] > CCP2UpdatedMat[s, 2]) {
        FakeData[(m - 1) * NumSimPeriods + t, 8] <- 1
      }
    } else if (s == 8) {
      FakeData[(m - 1) * NumSimPeriods + t, 4:6] <- c(2, 1, 1)

      if (RandomNumbers[m, t, 1] > CCP1UpdatedMat[s, 2]) {
        FakeData[(m - 1) * NumSimPeriods + t, 7] <- -1
      }

      if (RandomNumbers[m, t, 2] > CCP2UpdatedMat[s, 2]) {
        FakeData[(m - 1) * NumSimPeriods + t, 8] <- -1
      }
    }
  }
}

# 4. Aguirregabiria and Mira (2007)の方法によるパラメターの推定----

# 紙面における点推定値を同じものを得るためには、
# Matlabで生成したFakedataを利用する必要がある
# なおこの場合でも、Bootstrapにより計算した標準誤差については乱数の影響で若干の差が出ることには留意されたい

# isUseMatlabData = 0のときは上で生成したFakedataを利用する
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
  if (check_convergence(newCCP, cbind(ccp1, ccp2))) break

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

print("推定されたCCP")
print(newCCP)
print("真のCCP")
print(cbind(CCP1Updated, CCP2Updated))

# 5. Bootstrapによる標準誤差の計算 ----
# マーケット単位でリサンプリングを行う。マーケットは500個
# Bootstrapのための乱数を固定
set.seed(2023)

# Bootstrap のリサンプリング回数
# 発散してしまうbootstrap sampleがあるため、多めに設定しておく
numBootSample <- 110

# 各Bootstrap sampleで用いるマーケットのインデックスを乱数から発生させる
# 市場がNumSimMarkets個であるため、1からNumSimMarketsの整数について、重複を許してNumSimMarkets個ドローする
bootindex <- matrix(sample(1:NumSimMarkets, NumSimMarkets * numBootSample, replace = TRUE),
  ncol = numBootSample
)

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



