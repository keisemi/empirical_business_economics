# 疑似データの作成

# Step 1: CCPの初期値の設定（および設定したCCPをベクトルから行列へ変換する)

# 以下ではCCPを全て0.5とする
CCP1 <- matrix(0.5, 8, 1)
CCP2 <- matrix(0.5, 8, 1)
# ただし、この行列は同じ確率である必要はなく、例えば以下のような行列でも構わない
# CCP1 <- matrix(c(.5, .55, .6, .65, .4, .6, .45, .55))
# CCP2 <- matrix(c(.5, .6, .55, .65, .4, .45, .6, .55))

# 以下の関数を用いて上で与えた初期値のベクトルを行列へと変換する。
# 例えば、CCP1Matは
# CCP1Mat = [ 0   0.5 0.5
# 0   0.5 0.5
# :    :   :
# 0.5  0.5  0
# 0.5  0.5  0  ] のような行列になる
CCP1Mat <- CCP1Transform(CCP1)
CCP2Mat <- CCP2Transform(CCP2)


# Step 2: Step 1 で与えられたCCPの初期値を基に、事前の価値関数を計算する

# 基本的には（４）式の計算を行いたいため、以下で必要なパーツを計算する

# まず、F^P^sigmaを求める。これは109ページの囲みで与えられているように、
# 遷移行列とCCP1とCCP2から作られる３つの行列の、要素ごとの積（アダマール
# 積）として計算できる
fPsigma <- fP(TransitionMat, CCP1, CCP2)

# パラメターの下で、pi1Psigma と pi2Psigma を計算する

pi1Psigma <- pi1PsigmaGen(pi1, CCP2Mat)
pi2Psigma <- pi2PsigmaGen(pi2, CCP1Mat)

# e^P_i(a_i,s)を計算する
eP1 <- eulergamma - CCP1LogTransform(CCP1)
eP2 <- eulergamma - CCP2LogTransform(CCP2)

# 上で求めたものを（４）式に代入して、事前の価値関数を求める
ExanteV1 <- as.matrix(matlib::inv(diag(8) - beta * fPsigma) %*%
  ((CCP1Mat * (pi1Psigma + eP1)) %*% matrix(c(1, 1, 1))))
ExanteV2 <- matlib::inv(diag(8) - beta * fPsigma) %*%
  ((CCP2Mat * (pi2Psigma + eP2)) %*% matrix(c(1, 1, 1)))
# (pij + ePj) が何らかの値をとっていたとしても、CCP1MatやCCP2Matを
# かけることで、取ることのできない選択肢にはゼロが振られることになり、
# 分母の和は正しく計算できていることがわかる。例えば、
# CCP1Mat.*(pi1Psigma+eP1) を表示させてみると、実際にいくつかの
# 要素がゼロになっていることが確認できる。

## Step 3: 事前の価値関数を基に、CCPを計算する ----

# 各企業の行動が与えられたときの遷移行列を求める。これは、本文では
# f^P*(s'|s,a_i)として表現されている行列である。fP_a1/2はリスト形式になっているが
# fP_a1/2[[1]] がa1/2=-1であるときの、
# fP_a1/2[[2]] がa1/2=0 であるときの、
# fP_a1/2[[3]] がa1/2=1 であるときの遷移行列に相当します

fP_a1 <- fP_a1given(TransitionMat, CCP2)
fP_a2 <- fP_a2given(TransitionMat, CCP1)

# 以下では、求められた事前の価値関数をもとにCCPを更新する。最終的には、
# 計算の簡便化のため更新されたCCPを行列の形（CCP1UpdatedMat）とベクトル
# の形（CCP1Updated）の両方で求めておく。この計算は一見複雑だが、需要
# 関数の推定の時と似た計算であるため、説明は省略する。最初のブロックは企
# 業1について、次のブロックは企業2について更新されたCCPを求めている。
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
# ただし、 Step 2と同一の作業であるため、説明は省略する。

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
      # で、今期の景気の状態を4行目に格納する。
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
      # FakeDataの3列目に格納していく。
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
    # が、s=2以降については省略する。

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
