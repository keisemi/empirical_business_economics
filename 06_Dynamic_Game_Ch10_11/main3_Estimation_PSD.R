# 経済セミナー連載「実証ビジネス・エコノミクス」
# 第10回「出店戦略モデルの性質を見極める」
# 第11回「市場のダイナミクスから利益構造を見通す」
# 上武康亮・遠山祐太・若森直樹・渡辺安虎
# 最終更新：2023年5月20日

# 謝辞：MatlabからRコードへの書き換えに際して、
# 箕輪創太さん(東京大学大学院経済学研究科)及び山田渉生さん(早稲田大学政治経済学部)には多大なるご尽力を頂きました。
# この場を借りて御礼申し上げます。

# 1. はじめに----

# このRファイルは経済セミナー「実証ビジネス・エコノミクス」
# 第10回「出店戦略モデルの性質を見極める」
# 第11回「市場のダイナミクスから利益構造を見通す」
# の分析を再現します
# 乱数を利用しているため、matlab用のコードを利用した場合の数値
# 及び連載に掲載された数値と異なる結果が出力されます。

# 本連載の内容、およびサンプルコード等の資料は、情報の提供のみを目的として
# いますので、運用につきましては十分にご確認をいただき、お客様ご自身の責任
# とご判断によって行ってください。これらの情報の運用結果により損害が生じた
# 場合でも、日本評論社および著者はいかなる責任を負うことはできませんので、
# ご留意ください。

# 全体の流れは以下のようになっています
# 1. はじめに
# 2. 下準備
# 3. 疑似データの生成 (第10回の内容に相当)
# 4. Pesendorer and Schmidt-Dengler (以下P-SD) の方法によるパラメタの推定 (以下第11回) 
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
                         pattern="*.R$", full.names=TRUE, 
                         ignore.case=TRUE)
sapply(functionlt, source)

# 割引因子
beta <- 0.8

# 定数
eulergamma <- 0.5772

# 景気の遷移行列
TransitionMat <- matrix(c(0.7,0.3, 0.4, 0.6), nrow = 2, byrow = TRUE)
# デフォルトではbyrowがFALSEになっており、列方向に各要素が代入されます。

# パラメターの設定
Parameters <- matrix(c(
  0.3, # 企業1のベース利潤
  0.2, # 企業2のベース利潤
  -0.27, # 顧客収奪効果
  0.45, # 景気が良い時の追加的利潤
  -0.15, # 退出のためのコスト
  -2.10 # 参入のためのコスト
))
# matrix(c())はncolやnrowを指定しない場合、n行1列の行列を作成します。 

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
  Parameters[6]  # 企業2の出店のためのコスト
))

# 後々の計算を簡略化するための行列として、CCP1AdjusterとCCP2Adjusterを
# 定義する。これらは、以下のように、各状態変数で選択することができる選択肢には
# 1が、選択不可能な選択肢には0が割り振られている行列である
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

# 真のパラメターの下で、pi1とpi2を計算する。各々は8*3の行列になっていて、縦が
# 状態変数（G00,G01,G10,...,B11）を意味し、横がa_i = -1, 0, 1 という行動を
# 意味し、それぞれの状態変数で行動を選んだ時の利潤が行列として与えられている。
# 一般的には、企業1の利潤関数はpi_1(a_1,a_2,s)のように、自身の行動a_1だけで
# なく、企業2の行動a_2にも依存するが、今回の特定化ではお互いの企業の利潤は、
# ライバル企業の行動の影響を受けない関数形になっていることに注意。

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

# 以下の関数を用いて上で与えた初期値のベクトルを行列へと変換する。
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
# Rでは自動的に行列の次元が設定されます。
eP1 <- eulergamma - CCP1LogTransform(CCP1)
eP2 <- eulergamma - CCP2LogTransform(CCP2)
# 少々長くはなりますが、以下のように書くことも可能です。 
# eP1 <- eulergamma * matrix(rep(1, 8*3), ncol = 3) - CCP1LogTransform(CCP1)
# eP2 <- eulergamma * matrix(rep(1, 8*3), ncol = 3) - CCP2LogTransform(CCP2)

# 上で求めたものを（４）式に代入して、事前の価値関数を求める
ExanteV1 <- as.matrix(matlib::inv(diag(8) - beta*fPsigma) %*% 
                        ((CCP1Mat * (pi1Psigma+eP1)) %*% matrix(c(1,1,1))))
ExanteV2 <- matlib::inv(diag(8) - beta*fPsigma) %*% 
  ((CCP2Mat * (pi2Psigma+eP2)) %*% matrix(c(1,1,1)))
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
#fP_a1はlistであるため、この様にかっこをつけなければならない。
#[]を一つにするとclassがlistのままで行列の積が計算できない。
NewSigmaSeed1  <- (pi1Psigma + beta*c(fP_a1[[1]]%*% ExanteV1, fP_a1[[2]]%*%ExanteV1 , fP_a1[[3]]%*%ExanteV1 ))*CCP1Adjuster
NewSigmaDeno1  <- apply(exp(NewSigmaSeed1), 1, sum) - rep(1,8)
NewSigmaDeno1M <- rep(NewSigmaDeno1,3)
NewSigma1      <- exp(NewSigmaSeed1)/NewSigmaDeno1M
CCP1UpdatedMat <- NewSigma1*CCP1Adjuster
CCP1Updated    <- CCP1UpdatedMat[,2]

NewSigmaSeed2  <-  (pi2Psigma + beta*c(fP_a2[[1]]%*% ExanteV2, fP_a2[[2]]%*%ExanteV2 , fP_a2[[3]]%*%ExanteV2 ))*CCP2Adjuster
NewSigmaDeno2  <- apply(exp(NewSigmaSeed2), 1, sum) - rep(1,8)
NewSigmaDeno2M <- rep(NewSigmaDeno2,3)
NewSigma2      <- exp(NewSigmaSeed2)/NewSigmaDeno2M
CCP2UpdatedMat <- NewSigma2*CCP2Adjuster
CCP2Updated    <- CCP2UpdatedMat[,2]


## Step 4: アップデートされたCCPを基に、再び事前の価値関数を計算する----
# ただし、 Step 2と同一の作業であるため、説明は省略する。

fPsigma <- fP(TransitionMat, CCP1Updated, CCP2Updated)

pi1Psigma <- pi1PsigmaGen(pi1, CCP2UpdatedMat)
pi2Psigma <- pi2PsigmaGen(pi2, CCP1UpdatedMat)

eP1 <- eulergamma - CCP1LogTransform(CCP1Updated)
eP2 <- eulergamma - CCP2LogTransform(CCP2Updated)

ExanteV1Updated <- inv((diag(8)-beta*fPsigma)) %*% apply(CCP1UpdatedMat*(pi1Psigma+eP1),1,sum)
ExanteV2Updated <- inv((diag(8)-beta*fPsigma)) %*% apply(CCP2UpdatedMat*(pi2Psigma+eP2),1,sum)

## Step 5: Step 3 と Step 4 を事前の価値関数が一致するまで繰り返す----

# 前の事前の価値関数と更新された事前の価値関数の差を定義する
DiffExanteV <- apply((ExanteV1Updated - ExanteV1)^2 +(ExanteV2Updated - ExanteV2)^2 ,2 ,sum)

# 上で定義した差が1.0e-12よりも小さくなるまで、繰り返す
while ( DiffExanteV > 1.0e-12 ){
  
  # 更新されたCCPと事前の価値関数を、CCPとExanteVとして置き換える
  CCP1 <- CCP1Updated
  CCP2 <- CCP2Updated
  
  ExanteV1 <- ExanteV1Updated
  ExanteV2 <- ExanteV2Updated
  
  # Step 3 を再度実行する（上と全く一緒のため、説明は省略する）
  fP_a1 <- fP_a1given(TransitionMat, CCP2)
  fP_a2 <- fP_a2given(TransitionMat, CCP1)
  
  NewSigmaSeed1  <- (pi1Psigma + beta*c(fP_a1[[1]]%*% ExanteV1, fP_a1[[2]]%*%ExanteV1 , fP_a1[[3]]%*%ExanteV1 ))*CCP1Adjuster
  NewSigmaDeno1  <- apply(exp(NewSigmaSeed1), 1, sum) - rep(1,8)
  NewSigmaDeno1M <- rep(NewSigmaDeno1,3)
  NewSigma1      <- exp(NewSigmaSeed1)/NewSigmaDeno1M
  CCP1UpdatedMat <- NewSigma1*CCP1Adjuster
  CCP1Updated    <- CCP1UpdatedMat[,2]
  
  NewSigmaSeed2  <-  (pi2Psigma + beta*c(fP_a2[[1]]%*% ExanteV2, fP_a2[[2]]%*%ExanteV2 , fP_a2[[3]]%*%ExanteV2 ))*CCP2Adjuster
  NewSigmaDeno2  <- apply(exp(NewSigmaSeed2), 1, sum) - rep(1,8)
  NewSigmaDeno2M <- rep(NewSigmaDeno2,3)
  NewSigma2      <- exp(NewSigmaSeed2)/NewSigmaDeno2M
  CCP2UpdatedMat <- NewSigma2*CCP2Adjuster
  CCP2Updated    <- CCP2UpdatedMat[,2]
  
  # Step 4 を再度実行する（上と全く一緒のため、説明は省略する）
  fPsigma <- fP(TransitionMat, CCP1Updated, CCP2Updated)
  
  pi1Psigma <- pi1PsigmaGen(pi1, CCP2UpdatedMat)
  pi2Psigma <- pi2PsigmaGen(pi2, CCP1UpdatedMat)
  
  eP1 <- eulergamma - CCP1LogTransform(CCP1Updated)
  eP2 <- eulergamma - CCP2LogTransform(CCP2Updated)
  
  ExanteV1Updated <- inv((diag(8)-beta*fPsigma)) %*% apply(CCP1UpdatedMat*(pi1Psigma+eP1),1,sum)
  ExanteV2Updated <- inv((diag(8)-beta*fPsigma)) %*% apply(CCP2UpdatedMat*(pi2Psigma+eP2),1,sum)
  
  # Step 5 の冒頭部分の差分の計算を再度実行する
  DiffExanteV <- apply((ExanteV1Updated - ExanteV1)^2 +(ExanteV2Updated - ExanteV2)^2 ,2 ,sum)
}

# 上記のwhileループを実行することで、均衡のCCPを求めることができる。
# ここで、確認のため、どのようなCCPになっているかを表示させる
print( matrix(c(CCP1UpdatedMat ,CCP2UpdatedMat), ncol = 6)) # 均衡でのCCP


## Step 6: 均衡におけるCCPをもとに、疑似データをシミュレートする----

# まず、シミュレーションの設定を行う。
# 乱数の発生が必要なため、乱数発生のためのシードを設定する
set.seed(2023)

# 次に、何個の市場で何期のデータを生成するかを設定し、企業数は2とする
NumSimMarkets <- 500 #市場の数
NumSimPeriods <- 50  #シミュレーションの期間
NumSimFirms   <- 2   #企業の数

# さらに、各市場の1期の状態変数を、一様分布に従ってランダムに生成する
InitialState  <- sample(1:8, NumSimMarkets, replace = TRUE)  

# 企業1と企業2がどのような行動をとるのか、及び、次期の状態変数を決める際
# に必要となる一様乱数を RandomNumbers という行列に格納する
RandomNumbers <- array(runif(NumSimMarkets*NumSimPeriods*(NumSimFirms+1), min=0, max=1)
                       ,dim = c(NumSimMarkets, NumSimPeriods, NumSimFirms+1))

# 最後に疑似データを格納するために、FakeDataという行列を準備する
FakeData <- matrix(rep(0,NumSimMarkets*NumSimPeriods*8),ncol=8)
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
    FakeData[(m-1)*NumSimPeriods+t,1] <- m
    FakeData[(m-1)*NumSimPeriods+t,2] <- t
    
    # 1期目の状態変数は先の通り乱数で発生させているが、2期目以降は
    # 内生的に決まることから、ここでは1期目と2期目以降は別々のプロ
    # グラムとなっている
    if (t == 1) {
      
      # 先に決めたInitialState(m)から、該当する市場の初期状態を
      # FakeDataの3列目に格納する
      FakeData[(m-1)*NumSimPeriods+1,3] <- InitialState[m]
      
      # 格納した初期状態に応じて、景気の状態がどうなっているかを
      # 判断して、FakeDataの4列目に格納する
      if (InitialState[m] >= 1 && InitialState[m] <= 4){
        FakeData[(m-1)*NumSimPeriods+1,4] <- 1
      } else{
        FakeData[(m-1)*NumSimPeriods+1,4] <- 2
      }
    } else{
      
      # t=2,...,NumSimPeriodsの場合, まずは前期の状態変数を一行
      # 前を参照にして sprev (previous state)として持ってくる。
      # さらに、今期の状態変数は前期の各企業の行動にも依存するので、
      # 前期の各企業の行動をa1prev, a2prevとしてもってくる
      sprev  <-  FakeData[(m-1)*NumSimPeriods+t-1,3]
      a1prev <-  FakeData[(m-1)*NumSimPeriods+t-1,7]
      a2prev <-  FakeData[(m-1)*NumSimPeriods+t-1,8]
      
      # 以下では、まず前期の状態変数から、景気の状態についての情報
      # 取り出し、一様乱数を用いて今期の景気の状態を決める。その上
      # で、今期の景気の状態を4行目に格納する。
      if  (sprev >= 1 && sprev <= 4) {
        if (RandomNumbers[m,t,3] < TransitionMat[1,1]){
          FakeData[(m-1)*NumSimPeriods+t,4] <- 1
        } else{
          FakeData[(m-1)*NumSimPeriods+t,4] <- 2
        }
      } else{
        if (RandomNumbers[m,t,3] < TransitionMat[2,2]){
          FakeData[(m-1)*NumSimPeriods+t,4] <- 2
        } else{
          FakeData[(m-1)*NumSimPeriods+t,4] <- 1
        }
      } 
      
      # さらに、今期の景気の状態（snow）と、今期の各企業の出店数
      # （n1tとn2t）を基に、今期の状態変数が何かを判断して、
      # FakeDataの3列目に格納していく。
      snow <- FakeData[(m-1)*NumSimPeriods+t,4]
      n1t  <- FakeData[(m-1)*NumSimPeriods+t-1,5]+a1prev
      n2t  <- FakeData[(m-1)*NumSimPeriods+t-1,6]+a2prev
      
      if (snow == 1 && n1t==0 && n2t==0){
        FakeData[(m-1)*NumSimPeriods+t,3] <- 1
      } else if (snow == 1 && n1t==0 && n2t==1){
        FakeData[(m-1)*NumSimPeriods+t,3] <- 2
      } else if (snow == 1 && n1t==1 && n2t==0){
        FakeData[(m-1)*NumSimPeriods+t,3] <- 3
      } else if (snow == 1 && n1t==1 && n2t==1){
        FakeData[(m-1)*NumSimPeriods+t,3] <- 4
      } else if (snow == 2 && n1t==0 && n2t==0){
        FakeData[(m-1)*NumSimPeriods+t,3] <- 5
      } else if (snow == 2 && n1t==0 && n2t==1){
        FakeData[(m-1)*NumSimPeriods+t,3] <- 6
      } else if (snow == 2 && n1t==1 && n2t==0){
        FakeData[(m-1)*NumSimPeriods+t,3] <- 7
      } else if (snow == 2 && n1t==1 && n2t==1){
        FakeData[(m-1)*NumSimPeriods+t,3] <- 8
      }  
    }
    
    # 上記の手続きで決まった状態変数を、以下で用いやすくするために、
    # sと名づける
    s <- FakeData[(m-1)*NumSimPeriods+t,3]
    
    # 以下では、各状態変数だった時に、どのようなデータを格納していく
    # かを逐次決定していく。例としてs=1 (G00)の時について説明を行う
    # が、s=2以降については省略する。
    
    if (s==1) {
      
      # s=1の場合、状態変数は1であり、各企業は出店していないという
      # ことなので、[ 1, 0, 0 ]をFakeDataの４ー６行目に格納する
      FakeData[(m-1)*NumSimPeriods+t,4:6] <- c(1, 0, 0)
      
      # 企業1はCCP1UpdatedMat(s,2)という確率で、現状維持をして、
      # 1-CCP1UpdatedMat(s,2)で新規出店を行うので、もしある乱数
      # がCCP1UpdatedMat(s,2)より大きければa1=1をアサインし、
      # そうでなければa1=0を生成された疑似データとして格納する
      if (RandomNumbers[m,t,1] > CCP1UpdatedMat[s,2]){
        FakeData[(m-1)*NumSimPeriods+t,7] <- 1
      }  
      
      # 企業2も企業1と同様、引いてきた一様乱数が
      # CCP2UpdatedMat(s,2)よりも大きければa2=1を、小さければ
      # a2=0を生成された疑似データとして格納する
      if (RandomNumbers[m,t,2] > CCP2UpdatedMat[s,2]){
        FakeData[(m-1)*NumSimPeriods+t,8] <- 1
      }  
    } else if (s==2){
      
      FakeData[(m-1)*NumSimPeriods+t,4:6] <-c(1, 0, 1)
      
      if (RandomNumbers[m,t,1] > CCP1UpdatedMat[s,2]){
        FakeData[(m-1)*NumSimPeriods+t,7] <- 1
      }  
      
      if (RandomNumbers[m,t,2] > CCP2UpdatedMat[s,2]){
        FakeData[(m-1)*NumSimPeriods+t,8] <- -1
      }    
    } else if (s==3){
      FakeData[(m-1)*NumSimPeriods+t,4:6] <-  c(1, 1, 0)
      
      if (RandomNumbers[m,t,1] > CCP1UpdatedMat[s,2]){
        FakeData[(m-1)*NumSimPeriods+t,7] <- -1
      }
      
      if (RandomNumbers[m,t,2] > CCP2UpdatedMat[s,2]){
        FakeData[(m-1)*NumSimPeriods+t,8] <- 1
      }
    } else if (s==4){
      
      FakeData[(m-1)*NumSimPeriods+t,4:6] <- c(1, 1, 1)
      
      if (RandomNumbers[m,t,1] > CCP1UpdatedMat[s,2]){
        FakeData[(m-1)*NumSimPeriods+t,7] <- -1
      }
      
      if (RandomNumbers[m,t,2] > CCP2UpdatedMat[s,2]){
        FakeData[(m-1)*NumSimPeriods+t,8] <- -1
      }
    } else if (s==5){
      
      FakeData[(m-1)*NumSimPeriods+t,4:6] <- c(2, 0, 0)
      
      if (RandomNumbers[m,t,1] > CCP1UpdatedMat[s,2]){
        FakeData[(m-1)*NumSimPeriods+t,7] <- 1
      }
      
      if (RandomNumbers[m,t,2] > CCP2UpdatedMat[s,2]){
        FakeData[(m-1)*NumSimPeriods+t,8] <- 1
      }
    } else if (s==6){
      FakeData[(m-1)*NumSimPeriods+t,4:6] <- c(2, 0, 1)
      
      if (RandomNumbers[m,t,1] > CCP1UpdatedMat[s,2]){
        FakeData[(m-1)*NumSimPeriods+t,7] <- 1
      }
      
      if (RandomNumbers[m,t,2] > CCP2UpdatedMat[s,2]){
        FakeData[(m-1)*NumSimPeriods+t,8] <- -1
      }
    } else if (s==7){
      FakeData[(m-1)*NumSimPeriods+t,4:6] <- c(2, 1, 0)
      
      if (RandomNumbers[m,t,1] > CCP1UpdatedMat[s,2]){
        FakeData[(m-1)*NumSimPeriods+t,7] <- -1
      }
      
      if (RandomNumbers[m,t,2] > CCP2UpdatedMat[s,2]){
        FakeData[(m-1)*NumSimPeriods+t,8] <- 1
      }
    } else if (s==8){
      
      FakeData[(m-1)*NumSimPeriods+t,4:6] <- c(2, 1, 1)
      
      if (RandomNumbers[m,t,1] > CCP1UpdatedMat[s,2]){
        FakeData[(m-1)*NumSimPeriods+t,7] <- -1
      }
      
      if (RandomNumbers[m,t,2] > CCP2UpdatedMat[s,2]){
        FakeData[(m-1)*NumSimPeriods+t,8] <- -1
      }
    }
    
  }
  
}


# 4. パラメタの推定----

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

#Transition CCP
Summary_Transition <- matrix( c(c(0.7, 0.6), diag(EstimatedTransition) , apply(bootresult_transition,1,sd)) , ncol = 3)
colnames(Summary_Transition) <- c("Transition Probability (GG and BB): True", "Estimated", "SE") 
print(Summary_Transition)

#Payoff parameter
true <- initial
normalized_param <- c( 0.3 - (1-beta)/beta*(-0.15), 0.2 - (1-beta)/beta*(-0.15),-0.27, 0.45, -2.1 + (-0.15))
Summary_Payoff_param <- matrix(c( true, normalized_param , result_pick, apply(bootresult_payoff,1,sd)) ,ncol = 4)
colnames(Summary_Payoff_param) <- c("Payoff parameter: True", "Normalized true", "Estimated", "SE")
print(Summary_Payoff_param)
