# 割引因子
beta <- 0.8

# 定数
eulergamma <- 0.5772

# 景気の遷移行列
TransitionMat <- matrix(c(0.7, 0.3, 0.4, 0.6), nrow = 2, byrow = TRUE)
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
# ライバル企業の行動の影響を受けない関数形になっていることに注意。
pi1 <- pi1gen(TrueParameterValues) * CCP1Adjuster
pi2 <- pi2gen(TrueParameterValues) * CCP2Adjuster
