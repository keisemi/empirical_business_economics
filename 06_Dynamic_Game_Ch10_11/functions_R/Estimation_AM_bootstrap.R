Estimation_AM_bootstrap <- function(FakeData, beta) {
  # Aguirregabiria and Miraの推定方法による、ブートストラップを行う
  # FakeDataは、以下の推定に用いるデータ。ブートストラップにおいては、元のデータからリサンプリングしたデータになっている

  # Step 1: Estimating CCP and transition probabilities
  EstimatedCCP1 <- matrix(rep(0, 8), ncol = 1)
  EstimatedCCP2 <- matrix(rep(0, 8), ncol = 1)

  for (s in 1:8) {
    SubData <- FakeData[FakeData[, 3] == s, ]
    Subseta1iszero <- SubData[SubData[, 7] == 0, ]
    Subseta2iszero <- SubData[SubData[, 8] == 0, ]
    EstimatedCCP1[s, 1] <- dim(Subseta1iszero)[1] / dim(SubData)[1]
    EstimatedCCP2[s, 1] <- dim(Subseta2iszero)[1] / dim(SubData)[1]
  }

  EstimatedTransition <- matrix(rep(0, 4), ncol = 2)

  DataL1 <- matrix(c(0, FakeData[1:dim(FakeData)[1] - 1, 4]), ncol = 1)
  SubData <- cbind(FakeData, DataL1)
  SubData <- SubData[SubData[, 2] != 1, ]

  for (z in 1:2) {
    SubDataZ <- SubData[SubData[, 4] == z, ]
    SubDataZnext <- SubDataZ[SubDataZ[, 9] == z, ]
    EstimatedTransition[z, z] <- dim(SubDataZnext)[1] / dim(SubDataZ)[1]
  }

  EstimatedTransition[1, 2] <- 1 - EstimatedTransition[1, 1]
  EstimatedTransition[2, 1] <- 1 - EstimatedTransition[2, 2]

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

  # Step 3: AM estimator
  # iteration number
  i <- 1

  # データから観察された行動
  Actions <- FakeData[, c(7:8, 3)]

  # パラメターの初期値
  initial <- c(InitialParameters[1:4], InitialParameters[6])

  # Stack防止のためのゆらぎの乱数を準備する。
  set.seed(123456)
  random_draw <- runif(8*2*100, min = 0.75, max = 1.25)
  
  # 推定
  while (i > 0) {
    # 時間がかかりすぎる場合に強制的にスキップする
    if (i == 100) {
      print("CCP did not converge")
      return(NULL)
    }

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

    # CCPが収束していないか確認する
    if (check_convergence(newCCP, cbind(ccp1, ccp2), tol = 1e-6)) break

    # CCPが収束していない場合、新しいCCPでもとのCCPを更新する
    ccp1 <- newCCP[, 1]
    ccp2 <- newCCP[, 2]
    
    # Iterationが10の倍数のとき、CCPをRandomに揺らがせる。(Stack防止のため)
    if (mod(i,10) == 0){
      ccp1 <- ccp1*random_draw[ ((i/10 - 1)*8 + 1):((i/10 - 1)*8 + 8) ]
      ccp2 <- ccp2*random_draw[ ((i/10 - 1)*8 + 9):((i/10 - 1)*8 + 16) ]
    }

    # iを更新
    i <- i + 1
  }

  # 結果を出力
  output <- list(
    cbind(EstimatedCCP1, EstimatedCCP2),
    EstimatedTransition,
    sol$par
  )

  return(output)
}
