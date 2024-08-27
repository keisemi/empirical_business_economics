obj_lik <- function(param, GivenCCP1, GivenCCP2, TransitionProb, beta, Data){
  # パラメターの評価をするための対数尤度関数
  # input
  ## param: パラメター
  ## GivenCCP1: 企業1のCCP
  ## GivenCCP2: 企業2のCCP
  ## TransitionProb: 遷移確率行列
  ## beta: 割引率
  ## Data: データ
  # output
  ## ll: 対数尤度
  
  # パラメターと与えられたCCPから新たなCCPを計算
  CCPs <- CCP_to_Value_to_prediction(param, GivenCCP1, GivenCCP2, TransitionProb, beta)
  
  # 予測されたCCPとDataから観察された行動から尤度を計算
  ll <- sum((Data[, 1] == 0 & Data[, 3] == 1) * log(CCPs[1, 1]), # stateが0-8かつ企業1の選択が0であることの対数尤度
            (Data[, 1] == 0 & Data[, 3] == 2) * log(CCPs[2, 1]),
            (Data[, 1] == 0 & Data[, 3] == 3) * log(CCPs[3, 1]),
            (Data[, 1] == 0 & Data[, 3] == 4) * log(CCPs[4, 1]),
            (Data[, 1] == 0 & Data[, 3] == 5) * log(CCPs[5, 1]),
            (Data[, 1] == 0 & Data[, 3] == 6) * log(CCPs[6, 1]),
            (Data[, 1] == 0 & Data[, 3] == 7) * log(CCPs[7, 1]),
            (Data[, 1] == 0 & Data[, 3] == 8) * log(CCPs[8, 1]),
            (Data[, 2] == 0 & Data[, 3] == 1) * log(CCPs[1, 2]), # stateが0-8かつ企業2の選択が0であることの対数尤度
            (Data[, 2] == 0 & Data[, 3] == 2) * log(CCPs[2, 2]),
            (Data[, 2] == 0 & Data[, 3] == 3) * log(CCPs[3, 2]),
            (Data[, 2] == 0 & Data[, 3] == 4) * log(CCPs[4, 2]),
            (Data[, 2] == 0 & Data[, 3] == 5) * log(CCPs[5, 2]),
            (Data[, 2] == 0 & Data[, 3] == 6) * log(CCPs[6, 2]),
            (Data[, 2] == 0 & Data[, 3] == 7) * log(CCPs[7, 2]),
            (Data[, 2] == 0 & Data[, 3] == 8) * log(CCPs[8, 2]),
            (Data[, 1] != 0 & Data[, 3] == 1) * log(1 - CCPs[1, 1]), # stateが0-8かつ企業1の選択が0以外であることの対数尤度
            (Data[, 1] != 0 & Data[, 3] == 2) * log(1 - CCPs[2, 1]),
            (Data[, 1] != 0 & Data[, 3] == 3) * log(1 - CCPs[3, 1]),
            (Data[, 1] != 0 & Data[, 3] == 4) * log(1 - CCPs[4, 1]),
            (Data[, 1] != 0 & Data[, 3] == 5) * log(1 - CCPs[5, 1]),
            (Data[, 1] != 0 & Data[, 3] == 6) * log(1 - CCPs[6, 1]),
            (Data[, 1] != 0 & Data[, 3] == 7) * log(1 - CCPs[7, 1]),
            (Data[, 1] != 0 & Data[, 3] == 8) * log(1 - CCPs[8, 1]),
            (Data[, 2] != 0 & Data[, 3] == 1) * log(1 - CCPs[1, 2]), # stateが0-8かつ企業2の選択が0以外であることの対数尤度
            (Data[, 2] != 0 & Data[, 3] == 2) * log(1 - CCPs[2, 2]),
            (Data[, 2] != 0 & Data[, 3] == 3) * log(1 - CCPs[3, 2]),
            (Data[, 2] != 0 & Data[, 3] == 4) * log(1 - CCPs[4, 2]),
            (Data[, 2] != 0 & Data[, 3] == 5) * log(1 - CCPs[5, 2]),
            (Data[, 2] != 0 & Data[, 3] == 6) * log(1 - CCPs[6, 2]),
            (Data[, 2] != 0 & Data[, 3] == 7) * log(1 - CCPs[7, 2]),
            (Data[, 2] != 0 & Data[, 3] == 8) * log(1 - CCPs[8, 2])
            )
  return(ll)
}