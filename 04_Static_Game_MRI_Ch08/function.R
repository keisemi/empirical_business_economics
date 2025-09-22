# define functions for Ch08

# Bresnahan and Reiss (1991b)にしたがって目的関数を設定する関数 ----

## 以下で出てくるVやVVは、N_{max} = 6であれば、
# V = [
#   alpha_1  0  0  0  0  0
#   alpha_1  -alpha_2  0  0  0  0
#   alpha_1  -alpha_2  alpha_3  0  0  0 
#   alpha_1  -alpha_2  alpha_3  -alpha_4  0  0 
#   alpha_1  -alpha_2  alpha_3  -alpha_4  -alpha_5  0 
#   alpha_1  -alpha_2  alpha_3  -alpha_4  -alpha_5  -alpha_6
# ]
# VV = [
#   alpha_1 
#   alpha_1 - alpha_2 
#   alpha_1 - alpha_2 - alpha_3 
#   alpha_1 - alpha_2 - alpha_3 - alpha_4
#   alpha_1 - alpha_2 - alpha_3 - alpha_4 - alpha_5 
#   alpha_1 - alpha_2 - alpha_3 - alpha_4 - alpha_5 - \alpha_6 
# ]
#がM回（市場の個数分）繰り返される行列であり、病院がMRIスキャナーを購入した際の利潤を定義することに役立つ行列

obj <- function(params, 
                dataset, 
                N_max){
  
  # 本文に沿ってパラメターの名前をつけ直す
  # ただし、alpha2はマイナス1を乗じていることに注意
  alpha1 <- params[1]
  alpha2 <- -params[2:N_max]
  alpha <- c(alpha1, alpha2)
  gamma <- params[N_max+1]
  
  # 推定に用いる変数をデータから抽出する
  NumMRI <- dataset$NumMRI
  M <- nrow(dataset)
  pop <- as.matrix(dataset$Pop)
  pop <- matrix(pop, M, N_max)
  
  # 可変利潤を定義する際に便利な行列（上で説明済）を定義する
  V <- matrix(0, N_max, N_max)
  for (i in 1:length(alpha)){
    V[i:N_max, i] <- alpha[i]
  }
  
  # 可変利潤部分の行列（ただしPop_mを乗じる前）を定義する
  VV <- t(V %*% matrix(1, N_max, M))
  
  # 固定費用部分の行列を定義する
  F <- matrix(gamma, M, N_max)
  
  # 各市場mの企業の参入時の利潤を定義する
  pi <- pop*VV-F
  
  # 標準正規分布のCDFに各値を代入する
  phi <- pnorm(pi, mean=0, sd=1)
  
  # 端点（n=0 および n=N_max）に注意しながら、各企業数になる確率を計算する
  mat <- cbind(matrix(1, M, 1)-phi[,1], phi[,1:N_max-1]-phi[,2:N_max], phi[,N_max])
  
  # データで観測される企業数の確率を抽出する
  ml <- rep(0, M)
  for (i in 0:N_max){
    ml[NumMRI==i] <- log(mat[, i+1][NumMRI==i])
  }
  
  # 対数尤度が定義できない場合（確率が小さすぎて対数をとれない場合）の対処
  ml <- replace(ml, which(is.infinite(ml)), -10000)
  
  # 対数尤度の和をとることで目的関数の値を計算する
  val <- sum(ml)/M
  return(val)
}