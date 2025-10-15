
# パラメタを所与として走行距離の遷移行列を出力する関数 
gen_mileage_trans <- function(kappa) {
  # 走行距離が1,2段階上がる確率をパラメタとする
  kappa_1 <- kappa[1]
  kappa_2 <- kappa[2]
  # 購買しなかった場合の遷移行列を作成
  mileage_trans_mat_hat_not_buy <-
    matrix(0, ncol = num_mileage_states, nrow = num_mileage_states)
  for (i in 1:num_mileage_states) {
    for (j in 1:num_mileage_states) {
      if (i == j){
        mileage_trans_mat_hat_not_buy[i,j] <- 1 - kappa_1 - kappa_2
      } else if (i == j - 1) {
        mileage_trans_mat_hat_not_buy[i,j] <- kappa_1
      } else if (i == j - 2){
        mileage_trans_mat_hat_not_buy[i,j] <- kappa_2
      }
    }
  }
  mileage_trans_mat_hat_not_buy[num_mileage_states - 1, num_mileage_states] <- 
    kappa_1 + kappa_2
  mileage_trans_mat_hat_not_buy[num_mileage_states, num_mileage_states] <- 1
  # 購買した場合の遷移行列を作成
  # 購入した期では m=0 となるため次の期のmileageはそこから決まることに注意
  mileage_trans_mat_hat_buy <-
    matrix(1, nrow = num_mileage_states, ncol = 1) %*%
    mileage_trans_mat_hat_not_buy[1,]
  # 3次元のarrayとして出力
  return(array(c(mileage_trans_mat_hat_not_buy, 
                 mileage_trans_mat_hat_buy),
               dim = c(num_mileage_states, num_mileage_states, num_choice)))
}


# パラメタを所与として価格の遷移行列を出力する関数を作成 
gen_price_trans <- function(lambda){
  
  lambda_11 <- 1 - lambda[1] - lambda[2] - lambda[3] - lambda[4] - lambda[5]
  lambda_22 <- 1 - lambda[6] - lambda[7] - lambda[8] - lambda[9] - lambda[10]
  lambda_33 <- 1 - lambda[11] - lambda[12] - lambda[13] - lambda[14] - lambda[15]
  lambda_44 <- 1 - lambda[16] - lambda[17] - lambda[18] - lambda[19] - lambda[20]
  lambda_55 <- 1 - lambda[21] - lambda[22] - lambda[23] - lambda[24] - lambda[25]
  lambda_66 <- 1 - lambda[26] - lambda[27] - lambda[28] - lambda[29] - lambda[30]
  
  price_trans_mat_hat <- 
    c(lambda_11, lambda[1], lambda[2], lambda[3], lambda[4], lambda[5],
      lambda[6], lambda_22, lambda[7], lambda[8], lambda[9], lambda[10],
      lambda[11], lambda[12], lambda_33, lambda[13], lambda[14], lambda[15],
      lambda[16], lambda[17], lambda[18], lambda_44, lambda[19], lambda[20],
      lambda[21], lambda[22], lambda[23], lambda[24], lambda_55, lambda[25],
      lambda[26], lambda[27], lambda[28], lambda[29], lambda[30], lambda_66) %>% 
    matrix(ncol = num_price_states, nrow = num_price_states, byrow=T)
  
  return(price_trans_mat_hat)
}


# 状態変数、コントロール変数毎の今期の効用を返す関数
flow_utility <- function(theta,
                         state_df){
  theta_c <- theta[1]
  theta_p <- theta[2]
  
  U <- cbind(
      # その期における車を購入しない場合の効用
      U_not_buy = - theta_c * state_df$mileage, 
      # その期における車を購入する場合の効用
      U_buy = - theta_p * state_df$price
    ) 
  return(U)
}


# 価値関数反復法における縮小写像を定義する関数
contraction <- function(theta,
                        beta, 
                        trans_mat, 
                        state_df) {
    # 価値関数の初期値
    num_states <- nrow(state_df)
    V_old <- matrix(0, nrow = num_states, ncol = 1)
    
    # パラメタより今期の効用を計算
    U <- flow_utility(theta, state_df)
    
    # 価値関数の差の初期値
    diff <- 1000
    
    # 縮小写像の誤差範囲
    tol_level <- 1.0e-12
    
    while (diff > tol_level) {
      # 価値関数を計算
      V_new <- log(rowSums(exp(U + beta * cbind(trans_mat$not_buy %*% V_old,
                                                trans_mat$buy %*% V_old)))) + Euler_const
      # 価値関数の更新による差を評価
      diff <- max(abs(V_new - V_old))
      # 価値関数を次のループに渡す(価値関数の更新)
      V_old <- V_new
    }
    # 得られた事前の価値関数を出力
    return(V_old)
  }



# 各消費者についてデータを生成する関数を作成
generate_data <- function(df, 
                          V_CS, 
                          state_df, 
                          price_dist_steady) {
  
  # Step 1: 各消費者について、初期のstate_idを決める
  # 価格は定常分布に従うとし、走行距離は0とする
  
  # 価格の定常分布の累積値を計算
  price_dist_steady_cumsum <- cumsum(price_dist_steady)
  
  # 一様分布から生成した乱数が、定常分布の累積値を
  # 初めて下回ったところを1期の状態変数とする
  price_id_consumer <- 0
  exceed_trans_prob_price <- TRUE
  
  while(exceed_trans_prob_price) {
    price_id_consumer <- price_id_consumer + 1
    exceed_trans_prob_price <- 
      (df$eps_price_state_unif[1] >
         price_dist_steady_cumsum[price_id_consumer])
  }
  
  # state_idに変換し、各消費者の初期のstate_idを決める
  df$state_id[1] <-  state_df %>% 
    # mileageは0とする
    dplyr::filter(mileage_id == 1) %>% 
    dplyr::filter(price_id == price_id_consumer) %>% 
    dplyr::select(state_id) %>% 
    as.numeric()
  
  # Step 2: 各消費者について、状態変数、コントロール変数を逐次的に生成
  for (t in 1:(num_period - 1)) {
    # t期のstateを取得
    state_id_today <- df$state_id[t]
    
    # 価値関数に基づいて、購入するかどうかを決める
    if (V_CS[, 'V_not_buy'][state_id_today] + df$eps_type1_not_buy[t] > 
        V_CS[, 'V_buy'][state_id_today] + df$eps_type1_buy[t]){
      
      # 購入しない
      df$action[t] <- 0
      
      # 直面する遷移行列を定義
      trans_mat_cum_today <- trans_mat_cum$not_buy
      
    } else {
      # 購入する
      df$action[t] <- 1
      
      # 直面する遷移行列を定義
      trans_mat_cum_today <- trans_mat_cum$buy
    }
    
    # t+1期のstateを決める
    state_id_tomorrow <- 0
    exceed_trans_prob <- TRUE
    # 一様分布から生成した乱数が、遷移確率の累積分布の値を
    # 初めて下回ったところをt+1期の状態変数とする
    while (exceed_trans_prob) {
      state_id_tomorrow <- state_id_tomorrow + 1
      trans_prob <- trans_mat_cum_today[state_id_today, state_id_tomorrow]
      exceed_trans_prob <- (df$eps_unif[t] > trans_prob)
    }
    df$state_id[t + 1] <- state_id_tomorrow
  }
  return(df)
}


# 行列の(i,j)要素を出力する関数
mat_ij <- Vectorize(
  function(i, j, mat) {mat[i, j]},
  vectorize.args = c('i', 'j')
  )


# 静学的なロジットにおける対数尤度関数を定義する関数
logLH_stat <- function(theta, 
                       state_df, 
                       df) {
  
  # 選択毎の効用関数を求める
  U <- flow_utility(theta, state_df)
  # 選択確率を計算
  prob_C_stat <- exp(U) / rowSums(exp(U))
  # 対数尤度を計算
  sum(log(mat_ij(df$state_id, df$action + 1, prob_C_stat)))
}


# 不動点アルゴリズムにおける対数尤度関数を定義する関数
logLH_nfxp <- function(theta, 
                       beta, 
                       trans_mat, 
                       state_df, 
                       df){
  # パラメタより今期の効用を計算
  U <- flow_utility(theta, state_df)
  # 縮小写像アルゴリズムで得られた事前の価値関数をVとする
  V <- contraction(theta, beta, trans_mat, state_df)
  # 選択肢ごとの価値関数を計算
  V_CS <- U + beta * cbind(trans_mat$not_buy %*% V,
                           trans_mat$buy %*% V)
  # 選択確率を計算
  prob_C <- exp(V_CS) / rowSums(exp(V_CS))
  # 対数尤度を計算
  sum(log(mat_ij(df$state_id, df$action + 1, prob_C)))
}


# 行列インバージョンで期待価値関数を導出し、CCPを出力する関数を作成
policy_operator_mat_inv <- function(theta,
                                    CCP,
                                    beta,
                                    G, 
                                    state_df) {
  # 今期の効用を計算（本誌のu(i)に対応）
  U <- flow_utility(theta, state_df)
  # psiを計算（本誌のpsi(i)に対応）
  num_states <- nrow(state_df)
  psi <- Euler_const * matrix(1, num_states, num_choice) - log(CCP)
  # 行列計算により期待価値関数ベクトルを計算
  V <- 
    solve(diag(num_states) - 
            beta * (CCP[, 1] %*% matrix(1, 1, num_states) * G[,, 1] + 
                      CCP[, 2] %*% matrix(1, 1, num_states) * G[,, 2])) %*% 
    rowSums(CCP * (U + psi))
  # 選択肢ごとの価値関数を計算
  CV <- U + beta * cbind(G[,, 1] %*% V, G[,, 2] %*% V)
  # CCPを更新
  CCP <- exp(CV) / rowSums(exp(CV))
  
  return(CCP)
}

# 二段階推定のStep 2における尤度関数を定義する関数
likelihood_fun <- function(theta, 
                           CCP, 
                           df, 
                           beta, 
                           G, 
                           state_df, 
                           policy_operator) {
  # CCP operaterにより、CCPを更新する
  # 今回は行列インバージョンと有限依存性(finite dependence)アプローチのどちらか
  CCP <- policy_operator(theta, CCP, beta, G, state_df)
  # 疑似尤度を計算
  obj <- sum(log(mat_ij(df$state_id, df$action + 1, CCP)))
  return(obj)
}


# 有限依存性 (finite dependence) アプローチで期待価値関数を導出し、CCPを出力する関数
policy_operator_finite_dep <- function(theta, CCP, beta, G, state_df) {
  # 今期の効用を計算（本誌のu(i)に対応）
  U <- flow_utility(theta, state_df)
  # Finite dependenceに基づき、選択ごとの価値関数の差を計算
  CV_dif <- 
    U[, 2] - U[, 1] + 
    beta * (G[,, 2] %*% (-log(CCP[, 2])) - G[,, 1] %*% (-log(CCP[, 2])))
  # CCPを更新
  prob_buy <- exp(CV_dif) / (1 + exp(CV_dif))
  CCP <- cbind(1 - prob_buy, prob_buy)
  return(CCP)
}

# 入れ子不動点アルゴリズムにおける縮小写像を定義する関数
contraction_nfxp <- function(theta, 
                             beta, 
                             G, 
                             state_df) {
    # 価値関数の初期値
    num_states <- nrow(state_df)
    V_old <- matrix(0, nrow = num_states, ncol = 1)
    # パラメタより今期の効用を計算
    U <- flow_utility(theta, state_df)
    
    # 価値関数の差の初期値
    diff <- 1000
    # 縮小写像の誤差範囲
    tol_level <- 1.0e-12
    
    while (diff > tol_level) {
      # 価値関数を計算
      V_new <- log(rowSums(exp(U + beta * cbind(G[,, 1] %*% V_old, G[,, 2] %*% V_old)))) + Euler_const
      # 価値関数の更新による差を評価
      diff <- max(abs(V_new-V_old))
      # 価値関数を次のループに渡す(価値関数の更新)
      V_old <- V_new
    }
    # 得られた事前の価値関数を出力
    return(V_old)
}


# 入れ子不動点アルゴリズムにおいて期待価値関数を導出し、CCPを出力する関数
policy_operator_nfxp <- function(theta, 
                                 beta, 
                                 G, 
                                 state_df) {
    # パラメタより今期の効用を計算
    U <- flow_utility(theta, state_df)
    # 縮小写像アルゴリズムで得られた事前の価値関数をVとする
    V <- contraction_nfxp(theta, beta, G, state_df)
    # 選択肢ごとの価値関数を計算
    CV <- U + beta * cbind(G[,, 1] %*% V, G[,, 2] %*% V)
    # CCPを計算
    CCP <- exp(CV) / rowSums(exp(CV))
    return(CCP)
}

# 入れ子不動点アルゴリズムにおける尤度関数を定義する関数
likelihood_fun_nfxp <- function(theta,
                                df, 
                                beta, 
                                G, 
                                state_df){
  # CCP operaterにより、
  CCP <- policy_operator_nfxp(theta, beta, G, state_df)
  # 疑似尤度を計算
  obj <- sum(log(mat_ij(df$state_id, df$action + 1, CCP)))
  return(obj)
}
