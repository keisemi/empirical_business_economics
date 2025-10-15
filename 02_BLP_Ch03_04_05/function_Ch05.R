
# 市場シェアを計算する関数
calculate_mktshare <- function(
    datalist, 
    parameter, 
    delta, 
    theta2
) {
  # datalist, parameterの中身を読み込む
  X2 <- datalist$X2
  N <- length(datalist$marketindex)
  T <- length(unique(datalist$marketindex))
  draw_vec <- datalist$draw_vec
  marketindex <- datalist$marketindex
  mkt_denom_d <- datalist$mkt_denom_d
  Nsim <- parameter$Nsim
  K2 <- length(theta2)
  
  # 非線形の要素
  mu <- X2 %*% diag(x = theta2, nrow = K2) %*% matrix(draw_vec, nrow = K2)
  delta_mu <- delta %*% matrix(1, nrow = 1, ncol = Nsim) + mu
  exp_delta_mu <- exp(delta_mu)
  
  # ロジット確率の分母を計算する。これは、市場ごとに和を計算する。
  denom_temp <- t(t(exp_delta_mu) %*% mkt_denom_d)
  denom_outside <- exp(matrix(0, T, Nsim))
  denom_temp <- denom_temp + denom_outside
  
  # denom は (N × ns) 行列で、ロジット確率の分母を表している
  denom <- mkt_denom_d %*% denom_temp
  
  # 個々の選択確率を表す (N × ns) 行列
  s_jt_i <- exp_delta_mu / denom
  
  # シェアを表す (N × 1) 行列
  s_jt <- apply(s_jt_i, 1, mean)
  
  return(s_jt)
}

# Berryインバージョンによる平均効用の計算を行う関数
calculate_avg_utility_by_Berry_inversion <- function(
    datalist, 
    parameter, 
    delta_ini, 
    theta2
) {
  # 初期値の代入
  ShareVec <- datalist$ShareVec
  exp_delta_old <- exp(delta_ini)
  
  # 収束の設定
  tol <- 1e-11
  norm <- 1e+10
  iter <- 0
  while (norm > tol && iter < 1000){
    # マーケットシェアの計算
    pred_mkt_share <- calculate_mktshare(
      datalist, 
      parameter, 
      log(exp_delta_old), 
      theta2
    )
    
    # deltaを更新
    exp_delta <- exp_delta_old * ShareVec / pred_mkt_share
    # 距離を評価
    norm <- max(abs(exp_delta - exp_delta_old))
    # 計算結果を次のループに渡す
    exp_delta_old <- exp_delta 
    iter <- iter + 1
  }
  return(log(exp_delta))
}


# GMM目的関数を定義する関数
GMM_obj <- function(
    theta2, 
    parameter, 
    datalist, 
    option # option=0:最適化の際に指定, 1:推定値やmean utility(delta)を得る際に指定
) {
  # Contractionの初期値として、グローバル変数のdelta_globalを用いる。
  delta_ini <- delta_global
  
  # 縮小写像
  delta <- calculate_avg_utility_by_Berry_inversion(
    datalist, 
    parameter, 
    delta_ini, 
    theta2
  )
  
  # 平均効用を永続代入演算子でアップデートしている。 
  delta_global <<- delta
  
  # 目的関数に必要な行列を準備
  X1 <- datalist$X1
  Z <- datalist$Z
  
  # 重み行列 W
  if (datalist$weight_mat_option == "2SLS") {
    W <- solve(t(Z) %*% Z)
  } else if (datalist$weight_mat_option == "Ident") {
    W <- eye(dim(Z)[2])
  } else if (is.null(datalist$weight_mat_option)) {
    stop("SET 'datalist$weight_mat_option' !!")
  }
  
  # regressionにより線形パラメータを得る
  theta1 <- solve(t(X1) %*% Z %*% W %*% t(Z) %*% X1) %*% t(X1) %*% Z %*% W %*% t(Z) %*% delta
  # xiを計算
  Xi <- delta - X1 %*% theta1
  
  # 目的関数
  output <- t(Xi) %*% Z %*% W %*% t(Z) %*% Xi
  
  # optionによって戻り値を変える
  if (option == 0) {
    return(output = output)
  } else if (option == 1) {
    return(list(output = output, theta1 = theta1, delta = delta, Xi = Xi))
  } else {
    stop(
      "SET the argument 'option' to 0 or 1 !! \n 
      0: Specify when optimizing \n
      1: Specify when getting coefficient and mean utility (delta)"
    )
  }
}

# 標準誤差を計算する関数
calculate_standard_error <- function(
    theta2, 
    parameter, 
    datalist
) {
  param <- parameter
  delta_ini <- datalist$logitshare
  
  # 縮小写像
  delta <- calculate_avg_utility_by_Berry_inversion(
    datalist, 
    param, 
    delta_ini, 
    theta2
  )
  
  # 目的関数
  X1 <- datalist$X1
  
  # IV行列
  Z <- datalist$Z
  
  # 重み行列 W
  if (datalist$weight_mat_option == "2SLS") {
    W <- solve(t(Z) %*% Z)
  } else if (datalist$weight_mat_option == "Ident") {
    W <- eye(dim(Z)[2])
  } else if (is.null(datalist$weight_mat_option)) {
    STOP("SET 'datalist$weight_mat_option' !!")
  }
  
  # regressionにより線形パラメータを得る
  theta1 <- solve(t(X1) %*% Z %*% W %*% t(Z) %*% X1) %*% t(X1) %*% Z %*% W %*% t(Z) %*% delta
  
  # Residualを取得する。
  Xi <- delta - X1 %*% theta1
  
  # Omegaを計算する
  Omega_hat <- matrix(0, nrow = ncol(Z), ncol = ncol(Z))
  for (ii in 1:param$N) {
    Omega_hat <- Omega_hat + (matrix(Z[ii, ]) %*% t(matrix(Z[ii, ])) * Xi[ii]^2) / parameter$N
  }
  
  # Gradient of delta. 単純化のために数値微分で計算する。
  Ddelta <- matrix(0, nrow = parameter$N, ncol = length(theta2))
  
  for (k in 1:length(theta2)) {
    tempparam <- param
    temptheta2 <- theta2
    temptheta2[k] <- theta2[k] + 1e-6
    
    delta_add <- calculate_avg_utility_by_Berry_inversion(
      datalist, 
      tempparam, 
      delta_ini, 
      temptheta2
    )
    
    Ddelta[, k] <- (delta_add - delta) / 1e-6
  }
  
  G <- (parameter$N^(-1)) * t(Z) %*% cbind(-X1, Ddelta)
  
  # 漸近分散共分散行列
  AsyVarMat <- solve(t(G) %*% W %*% G) %*% t(G) %*% W %*% Omega_hat %*% W %*% G %*% solve(t(G) %*% W %*% G)
  # 漸近標準誤差
  Ase <- sqrt(diag(AsyVarMat) / parameter$N)
  
  return(Ase)
}


# 価格弾力性を計算する関数
calculate_elasticity <- function(
    datalist, 
    parameter, 
    theta1, 
    theta2, 
    delta
) {
  # データをunpackする
  X2 <- datalist$X2
  Nsim <- parameter$Nsim
  N <- parameter$N
  price <- datalist$X1[, colnames(datalist$X1) == "price"] %>% as.matrix()
  K2 <- length(theta2)
  draw_vec <- datalist$draw_vec
  uniquemarketindex <- datalist$uniquemarketindex
  
  # 以下はcalculate_mktshareとほぼ同じ
  mu <- X2 %*% diag(x = theta2, nrow = K2) %*% matrix(draw_vec, nrow = K2)
  
  denom_outside <- exp(matrix(0, T, Nsim))
  delta_mu <- delta %*% matrix(1, nrow = 1, ncol = Nsim) + mu
  exp_delta_mu <- exp(delta_mu)
  mkt_denom_d <- datalist$mkt_denom_d
  denom_temp <- t(t(exp_delta_mu) %*% mkt_denom_d)
  denom_temp <- denom_temp + denom_outside
  denom <- mkt_denom_d %*% denom_temp
  
  s_jt_i <- exp_delta_mu / denom
  
  # 価格パラメータ alpha_i
  draw_for_price <- matrix(draw_vec, nrow = K2)[1, ]
  
  # theta1のなかで2つ目、theta2のなかで1つ目に価格係数がある
  alpha_i <- theta1[2] + theta2[1] * draw_for_price
  
  # 弾力性行列を保存する空のリストを作成する
  elaslist <- as.list(rep(NA, T))
  
  # 各マーケットごとに弾力性行列を作成する．
  for (t in 1:T) {
    year_beg <- 2016 - T 
    # 製品の数
    J_t <- sum(marketindex == year_beg + t)
    
    # それぞれのマーケットに注目した(J*ns)の選択確率の行列を作成する
    # あるマーケットに対して，製品ごとに集計した選択確率の行列(J*ns)を作成する
    ag_model_s_i <- s_jt_i[marketindex == year_beg + t,]
    
    # 製品ごとの選択確率のベクトル(J * 1)
    ag_model_s <- apply(ag_model_s_i, 1, mean) %>% as.matrix()
    
    # 空の弾力性行列(J * J)を定義
    elasmat <- matrix(0, nrow = J_t, ncol = J_t)
    
    # あるマーケットに注目した，製品ごとの価格のベクトル(J*1)を作成する
    price_t <- price[marketindex == year_beg + t]
    
    # 弾力性を計算し，弾力性行列を作成する
    for (k in 1:J_t) {
      for (j in 1:J_t) {
        if (k != j) {
          # 交差価格弾力性
          elasmat[k, j] <- (-1)*price_t[k] * ag_model_s[j, ]^(-1) * 
            mean(alpha_i * ag_model_s_i[j, ] * ag_model_s_i[k, ])
        } else if (k == j) {
          # 自己価格弾力性
          elasmat[k, j] <- price_t[j] * ag_model_s[j, ]^(-1) * 
            mean(alpha_i * ag_model_s_i[j, ] * (1 - ag_model_s_i[j, ]))
        } 
      }
    }
    # 弾力性行列をマーケットごとに名前を付けて保存する．
    elaslist[[t]] <- elasmat
    names(elaslist)[[t]] <- sprintf("elasmat_%d", t + year_beg)
  }
  return(elaslist)
}

# 限界費用を推定する関数
estimate_merginal_cost <- function(
    data, 
    marketindex, 
    elasticity
) {
  # マーケットの数を取り出す
  T <- length(unique(data$year))
  
  # 限界費用を保存する空のリストを作成する
  MC_list <- as.list(rep(NA, T))
  
  year_beg <- min(data$year) - 1 
  
  # 各マーケットごとに限界費用計算する
  for (t in 1:T) {
    # あるマーケットの弾力性行列を取り出す
    elasmat_t <- elasticity[[t]] 
    
    # あるマーケットの大きさ(車種の数)
    J_t <- sum(marketindex == year_beg + t)
    
    # あるマーケットのデータを抽出
    data_t <- data[data$year == year_beg + t, ]
    
    # あるマーケットの価格とマーケットシェアを取り出す
    Pricevec_t <- data_t$price
    Sharevec_t <- data_t$share
    
    # 空の所有構造行列と弾力性の微分部分を用意する
    Ownership_t <- matrix(0, J_t, J_t)
    # 各弾力性に対して計算をする
    for (j in 1:J_t) {
      for (k in 1:J_t) {
        
        # 車種jと車種kが同じ企業で生産されていれば1を入れる
        if (data_t$Maker[j] == data_t$Maker[k]) {
          Ownership_t[j, k] <- 1
        }
      }
    }
    Derivative_t <- - (elasmat_t) * (kronecker(matrix(1, J_t, 1), t(Sharevec_t))) / 
      (kronecker(matrix(1, 1, J_t), Pricevec_t))
    # 所有構造行列に弾力性の微分部分を掛けることでDeltaを作成する
    Delta_t <- Derivative_t * Ownership_t
    
    # 限界費用を計算する
    Marginal_Cost_t <- Pricevec_t - (solve(Delta_t) %*% Sharevec_t)
    
    # データフレーム化してNameと企業名を追加する(NameIDはleft_joinに必要)
    Marginal_Cost_t_D <- data.frame(Name = data_t$Name, NameID = data_t$NameID, 
                                    Maker = data_t$Maker, 
                                    price = Pricevec_t, mc = Marginal_Cost_t) %>%
      mutate(margin = (price - mc) / price)
    # 限界費用をマーケットごとに名前を付けて保存する．
    MC_list[[t]] <- Marginal_Cost_t_D
    names(MC_list)[[t]] <- sprintf("Marginal_Cost_%d", t + year_beg)
  }
  return(MC_list)
}

# 価格を更新する関数
update_price <- function(
    datalist, 
    p_old, 
    Ownership, 
    parameter, 
    theta1, 
    theta2, 
    mc, 
    Xi
) {
  # datalist内の価格を更新
  datalist$X1[, 'price'] <- as.matrix(p_old)
  datalist$X2 <- as.matrix(p_old)
  
  X1 <- datalist$X1
  
  # 市場シェアを求める
  delta <- (X1 %*% theta1) + Xi
  Sharevec <- calculate_mktshare(datalist, parameter, delta, theta2)
  
  # 価格弾力性を求める
  elas <- calculate_elasticity(datalist, parameter, theta1, theta2, delta)
  elasmat <- elas$elasmat_2016
  
  J <- length(p_old)
  Derivative <- - (elasmat) * (kronecker(matrix(1, J, 1), t(Sharevec))) / (kronecker(matrix(1, 1, J), p_old))
  Delta <- Ownership * Derivative
  p_new <- mc + (solve(Delta) %*% Sharevec)
  
  return(p_new)
}


# iterationで均衡価格を求める関数
solve_equilibrium_price <- function(
    datalist, 
    p_ini, 
    Ownership, 
    parameter, 
    theta1, 
    theta2, 
    mc, 
    Xi
) {
  # 閾値を設定する
  lambda <- 1e-6
  p_old <- p_ini
  distance <- 10000
  while (distance > lambda) {
    p_new <- update_price(
      datalist, 
      p_old, 
      Ownership, 
      parameter, 
      theta1, 
      theta2, 
      mc, 
      Xi
      )
    distance <- max(abs(p_new - p_old))
    p_old <- p_new
  }
  return(p_new)
}

# シミュレートされた価格でのシェアを与える関数
calculate_mktshare_sim <- function(
    datalist, 
    p, 
    parameter, 
    theta1, 
    theta2, 
    Xi
) {
  # datalist内の価格を更新
  datalist$X1[, 'price'] = p
  datalist$X2 = p
  X1 = datalist$X1
  
  # 市場シェアを求める
  delta = (X1 %*% theta1) + Xi
  Sharevec = calculate_mktshare(datalist, parameter, delta, theta2)
  
  return(Sharevec)
}

# 消費者余剰を計算する関数
calculate_consumer_surplus = function(
    datalist, 
    p, 
    parameter, 
    theta1, 
    theta2, 
    Xi, 
    HH
) {
  # datalist内の価格を更新
  datalist$X1[, 'price'] <- p
  datalist$X2 <- as.matrix(p)
  X1 <- datalist$X1
  
  # deltaを求める
  delta <- (X1 %*% theta1) + Xi
  
  # calculate_mktshareの分子と同じ計算
  X2 <- datalist$X2
  N <- length(datalist$marketindex)
  draw_vec <- datalist$draw_vec
  marketindex <- datalist$marketindex
  mkt_denom_d <- datalist$mkt_denom_d
  Nsim <- parameter$Nsim
  K2 <- length(theta2)
  
  # 非線形の要素
  mu <- X2 %*% diag(x = theta2, nrow = K2) %*% matrix(draw_vec, nrow = K2)
  
  # prepare exp_V
  V <- delta %*% matrix(1, nrow = 1, ncol = Nsim) + mu
  exp_V <- exp(V)
  # 分子
  # 第2項はアウトサイド・グッズに対する
  numerator <- log(colSums(exp_V) + exp(matrix(0, 1, Nsim)))
  
  # 価格パラメータalpha_i
  draw_for_price <- matrix(draw_vec, nrow = K2)[1, ]
  # theta1のなかで2つ目、theta2のなかで1つ目に価格係数がある
  alpha_i <- - (theta1[2] + theta2[1] * draw_for_price)
  
  # 消費者余剰
  CS <- mean(numerator / alpha_i) * HH
  return(CS)
}

# 利潤を計算する関数
calculate_profit <- function(
    Maker, 
    price,
    mc, 
    share, 
    HH
) {
  
  pro_rev_dt <- data.frame(
    Maker, 
    profit_each = (price - mc) * share * HH,
    revenue_each = price * share * HH
  )
  
  pro_rev <- pro_rev_dt %>%
    dplyr::group_by(Maker) %>%
    dplyr::summarise(profit = sum(profit_each), revenue = sum(revenue_each)) %>%
    dplyr::ungroup() %>%
    as.data.frame()
  return(pro_rev)
}

# 利潤と収入の表を作成する関数
generate_profit_revenue_table <- function(
    pro_rev_2016, 
    pro_rev_NippyoA, 
    pro_rev_NippyoB
) {
  
  result_df <- dplyr::tibble(
    Maker = pro_rev_2016$Maker,
    col2 = pro_rev_NippyoA$profit - pro_rev_2016$profit,
    col3 = pro_rev_NippyoA$revenue - pro_rev_2016$revenue,
    col4 = pro_rev_NippyoB$profit - pro_rev_2016$profit,
    col5 = pro_rev_NippyoB$revenue - pro_rev_2016$revenue
  ) 
  
  Total <- dplyr::tibble(
    Maker = "Total",
    col2 = sum(result_df$col2),
    col3 = sum(result_df$col3),
    col4 = sum(result_df$col4),
    col5 = sum(result_df$col5)
  )
  result_df <- rbind(result_df, Total)
  
  colnames(result_df) <- c("Maker", "Profits", "Revenues", "Profits", "Revenues")
  
  return(result_df)
}


# 特定の企業のみ限界費用を定数倍した時の、総余剰の変化を計算する関数
calculate_surplus_change <- function(
    cost_red, 
    cost_red_firm, 
    Ownership, 
    data, 
    datalist, 
    parameter, 
    theta1, 
    theta2, 
    HH, 
    p_pre, 
    pro_rev_pre, 
    CS_pre
) {
  # 特定の企業のみ、限界費用を一律に係数倍
  data <- data %>%
    dplyr::mutate(mc = if_else(Maker %in% cost_red_firm, mc * cost_red, mc))
  mc <- data$mc
  
  # 均衡価格
  p_post <- solve_equilibrium_price(
    datalist, 
    p_pre, 
    Ownership, 
    parameter, 
    theta1, 
    theta2, 
    mc, 
    Xi
  )
  
  # CVの計算
  CV <- (calculate_consumer_surplus(
    datalist, 
    p_post, 
    parameter, 
    theta1, 
    theta2, 
    Xi, 
    HH
    ) - CS_pre)
  
  # 利益と売上の計算
  share_post <- calculate_mktshare_sim(
    datalist, 
    p_post, 
    parameter, 
    theta1, 
    theta2, 
    Xi
  )
  
  pro_rev_post <- calculate_profit(
    data$Maker, 
    p_post, 
    data$mc, 
    share_post, 
    HH
  ) 
  
  # 総余剰の変化
  TS_Chage <- CV + sum(pro_rev_post$profit - pro_rev_pre$profit)
  obj <- TS_Chage
  
  return(obj)
}
