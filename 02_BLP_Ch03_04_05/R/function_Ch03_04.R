
# 価格をインプットとして、販売台数を返す関数
calculate_sales <- function(
    price_cand, 
    year, 
    modelID_target, 
    dt, 
    est_param
) {
  
  dt_result <- dt %>% 
    dplyr::filter(year == 2016) %>% 
    dplyr::mutate(temp_price = ifelse(NameID == modelID_target, 
                                      price_cand, 
                                      price)) %>% 
    dplyr::mutate(
      delta = estparam[1] + estparam[2] * temp_price + estparam[3] * hppw + 
        estparam[4] * FuelEfficiency + estparam[5] * size + xi_fit
    ) %>% 
    dplyr::mutate(denom = 1 + sum(exp(delta))) %>%
    dplyr::mutate(pred_sales = (exp(delta) / denom) * HH) %>%
    dplyr::filter(NameID == modelID_target)
  
  return(dt_result$pred_sales)
}

# 収入の最大化を実行する関数
maximize_revenue <- function(
    price_cand, 
    year, 
    modelID_target, 
    dt, 
    est_param
) {
  
  quantity <- calculate_sales(
    price_cand,
    year, 
    modelID_target, 
    dt, 
    est_param
  )
  
  revenue <- price_cand * quantity
  
  return(revenue)
}


# 市場シェアを計算する関数
calculate_mktshare <- function(
    datalist, 
    parameter, 
    delta
) {
  # データをunpackする
  X2 <- datalist$X2
  ns <- parameter$Nsim
  T <- parameter$T
  N <- parameter$N
  
  draw_vec <- datalist$draw_vec
  
  theta2 <- parameter$theta2
  K2 <- length(theta2)
  
  draw_vec <- datalist$draw_vec[1:(ns * K2)]
  
  marketindex <- datalist$marketindex
  
  # 非線形の要素
  mu <- X2 %*% diag(x = theta2, nrow = length(theta2)) %*% matrix(draw_vec, nrow = K2)
  
  denom_outside <- exp(matrix(0, T, ns))
  
  delta_mu <-  delta %*% matrix(1, nrow = 1, ncol = ns) + mu
  exp_delta_mu <- exp(delta_mu)
  
  # ロジット確率の分母として、市場ごとの和を計算
  tempmat <- datalist$tempmat
  
  denom_temp <- t(t(exp_delta_mu) %*% tempmat)
  denom_temp <- denom_temp + denom_outside
  
  # denom は (N × ns) 行列で、ロジット確率の分母を表している
  denom <- tempmat %*% denom_temp
  
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
    delta_ini
) {
  tol <- 1e-11
  share_obs <- datalist$ShareVec
  norm <- 1e+10
  
  delta_old <- delta_ini
  exp_delta_old <- exp(delta_old)
  
  iter <- 0
  
  while (norm > tol && iter < 1000) {
    
    pred_mkt_share <- calculate_mktshare(
      datalist, 
      parameter, 
      delta_old
    )
    
    exp_delta <- exp_delta_old * share_obs / pred_mkt_share
    
    t <- abs(exp_delta - exp_delta_old)
    norm <- max(t)
    
    # 更新
    exp_delta_old <- exp_delta 
    
    delta_old <- log(exp_delta)
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
  param <- parameter
  param$theta2 <- theta2
  
  # Contractionの初期値として、log(s_j)-log(s_0) を用いる。
  delta_ini <- datalist$delta_ini
  
  # 縮小写像
  delta <- calculate_avg_utility_by_Berry_inversion(
    datalist, 
    param, 
    delta_ini
  )
  
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
  beta_hat <- solve(t(X1) %*% Z %*% W %*% t(Z) %*% X1) %*% t(X1) %*% Z %*% W %*% t(Z) %*% delta
  
  Xi <- delta - X1 %*% beta_hat
  
  # 目的関数
  output <- t(Xi) %*% Z %*% W %*% t(Z) %*% Xi
  
  # optionによって戻り値を変える
  if (option == 0) {
    return(output = output)
  } else if (option == 1) {
    return(list(output = output, beta_hat = beta_hat, delta = delta))
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
  param$theta2 <- theta2
  
  delta_ini <- datalist$logitshare
  
  # 縮小写像
  delta <- calculate_avg_utility_by_Berry_inversion(
    datalist, 
    param, 
    delta_ini
  )
  
  # 目的関数に必要な行列を準備
  X1 <- datalist$X1
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
  beta_hat <- solve(t(X1) %*% Z %*% W %*% t(Z) %*% X1) %*% t(X1) %*% Z %*% W %*% t(Z) %*% delta
  
  Xi <- delta - X1 %*% beta_hat
  
  # Omega
  Omega_hat <- matrix(0, nrow = ncol(Z), ncol = ncol(Z))
  
  for (ii in 1:param$N) {
    Omega_hat <- Omega_hat + 
      (matrix(Z[ii, ]) %*% t(matrix(Z[ii, ])) * Xi[ii]^2) / parameter$N
  }
  
  # Gradient of delta. 単純化のために数値微分で計算
  Ddelta <- matrix(0, nrow = parameter$N, ncol = length(theta2))
  
  for (k in 1:length(theta2)) {
    tempparam <- param
    tempparam$theta2[k] <- theta2[k] + 1e-6
    
    delta_add <- calculate_avg_utility_by_Berry_inversion(
      datalist, 
      tempparam, 
      delta_ini
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
    beta_hat, 
    delta
) {
  
  X2 <- datalist$X2
  ns <- parameter$Nsim
  N <- parameter$N
  price <- datalist$X1[, colnames(datalist$X1) == "price"] %>% as.matrix()
  theta2 <- parameter$theta2
  K2 <- length(theta2)
  draw_vec <- datalist$draw_vec[1:(ns * K2)]
  uniquemarketindex <- datalist$uniquemarketindex
  
  # 以下はcalculate_mktshareとほとんど同じ
  mu <- X2 %*% diag(x = theta2, nrow = length(theta2)) %*% matrix(draw_vec, nrow = K2)
  
  denom_outside <- exp(matrix(0, T, ns))
  delta_mu <- delta %*% matrix(1, nrow = 1, ncol = ns) + mu
  exp_delta_mu <- exp(delta_mu)
  tempmat <- datalist$tempmat
  denom_temp <- t(t(exp_delta_mu) %*% tempmat)
  denom_temp <- denom_temp + denom_outside
  
  denom <- tempmat %*% denom_temp
  
  # 選択確率
  s_jt_i <- exp_delta_mu / denom
  
  # 価格パラメータ
  draw_for_price <- matrix(draw_vec, nrow = K2)[1, ]
  alpha_i <- beta_hat[rownames(beta_hat) == "価格：平均"] +
    beta_hat[rownames(beta_hat) == "価格：標準偏差"] * draw_for_price
  
  # 弾力性行列を保存する空のリストを作成する
  elaslist <- as.list(rep(NA, T))
  
  # 各マーケットごとに弾力性行列を作成する．
  for (t in 1:T) {
    
    # 製品の数
    J_t <- sum(marketindex == 2005 + t)
    
    # それぞれのマーケットに注目した(J*ns)の選択確率の行列を作成する
    
    # あるマーケットに対して，製品ごとに集計した選択確率の行列(J*ns)を作成する
    ag_model_s_i <- s_jt_i[marketindex == 2005 + t,]
    
    # 製品ごとの選択確率のベクトル(J * 1)
    ag_model_s <- apply(ag_model_s_i, 1, mean) %>% as.matrix()
    
    # 空の弾力性行列(J * J)を定義
    elasmat <- matrix(0, nrow = J_t, ncol = J_t)
    
    # あるマーケットに注目した，製品ごとの価格のベクトル(J*1)を作成する
    price_t <- price[marketindex == 2005 + t]
    
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
            mean(alpha_i * ag_model_s_i[j, ] * (1 - ag_model_s_i[j, ]) )
        } 
      }
    }
    
    # 弾力性行列をマーケットごとに名前を付けて保存する．
    elaslist[[t]] <- elasmat
    names(elaslist)[[t]] <- sprintf("elasmat_%d", t + 2005)
  }
  return(elaslist)
}



# 価格に対する収入を計算する関数
calculate_revenue <- function(
    price_cand, 
    data, 
    datalist, 
    result2, 
    parameter, 
    option
) {
  
  # ベータードの限界費用（今回はラーナーの公式から簡便的に定義する）
  mc_betado <- 3.198 * (1 - 1 / abs(-2.16720791))
  
  tempprice <- data$price
  tempprice[data$NameID == betard & data$year == 2016] <- price_cand
  
  datalist_temp <- datalist
  datalist_temp$X1[, "price"] <- tempprice
  datalist_temp$X2[, "price"] <- tempprice
  
  # 新たな価格を用いて平均効用を計算し直す
  org_xi <- result2$delta - datalist$X1 %*% result2$beta_hat
  new_delta <- datalist_temp$X1 %*% result2$beta_hat + org_xi
  
  mktshare <- calculate_mktshare(
    datalist = datalist_temp, 
    parameter = parameter, 
    delta = new_delta
  )
  
  quant <- mktshare * data$HH
  revenue <- tempprice * quant
  
  # ベータードのみの収入
  revenuevec <- revenue[data$NameID == betard & data$year == 2016]
  
  # 日評自動車全体の収入
  revenuevec2 <- sum(revenue[data$NameID %in% NIPPYOautoIDvec & data$year == 2016])
  
  # ベータードのみの利潤
  pivec <- revenuevec - mc_betado * quant[data$NameID == betard & data$year == 2016]
  
  # ベータードのみの利潤+日評自動車の他の車種の収入
  pivec2 <- revenuevec2 - mc_betado * quant[data$NameID == betard & data$year == 2016]
  
  if (option == "own"){
    return(revenuevec)
  } else if (option == "total") {
    return(revenuevec2)
  } else if (option == "ownpi") {
    return(pivec)
  } else if (option == "totalpi") {
    return(pivec2)
  }
}


