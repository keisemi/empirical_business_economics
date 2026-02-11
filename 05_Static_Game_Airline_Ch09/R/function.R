
# 予測スコアと真のラベルからROC曲線のAUC（Area Under the Curve）を計算する関数
calculate_AUC <- function(s, Y){
  
  pred <- ROCR::prediction(s, Y)
  auc <- performance(pred, "auc")@y.values[[1]]
  
  return(auc)
}


# B-splineを使った回帰式を生成する関数
generate_bspline_formula = function(var, comp, degree){
  
  formula = paste0("y_", comp, " ~Train*(")
  
  for(i in var){
    formula = paste0(formula, 
                      paste0("bs(", i, ",knots = c(quantile(", i, ",.5)), degree = ", degree, ")+"))
  }
  
  formula = str_sub(formula, end = -2)
  formula = paste0(formula, ")")
  
  formula = as.formula(formula)
  
  return(formula)
}


# モーメント条件を計算する関数
calculate_moment = function(par, X, Z, y){
  
  # 予測確率の計算
  p =  1 / (1 + exp(- X %*% par))
  
  # モーメント条件の計算
  g = t(Z) %*% (y - p)
  
  return(g)
}


# 係数の標準誤差を計算する関数
calculate_standard_error = function(par, X, Z, y){
  
  # 予測確率の計算
  p = 1 / (1 + exp(- X %*% par))
  
  N = nrow(X)
  g = kronecker(t(rep(1,ncol(X))),(p-y)) * Z
  Omega = t(g) %*% g/N
  G = t(Z) %*% (kronecker(t(rep(1,ncol(X))),p*(1-p)) * X )/N
  
  #just identifiedを仮定した分散
  V = solve(G) %*% Omega %*% t(solve(G)) / N
  
  return(sqrt(diag(V)))
}


# 予測確率を計算する関数
calculate_predict_probability = function(par, f_spec = NA){
  
  # "p_"で始まる要素を抽出
  which_p = grepl("^p_",　names(par))
  
  # f_specが指定されていれば、それを p_name（最後に選択する列名）として使用
  # f_spec = NAであれば、par の中に含まれる "p_" で始まる要素を自動で採用
  if(!is.na(f_spec)){
    p_name = f_spec
  } else{
    p_name = names(par)[which_p]
  }
  
  par = c(par["(Intercept)"], 
          par["Distance"], 
          par["Population"], 
          par["Pop_Square"], 
          par["Train"], 
          par["Flight"],
          par["Train:Distance"], 
          par["Train:Population"], 
          par["Train:Flight"], 
          par["Train:Pop_Square"], 
          par[which_p])
  
  
  X = df_long %>% 
    select(Constant, Distance, Population, Pop_Square, Train, Flight, 
           Distance_Train, Population_Train, Flight_Train, Pop_Square_Train, any_of(p_name)) %>%
    as.matrix()
  
  # 確率の計算
  p = 1 / (1 + exp(- X %*% par))
  p = p %>% as.vector
  
  return(p)
}


# 対数尤度を計算する関数
caluculate_log_likelihood = function(par, f_spec = NA){
  
  p = calculate_predict_probability(par, f_spec)
  y = df_long$y
  
  likelihood = y*log(p) + (1 - y)*log(1 - p)
  
  return(sum(likelihood))
}


# 最適反応を計算する関数
calculate_best_response = function(p_opp, par, X){
  
  X = cbind(X, p_opp)
  p = 1 / (1 + exp(- X %*% par))
  p = p %>% as.vector
  
  return(p)
}


# 相手プレイヤーの最適反応を計算する関数
calculate_best_response_opponent = function(p_opp, par, X){
  
  n_market = nrow(X) / 2
  p = calculate_best_response(p_opp, par, X)
  p_opp = c(p[(n_market + 1):(2*n_market)], p[1:n_market])
  
  return(p_opp)
}


# 最適反応を繰り返し計算することで均衡の確率ベクトルを求める関数
find_equilibrium = function(p_init, par, X, tol = 1e-6, max_iter = 1000){
  
  n_market = nrow(X) / 2
  p_opp = c(p_init[(n_market + 1):(2*n_market)], p_init[1:n_market])
  
  for (i in 1:max_iter){
    p_opp_new = calculate_best_response_opponent(p_opp,par,X)
    # 収束判定
    if (sum((p_opp_new - p_opp)^2) < tol){
      break
    }
    p_opp = p_opp_new
  }
  
  # 均衡確率ベクトル
  p_eq = c(p_opp[(n_market + 1):(2*n_market)], p_opp[1:n_market])
  
  return(p_eq)
}


# ANAとJAL2社の静学参入ゲームの均衡確率ベクトルを計算する関数
find_equilibrium_ANA_JAL = function(ids, optimal, par, df_long = df_long){
  
  X = df_long %>% 
    filter(id %in% ids) %>% 
    select(Constant, Distance, Population, Pop_Square, Train, Flight,
           Distance_Train, Population_Train, Flight_Train, Pop_Square_Train) %>%
    as.matrix()
  
  par = c(par["(Intercept)"],
          par["Distance"],
          par["Population"],
          par["Pop_Square"],
          par["Train"],
          par["Flight"],
          par["Train:Distance"],
          par["Train:Population"],
          par["Train:Flight"],
          par["Train:Pop_Square"],
          par[grepl("^p_",names(par))])
  
  n_market = nrow(X) / 2
  
  if (optimal == "ANA"){
    p_init = c(rep(1, n_market), rep(0, n_market)) # ANA優位の初期値を設定 (ANA:p = 1, JAL:p = 0) 
  }else{
    p_init = c(rep(0, n_market), rep(1, n_market)) # JAL優位の初期値を設定 (ANA:p = 0, JAL:p = 1) 
  }
  
  p_eq = find_equilibrium(p_init, par, X)
  
  return(p_eq)
}
