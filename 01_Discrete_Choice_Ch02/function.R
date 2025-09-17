# define functions for Ch02


# 消費者属性を入れない多項ロジットにおいて、選択確率を算出する関数
f_logit_prob <- function(alpha_Kinoko, 
                         alpha_Takenoko, 
                         beta, 
                         price_Kinoko, 
                         price_Takenoko){
  # きのこ、たけのこから得られるそれぞれの効用の算出
  util_Kinoko <- alpha_Kinoko - beta*price_Kinoko
  util_Takenoko <- alpha_Takenoko - beta*price_Takenoko
  
  # 選択確率の算出
  prob_Kinoko <- exp(util_Kinoko) / (1 + exp(util_Kinoko) + exp(util_Takenoko))
  prob_Takenoko <- exp(util_Takenoko) / (1 + exp(util_Kinoko) + exp(util_Takenoko))
  prob_Other <- 1 - (prob_Kinoko + prob_Takenoko)
  
  return(cbind(prob_Kinoko, prob_Takenoko, prob_Other))
}



# 各個人の選択確率を全て計算した後に、足し合わせて需要を得る関数
f_demand <- function(alpha_Kinoko_vec,
                     alpha_Takenoko_vec,
                     beta_vec,
                     price_Kinoko,
                     price_Takenoko){
  
  R = length(alpha_Kinoko_vec)  # 設定する潜在消費者数
  # 236人のそれぞれの消費者属性を代表として、観察されない異質性を1000回ずつドローして得られた潜在消費者236000人
  
  # 結果を保存するベクトルを事前に準備
  prob_Kinoko = numeric(R)
  prob_Takenoko = numeric(R)
  prob_Other = numeric(R)
  
  for (r in 1:R){
    result = f_logit_prob(alpha_Kinoko_vec[r], 
                          alpha_Takenoko_vec[r], 
                          beta_vec[r], 
                          price_Kinoko, 
                          price_Takenoko)
    prob_Kinoko[r] = result[1]
    prob_Takenoko[r] = result[2]
    prob_Other[r] = result[3]
  }
  
  # 選択確率をすべての消費者について足し合わせて、各オプションの需要を得る
  return(c(sum(prob_Kinoko), sum(prob_Takenoko), sum(prob_Other)))
}



# 消費者属性を入れない多項ロジットにおける対数尤度関数を定義する関数
f_likelihood_logit <- function(param,
                               data){
  
  alpha_Kinoko <- param[1]
  alpha_Takenoko <- param[2]
  beta   <- param[3]
  
  # 各選択フェーズでの選択確率を予測
  P1 <- f_logit_prob(alpha_Kinoko, alpha_Takenoko, beta, 200, 200)
  P2 <- f_logit_prob(alpha_Kinoko, alpha_Takenoko, beta, 180, 200)
  P3 <- f_logit_prob(alpha_Kinoko, alpha_Takenoko, beta, 200, 170)
  P4 <- f_logit_prob(alpha_Kinoko, alpha_Takenoko, beta, 220, 200)
  P5 <- f_logit_prob(alpha_Kinoko, alpha_Takenoko, beta, 190, 210)
  
  pred <- rbind(P1, P2, P3, P4, P5)
  pred <- as_tibble(pred)
  pred$occasion <- c("Q1", "Q2", "Q3", "Q4", "Q5")
  
  # 各選択フェーズに対する実際の選択と、予測された選択確率を対応させる
  result <- data %>%
    left_join(pred, by = "occasion") %>%
    mutate(choice_prob = case_when(choice == 0 ~ prob_Other,
                                    choice == 1 ~ prob_Kinoko,
                                    choice == 2 ~ prob_Takenoko)) %>%
    mutate(log_choice_prob = log(choice_prob)) %>%
    select(log_choice_prob)
  
  # 和を取って対数尤度を計算する
  likelihood <- sum(result)
  
  return(likelihood)
  
}



# 消費者属性を入れた多項ロジットにおいて、選択確率を算出する関数
f_logit_prob_with_atr <- function(alpha_Kinoko, 
                                  alpha_Takenoko, 
                                  beta, 
                                  price_Kinoko, 
                                  price_Takenoko, 
                                  param_atr_Kinoko, 
                                  param_atr_Takenoko, 
                                  param_atr_price,atr){
  # きのこ、たけのこから得られるそれぞれの効用の算出
  # param_atr_〇〇が各属性との交差項の係数のベクトル、atrが消費者属性のベクトル
  util_Kinoko <- alpha_Kinoko + sum(param_atr_Kinoko*atr) - price_Kinoko*(beta + sum(param_atr_price*atr))
  util_Takenoko <- alpha_Takenoko + sum(param_atr_Takenoko*atr) - price_Takenoko*(beta + sum(param_atr_price*atr))
  
  # 選択確率の算出
  prob_Kinoko <- exp(util_Kinoko) / (1 + exp(util_Kinoko) + exp(util_Takenoko))
  prob_Takenoko <- exp(util_Takenoko) / (1 + exp(util_Kinoko) + exp(util_Takenoko))
  prob_Other <- 1 - (prob_Kinoko + prob_Takenoko)
  
  return(cbind(prob_Kinoko, prob_Takenoko, prob_Other))

}



# 消費者属性を入れた多項ロジットにおける対数尤度関数を定義する関数
f_likelihood_logit_with_atr <- function(param,
                                        data){
  # 空のデータフレーム
  pred <- matrix(ncol = 3, nrow = 1180)
  colnames(pred) <- c("prob_Kinoko", "prob_Takenoko","prob_Other")
  ID <- numeric(236)
  atr <- numeric(5)
  beta <- param[1]
  alpha_Kinoko <- param[2]
  alpha_Takenoko <- param[3]
  param_atr_Kinoko <- param[c(6, 9, 12, 17, 18)]
  param_atr_Takenoko <- param[c(5, 8, 11, 15, 16)]
  param_atr_price <- param[c(4, 7, 10, 13, 14)]
  
  # 各回答者の選択確率を予測する
  for(i in 1:236){
    # 消費者属性を取得
    gender <- data$gender[(5*(i-1)+1)]
    familyhouse <- data$familyhouse[(5*(i-1)+1)]
    adult <- data$adult[(5*(i-1)+1)]
    kansai <- data$kansai[(5*(i-1)+1)]
    oversea <- data$oversea[(5*(i-1)+1)]
    atr[1:5] <- c(gender, familyhouse, adult, kansai, oversea)
    
    # 各選択フェーズでの選択確率を予測
    P1 <- f_logit_prob_with_atr(alpha_Kinoko,
                                alpha_Takenoko,
                                beta,
                                200,
                                200,
                                param_atr_Kinoko,
                                param_atr_Takenoko,
                                param_atr_price,atr)
    
    P2 <- f_logit_prob_with_atr(alpha_Kinoko, 
                                alpha_Takenoko, 
                                beta, 
                                180, 
                                200, 
                                param_atr_Kinoko,
                                param_atr_Takenoko, 
                                param_atr_price,atr)
    
    P3 <- f_logit_prob_with_atr(alpha_Kinoko, 
                                alpha_Takenoko, 
                                beta, 
                                200, 
                                170, 
                                param_atr_Kinoko,
                                param_atr_Takenoko, 
                                param_atr_price,atr)
    
    P4 <- f_logit_prob_with_atr(alpha_Kinoko, 
                                alpha_Takenoko, 
                                beta, 
                                220, 
                                200, 
                                param_atr_Kinoko,
                                param_atr_Takenoko, 
                                param_atr_price,atr)
    
    P5 <- f_logit_prob_with_atr(alpha_Kinoko, 
                                alpha_Takenoko, 
                                beta, 
                                190, 
                                210, 
                                param_atr_Kinoko,
                                param_atr_Takenoko, 
                                param_atr_price,atr)
    # 格納
    pred[(5*(i-1)+1),] <- P1
    pred[(5*(i-1)+2),] <- P2
    pred[(5*(i-1)+3),] <- P3
    pred[(5*(i-1)+4),] <- P4
    pred[(5*(i-1)+5),] <- P5
    ID[(5*(i-1)+1):(5*(i-1)+5)] <- rep(data$ID[(5*(i-1)+1)], 5)
  }
  pred <- as.data.frame(pred)
  
  # IDと設問(マージ用)
  pred$ID<- ID
  pred$occasion <- rep(c("Q1", "Q2", "Q3", "Q4", "Q5"), 236)
  
  # 各選択フェーズに対する実際の選択と、予測された選択確率を対応させる
  result <- data %>%
    left_join(pred, by = c("ID","occasion")) %>%
    mutate(choice_prob = case_when(choice == 0 ~ prob_Other,
                                   choice == 1 ~ prob_Kinoko,
                                   choice == 2 ~ prob_Takenoko)) %>%
    mutate(log_choice_prob = log(choice_prob)) %>%
    select(log_choice_prob)
  
  # 和を取って対数尤度を計算する
  likelihood <- sum(result)
  
  return(likelihood)
  
}



# 消費者属性を入れないランダム係数ロジットにおいて、選択確率を算出する関数
f_rcdclogit_prob <- function(alpha_Kinoko, 
                             alpha_Takenoko, 
                             beta, 
                             price_Kinoko, 
                             price_Takenoko, 
                             sigma_Kinoko, 
                             sigma_Takenoko, 
                             sigma_price, 
                             tau){
  # きのこ、たけのこから得られるそれぞれの効用の算出
  # sigma_〇〇は各パラメタの標準偏差、tauは標準正規分布に従う乱数を(ドロー回数)行×3列発生させたもの
  alpha_Kinoko_draw <- alpha_Kinoko + sigma_Kinoko*tau[, 1]
  alpha_Takenoko_draw <- alpha_Takenoko + sigma_Takenoko*tau[, 2]
  beta_draw <- beta + sigma_price*tau[, 3]
  util_Kinoko <- alpha_Kinoko_draw - price_Kinoko*beta_draw
  util_Takenoko <- alpha_Takenoko_draw - price_Takenoko*beta_draw
  
  # 選択確率の算出
  prob_Kinoko <- exp( util_Kinoko ) / ( 1 + exp( util_Kinoko ) + exp( util_Takenoko ))
  prob_Takenoko <- exp( util_Takenoko ) / ( 1 + exp( util_Kinoko ) + exp( util_Takenoko ))
  prob_Other <- 1 - (prob_Kinoko + prob_Takenoko)
  
  return( cbind(prob_Kinoko, prob_Takenoko, prob_Other))
}



# 消費者属性を入れないランダム係数ロジットにおける対数尤度関数を定義する関数
f_likelihood_rcdclogit <- function(param,
                                   data){
  # 空のデータフレーム
  pred <- matrix(ncol = 3, nrow = 1180*times_draw)
  colnames(pred) <- c("prob_Kinoko", "prob_Takenoko","prob_Other")
  ID <- numeric(236*times_draw)
  beta <- param[1]
  alpha_Kinoko <- param[2]
  alpha_Takenoko <- param[3]
  sigma_price <- param[4]
  sigma_Kinoko <- param[5]
  sigma_Takenoko <- param[6]
  
  repID <- c(data$ID)
  
  # 各回答者の選択確率を予測する
  for(i in 1:236){
    
    # 各選択フェーズでのドロー回数分の選択確率を予測
    P1 <- f_rcdclogit_prob(alpha_Kinoko, 
                           alpha_Takenoko, 
                           beta, 
                           200, 
                           200, 
                           sigma_Kinoko, 
                           sigma_Takenoko, 
                           sigma_price, 
                           tau[(times_draw*(i-1)+1):(times_draw*i),])
    
    P2 <- f_rcdclogit_prob(alpha_Kinoko, 
                           alpha_Takenoko, 
                           beta, 
                           180, 
                           200, 
                           sigma_Kinoko, 
                           sigma_Takenoko, 
                           sigma_price, 
                           tau[(times_draw*(i-1)+1):(times_draw*i),])
    
    P3 <- f_rcdclogit_prob(alpha_Kinoko, 
                           alpha_Takenoko, 
                           beta, 
                           200, 
                           170, 
                           sigma_Kinoko, 
                           sigma_Takenoko, 
                           sigma_price, 
                           tau[(times_draw*(i-1)+1):(times_draw*i),])
    
    P4 <- f_rcdclogit_prob(alpha_Kinoko, 
                           alpha_Takenoko, 
                           beta, 
                           220, 
                           200, 
                           sigma_Kinoko, 
                           sigma_Takenoko, 
                           sigma_price, 
                           tau[(times_draw*(i-1)+1):(times_draw*i),])
    
    P5 <- f_rcdclogit_prob(alpha_Kinoko, 
                           alpha_Takenoko, 
                           beta, 
                           190, 
                           210, 
                           sigma_Kinoko, 
                           sigma_Takenoko, 
                           sigma_price, 
                           tau[(times_draw*(i-1)+1):(times_draw*i),])
    
    # 各設問ごとのドロー回数分の選択確率を格納
    pred[(5*(i-1)*(times_draw)+1):(5*(i-1)*(times_draw)+times_draw),] <- P1
    pred[(5*(i-1)*(times_draw)+times_draw+1):(5*(i-1)*(times_draw)+2*times_draw),] <- P2
    pred[(5*(i-1)*(times_draw)+2*times_draw+1):(5*(i-1)*(times_draw)+3*times_draw),] <- P3
    pred[(5*(i-1)*(times_draw)+3*times_draw+1):(5*(i-1)*(times_draw)+4*times_draw),] <- P4
    pred[(5*(i-1)*(times_draw)+4*times_draw+1):(5*(i-1)*(times_draw)+5*times_draw),] <- P5
    
    # 各回答者の固有のIDをあらかじめ用意しておく。
    ID[(5*(i-1)*(times_draw)+1):(5*(i-1)*(times_draw)+5*times_draw)] <- rep(repID[(5*(i-1)+1)], 5*times_draw)
    
  }
  
  pred<- as.data.frame(pred)
  
  # IDと設問、ドローフェーズ(マージ用)
  pred$ID<- ID
  pred$occasion <- rep(rep(c("Q1", "Q2", "Q3", "Q4", "Q5"), each = times_draw), 236)
  pred$draw_phase <- rep(c(1:times_draw), 5*236)
  
  # 各選択フェーズに対する実際の選択と、予測された選択確率を対応させる
  result <- data %>%
    left_join(pred, by = c("ID","occasion")) %>%
    mutate( choice_prob = case_when( choice == 0 ~ prob_Other,
                                     choice == 1 ~ prob_Kinoko,
                                     choice == 2 ~ prob_Takenoko) ) %>%
    group_by(ID, draw_phase) %>%　# 各回答者の設問1～5への選択の組み合わせが実現する確率を各ドローフェーズごとに算出
    summarise(joint_by_draw = prod(choice_prob),.groups = "drop") %>%
    group_by(ID) %>% # 算出した確率を各回答者ごとに平均する
    summarise(mean_joint_prob=mean(joint_by_draw), .groups = "drop") %>%　
    mutate(log_choice_joint_prob = log(mean_joint_prob)) %>%
    select(log_choice_joint_prob)
  
  
  # 和を取って対数尤度を計算する
  likelihood <- sum(result)
  
  return(likelihood)
  
}



# 消費者属性を入れたランダム係数ロジットにおいて、選択確率を算出する関数
f_rcdclogit_prob_with_atr <- function(alpha_Kinoko, 
                                      alpha_Takenoko, 
                                      beta, 
                                      price_Kinoko, 
                                      price_Takenoko, 
                                      param_atr_Kinoko,
                                      param_atr_Takenoko, 
                                      param_atr_price, 
                                      atr, 
                                      sigma_Kinoko, 
                                      sigma_Takenoko, 
                                      sigma_price, 
                                      tau){
  # きのこ、たけのこから得られるそれぞれの効用の算出
  # sigma_〇〇は各パラメタの標準偏差、tauは標準正規分布に従う乱数を(ドロー回数)行×3列発生させたもの
  alpha_Kinoko_draw <- alpha_Kinoko + sigma_Kinoko*tau[, 1]
  alpha_Takenoko_draw <- alpha_Takenoko + sigma_Takenoko*tau[, 2]
  beta_draw <- beta + sigma_price*tau[, 3]
  util_Kinoko <- alpha_Kinoko_draw + sum(param_atr_Kinoko*atr) - price_Kinoko*(beta_draw + sum(param_atr_price*atr))
  util_Takenoko <- alpha_Takenoko_draw + sum(param_atr_Takenoko*atr) - price_Takenoko*(beta_draw + sum(param_atr_price*atr))
  
  # 選択確率の算出
  prob_Kinoko <- exp(util_Kinoko) / (1 + exp(util_Kinoko) + exp(util_Takenoko))
  prob_Takenoko <- exp(util_Takenoko) / (1 + exp(util_Kinoko) + exp(util_Takenoko))
  prob_Other <- 1 - (prob_Kinoko + prob_Takenoko)
  
  return(cbind(prob_Kinoko, prob_Takenoko, prob_Other))
}



# 消費者属性を入れないランダム係数ロジットにおける対数尤度関数を定義する関数
f_likelihood_rcdclogit_with_atr <- function(param, 
                                            data){
  # 空のデータフレーム
  pred <- matrix(ncol = 3, nrow = 1180*times_draw)
  colnames(pred) <- c("prob_Kinoko", "prob_Takenoko","prob_Other")
  ID <- numeric(236*times_draw)
  atr <- numeric(5)
  beta <- param[1]
  alpha_Kinoko <- param[2]
  alpha_Takenoko <- param[3]
  sigma_price <- param[4]
  sigma_Kinoko <- param[5]
  sigma_Takenoko <- param[6]
  param_atr_Kinoko <- param[c(9, 12, 15, 20, 21)]
  param_atr_Takenoko <- param[c(8, 11, 14, 18, 19)]
  param_atr_price <- param[c(7, 10, 13, 16, 17)]
  repID <- c(data$ID)
  
  
  # 各回答者の選択確率を予測する
  for(i in 1:236){
    gender <- data$gender[(5*(i-1)+1)]
    familyhouse <- data$familyhouse[(5*(i-1)+1)]
    adult <- data$adult[(5*(i-1)+1)]
    kansai <- data$kansai[(5*(i-1)+1)]
    oversea <- data$oversea[(5*(i-1)+1)]
    atr[1:5] <- c(gender, familyhouse, adult, kansai, oversea)
    
    # 各選択フェーズでのドロー回数分の選択確率を予測
    P1 <- f_rcdclogit_prob_with_atr(alpha_Kinoko, 
                                    alpha_Takenoko, 
                                    beta, 
                                    200, 
                                    200, 
                                    param_atr_Kinoko,
                                    param_atr_Takenoko, 
                                    param_atr_price, 
                                    atr, 
                                    sigma_Kinoko, 
                                    sigma_Takenoko, 
                                    sigma_price, 
                                    tau[(times_draw*(i-1)+1):(times_draw*i),])
    
    P2 <- f_rcdclogit_prob_with_atr(alpha_Kinoko, 
                                    alpha_Takenoko, 
                                    beta, 
                                    180, 
                                    200, 
                                    param_atr_Kinoko,
                                    param_atr_Takenoko, 
                                    param_atr_price, 
                                    atr, 
                                    sigma_Kinoko, 
                                    sigma_Takenoko, 
                                    sigma_price, 
                                    tau[(times_draw*(i-1)+1):(times_draw*i),])
    
    P3 <- f_rcdclogit_prob_with_atr(alpha_Kinoko, 
                                    alpha_Takenoko, 
                                    beta, 
                                    200, 
                                    170, 
                                    param_atr_Kinoko,
                                    param_atr_Takenoko, 
                                    param_atr_price, 
                                    atr, 
                                    sigma_Kinoko, 
                                    sigma_Takenoko, 
                                    sigma_price, 
                                    tau[(times_draw*(i-1)+1):(times_draw*i),])
    
    P4 <- f_rcdclogit_prob_with_atr(alpha_Kinoko, 
                                    alpha_Takenoko, 
                                    beta, 
                                    220, 
                                    200, 
                                    param_atr_Kinoko,
                                    param_atr_Takenoko, 
                                    param_atr_price, 
                                    atr, 
                                    sigma_Kinoko, 
                                    sigma_Takenoko, 
                                    sigma_price, 
                                    tau[(times_draw*(i-1)+1):(times_draw*i), ])
    
    P5 <- f_rcdclogit_prob_with_atr(alpha_Kinoko, 
                                    alpha_Takenoko, 
                                    beta, 
                                    190, 
                                    210, 
                                    param_atr_Kinoko,
                                    param_atr_Takenoko, 
                                    param_atr_price, 
                                    atr, 
                                    sigma_Kinoko, 
                                    sigma_Takenoko, 
                                    sigma_price, 
                                    tau[(times_draw*(i-1)+1):(times_draw*i),])
    
    # 各設問ごとのドロー回数分の選択確率を格納
    pred[(5*(i-1)*(times_draw)+1):(5*(i-1)*(times_draw)+times_draw),] <- P1
    pred[(5*(i-1)*(times_draw)+times_draw+1):(5*(i-1)*(times_draw)+2*times_draw),] <- P2
    pred[(5*(i-1)*(times_draw)+2*times_draw+1):(5*(i-1)*(times_draw)+3*times_draw),] <- P3
    pred[(5*(i-1)*(times_draw)+3*times_draw+1):(5*(i-1)*(times_draw)+4*times_draw),] <- P4
    pred[(5*(i-1)*(times_draw)+4*times_draw+1):(5*(i-1)*(times_draw)+5*times_draw),] <- P5
    
    # 各回答者の固有のIDをあらかじめ用意しておく。
    ID[(5*(i-1)*(times_draw)+1):(5*(i-1)*(times_draw)+5*times_draw)] <- rep(repID[(5*(i-1)+1)], 5*times_draw)
    
  }
  
  pred<- as.data.frame(pred)
  
  # IDと設問、ドローフェーズ(マージ用)
  pred$ID<- ID
  pred$occasion <- rep(rep(c("Q1", "Q2", "Q3", "Q4", "Q5"), each = times_draw), 236)
  pred$draw_phase <- rep(c(1:times_draw), 5*236)
  
  # 各選択フェーズに対する実際の選択と、予測された選択確率を対応させる
  result <- data %>%
    left_join(pred, by = c("ID","occasion")) %>%
    mutate( choice_prob = case_when( choice == 0 ~ prob_Other,
                                     choice == 1 ~ prob_Kinoko,
                                     choice == 2 ~ prob_Takenoko) ) %>%
    group_by(ID, draw_phase) %>%　# 各回答者の設問1～5への選択の組み合わせが実現する確率を各ドローフェーズごとに算出。
    summarise(joint_by_draw = prod(choice_prob),.groups = "drop") %>%
    group_by(ID) %>% # 算出した確率を各回答者ごとに平均する。
    summarise(mean_joint_prob=mean(joint_by_draw), .groups = "drop") %>%　
    mutate(log_choice_joint_prob = log(mean_joint_prob)) %>%
    select(log_choice_joint_prob)
  
  
  # 和を取って対数尤度を計算する
  likelihood <- sum(result)
  
  return(likelihood)
  
}
