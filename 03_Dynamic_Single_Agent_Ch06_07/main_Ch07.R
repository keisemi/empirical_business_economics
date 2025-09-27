
# パッケージを読み込む
library(tidyverse)
library(skimr)
library(fixest)
library(numDeriv)
library(kableExtra)
library(here)

source(here("03_Dynamic_Single_Agent_Ch06_07/function.R"))

# データの読み込み ----

# 第6章で生成した自動車の買い替え行動に関するデータを引き続き用いる
data_gen <- readr::read_csv(here('03_Dynamic_Single_Agent_Ch06_07/intermediate/Chap6_data.csv'))

data_gen %>% head(3)

# Stateの作成 ----

# 価格の状態変数
price_states <- seq(2, 2.5, by = 0.1)

# 走行距離の状態変数
mileage_states <- seq(0, 0.1, by = 0.005)

# 価格の状態変数の数
num_price_states <- length(price_states)

# 走行距離の状態変数の数
num_mileage_states <- length(mileage_states)

# 状態変数は価格と走行距離の状態変数のペア
# したがって状態変数の数は価格の状態変数の数と走行距離の状態変数の数の積となる
num_states <- num_price_states * num_mileage_states

# 価格、走行距離の状態変数の組み合わせ(p,m)を1つのデータフレームで表す
state_df <- dplyr::tibble(
    # stateにidを番号付ける
    state_id = 1:num_states,
    # 順番は (p,m) = (2000,0), (2100,0), ..., (2500,0), (2000,5), (2100,5), ...
    price_id = rep(1:num_price_states, times = num_mileage_states),
    mileage_id = rep(1:num_mileage_states, each = num_price_states),
    price = rep(price_states, times = num_mileage_states),
    mileage = rep(mileage_states, each = num_price_states)
  )

# 下3行を表示
state_df %>% tail(3)

# 分析のためにデータを加工する
data_gen <- data_gen %>% 
  dplyr::group_by(consumer) %>% 
  # 遷移行列の推定で使うため、ラグ変数を追加
  dplyr::mutate(lag_price_id = lag(price_id),
                lag_mileage_id = lag(mileage_id),
                lag_action = lag(action)) %>% 
  dplyr::ungroup() 

# 記述統計の確認 ----

# 生成したデータの要約統計（Ch06と同じ）
data_gen %>% 
  dplyr::select(price, mileage, action) %>%
  skimr::skim() %>% 
  skimr::yank('numeric') %>% 
  dplyr::select(skim_variable, mean, sd, p0, p100) %>% { 
    txt <- capture.output(knitr::kable(., format = "simple"))
    writeLines(txt, here("03_Dynamic_Single_Agent_Ch06_07/output/tbl_descr_for_Ch07.txt"), useBytes = TRUE)
  }

# 二段階推定方法の実践 ----

# 推定に必要なパラメタの設定

# 時間割引率
beta <- 0.95

# オイラー定数
Euler_const <- - digamma(1)

# 消費者は車を購入するかどうか決めるため選択肢の数は2つ
num_choice <- 2

## Step 1: CCP の推定、及び State transition の推定 ----

# 走行距離の遷移行列を推定

# それぞれの確率が実現した観察の数を数える
num_cond_obs_mileage <- data_gen %>% 
  # 1期目は推定に使えないため落とす
  dplyr::filter(period != min(data_gen$period)) %>% 
  # t期の走行距離、t+1期の走行距離、t期の購買ごとにグループ化して、観察数を数える
  dplyr::group_by(lag_mileage_id, mileage_id, lag_action) %>% 
  dplyr::summarize(num_cond_obs = n(),
                   .groups = 'drop') %>% 
  # 確率ごとに名前を割り当てる
  dplyr::mutate(
    cond_obs_mileage = case_when(
      # 1 - kappa_1 - kappa_2 の場合
      (
        (lag_action == 0 &
           between(lag_mileage_id, 1, 20) &
           (lag_mileage_id == mileage_id)) |
          (lag_action == 1 & 
             mileage_id == 1)
      ) ~ 'cond_obs_mileage1',
      # kappa_1 の場合
      (
        (lag_action == 0 &
           between(lag_mileage_id, 1, 19) &
           (lag_mileage_id == mileage_id - 1)) |
          (lag_action == 1 & 
             mileage_id == 2)
      ) ~ 'cond_obs_mileage2',
      # kappa_2 の場合
      (
        (lag_action == 0 &
           between(lag_mileage_id, 1, 19) &
           (lag_mileage_id == mileage_id - 2)) |
          (lag_action == 1 & 
             mileage_id == 3)
      ) ~ 'cond_obs_mileage3',
      # kappa_1 + kappa_2 の場合
      (
        lag_action == 0 &
          lag_mileage_id == 20 &
          mileage_id == 21
      ) ~ 'cond_obs_mileage4',
      TRUE ~ 'other'
    )) %>% 
  # 'other' は推定には使わないため落とす
  dplyr::filter(cond_obs_mileage != 'other') %>% 
  # 確率ごとにグループ化し、再度、観察の数を数える
  dplyr::group_by(cond_obs_mileage) %>% 
  dplyr::summarize(num_cond_obs = as.numeric(sum(num_cond_obs)),
                   .groups = 'drop') %>% 
  dplyr::select(num_cond_obs) %>% 
  as.matrix() 

# 最尤法の解析解により推定値を求める
kappa_est <- c()

kappa_est[1] <- 
  (num_cond_obs_mileage[2] * 
     (num_cond_obs_mileage[2] + num_cond_obs_mileage[3] + num_cond_obs_mileage[4])) /
  ((num_cond_obs_mileage[2] + num_cond_obs_mileage[3]) * 
     (num_cond_obs_mileage[1] + num_cond_obs_mileage[2] + num_cond_obs_mileage[3] + num_cond_obs_mileage[4]))

kappa_est[2] <- 
  (num_cond_obs_mileage[3] * 
     (num_cond_obs_mileage[2] + num_cond_obs_mileage[3] + num_cond_obs_mileage[4])) /
  ((num_cond_obs_mileage[2] + num_cond_obs_mileage[3]) * 
     (num_cond_obs_mileage[1] + num_cond_obs_mileage[2] + num_cond_obs_mileage[3] + num_cond_obs_mileage[4]))

# 価格の遷移行列を推定

# それぞれの確率が実現した観察の数を数える
num_cond_obs_price <- data_gen %>% 
  # 1期目は推定に使えないため落とす
  dplyr::filter(period != min(data_gen$period)) %>% 
  # t期の価格、t+1期の価格ごとにグループ化して、観察数を数える
  dplyr::group_by(lag_price_id, price_id) %>% 
  dplyr::summarize(num_cond_obs = n(),
                   .groups = 'drop') %>% 
  # 観察数を行列（num_price_states行の正方行列）に変換
  # price_id (t+1期の価格) を横に広げる
  tidyr::pivot_wider(names_from = 'price_id',
                     values_from = 'num_cond_obs') %>%
  dplyr::select(!lag_price_id) %>% 
  as.matrix()


# 最尤法の解析解により推定値を求める
lambda_est_mat <- num_cond_obs_price / rowSums(num_cond_obs_price)
lambda_est_mat

lambda_est <- as.vector(t(lambda_est_mat))[c(-1, -8, -15, -22, -29, -36)]

# 推定された遷移行列を取得（本誌のG(i)に対応）
# 配列の第三次元が購買を意味しており、i=0,1という順番で並んでいる。
G <- array(c(gen_mileage_trans(kappa_est)[,,1] %x% gen_price_trans(lambda_est), 
             gen_mileage_trans(kappa_est)[,,2] %x% gen_price_trans(lambda_est)),
           dim = c(num_states, num_states, num_choice))

# CCPの推定 

# 単純なロジットモデルを推定する
logit_model <- fixest::feglm(action ~ price + price^2 + mileage + mileage^2 + price:mileage,
                             data_gen,
                             family = binomial('logit'))

# 推定値に基づいてCCPを予測する。列はi=0,1という順番。
CCP_1st <- cbind(1 - predict(logit_model, state_df),
                 predict(logit_model, state_df))

## Step 2: パラメタの推定 ----

# 以下のコードでは、二段階目の方法として紹介のあった (1) 行列形式によるインバージョン、と (2) 有限依存性 (finite dependence) アプローチの２つの手法を実装する

### 行列形式によるインバージョン ----

# パラメタの真の値
theta_true <- c(theta_c = 15, theta_p = 6)

start_time_mat_inv <- proc.time()

# 最適化
mat_inv_opt_mat_inv <- optim(theta_true,
                             likelihood_fun,
                             CCP = CCP_1st,
                             df = data_gen, 
                             beta = beta,
                             G = G, 
                             state_df = state_df, 
                             policy_operator = policy_operator_mat_inv,
                             control = list(fnscale = -1), 
                             method = 'Nelder-Mead')

end_time_mat_inv <- proc.time()
run_time_mat_inv <- (end_time_mat_inv - start_time_mat_inv)[[3]]
print(stringr::str_c('Runtime: ', round(run_time_mat_inv, 2)))

# 結果の確認
theta_mat_inv <- mat_inv_opt_mat_inv$par
theta_mat_inv


# 標準誤差を推定する
hessian <- numDeriv::hessian(func = likelihood_fun,
                             x = theta_mat_inv,
                             CCP = CCP_1st, 
                             df = data_gen, 
                             beta = beta, 
                             G = G, 
                             state_df = state_df, 
                             policy_operator = policy_operator_mat_inv)

theta_se_mat_inv <- sqrt(diag(solve(-hessian)))

dplyr::tibble(theta_mat_inv, theta_se_mat_inv) %>% { 
  txt <- capture.output(knitr::kable(., format = "simple"))
  writeLines(txt, here("03_Dynamic_Single_Agent_Ch06_07/output/tbl_theta_est_mat.txt"), useBytes = TRUE)
  }


### 有限依存性 (finite dependence) アプローチ ----

start_time_finite_dep <- proc.time()

# 最適化
finite_dep_opt <- optim(theta_true, 
                        likelihood_fun,
                        CCP = CCP_1st, 
                        df = data_gen, 
                        beta = beta, 
                        G = G, 
                        state_df = state_df, 
                        policy_operator = policy_operator_finite_dep,
                        control = list(fnscale = -1), 
                        method = 'Nelder-Mead')

end_time_finite_dep <- proc.time()
run_time_finite_dep <- (end_time_finite_dep - start_time_finite_dep)[[3]]
print(stringr::str_c('Runtime: ', round(run_time_finite_dep, 2)))

# 結果の確認
theta_finite_dep <- finite_dep_opt$par
theta_finite_dep

# 標準誤差を推定する
hessian <- numDeriv::hessian(func = likelihood_fun,
                             x = theta_finite_dep,
                             CCP = CCP_1st, 
                             df = data_gen, 
                             beta = beta, 
                             G = G, 
                             state_df = state_df, 
                             policy_operator = policy_operator_finite_dep)

theta_se_finite_dep <- sqrt(diag(solve(-hessian)))

dplyr::tibble(theta_finite_dep, theta_se_finite_dep) %>% { 
  txt <- capture.output(knitr::kable(., format = "simple"))
  writeLines(txt, here("03_Dynamic_Single_Agent_Ch06_07/output/tbl_theta_est_finite_dep.txt"), useBytes = TRUE)
}

# 推定方法の比較 ----

# 第6章の復習：入れ子不動点アルゴリズム

start_time_nfxp <- proc.time()

# 最適化
nfxp_opt <- optim(theta_true, 
                  likelihood_fun_nfxp,
                  df = data_gen, 
                  beta = beta, 
                  G = G, 
                  state_df = state_df, 
                  control = list(fnscale = -1), 
                  method = 'Nelder-Mead')

end_time_nfxp <- proc.time()

run_time_nfxp <- (end_time_nfxp - start_time_nfxp)[[3]]

print(stringr::str_c('Runtime: ', round(run_time_nfxp, 2)))

theta_nfxp <- nfxp_opt$par
theta_nfxp

# 標準誤差を推定する
hessian <- numDeriv::hessian(func = likelihood_fun_nfxp, 
                             x = theta_nfxp,
                             df = data_gen, 
                             beta = beta, 
                             G = G, 
                             state_df = state_df)

theta_se_nfxp <- sqrt(diag(solve(-hessian)))

dplyr::tibble(theta_nfxp, theta_se_nfxp) %>% { 
    txt <- capture.output(knitr::kable(., format = "simple"))
    writeLines(txt, here("03_Dynamic_Single_Agent_Ch06_07/output/tbl_theta_est_nfxp_ch7.txt"), useBytes = TRUE)
  }

# 推定方法の比較表
dplyr::tibble(
  algorithm = c('NFXP', 'Matrix', 'Finite', 'True'),
  theta_c = c(theta_nfxp[1], theta_mat_inv[1], theta_finite_dep[1], theta_true[1]),
  theta_se_c = c(theta_se_nfxp[1], theta_se_mat_inv[1], theta_se_finite_dep[1], NA),
  theta_p = c(theta_nfxp[2], theta_mat_inv[2], theta_finite_dep[2], theta_true[2]),
  theta_se_p = c(theta_se_nfxp[2], theta_se_mat_inv[2], theta_se_finite_dep[2], NA),
  run_time = c(run_time_nfxp, run_time_mat_inv, run_time_finite_dep, NA)
) %>% { 
  txt <- capture.output(knitr::kable(., format = "simple"))
  writeLines(txt, here("03_Dynamic_Single_Agent_Ch06_07/output/tab7_1_compare_algo.txt"), useBytes = TRUE)
}

# 反実仮想分析1: EDLP(Every Day Low Price) ----

# 以下の二つのシナリオについてそれぞれCCPを計算し、結果を比較する。
# シナリオ1: 推定した際と全く同じ状況（以降、ベースラインと呼ぶ）
# シナリオ2: 価格が最も低い価格に永続的に固定される状況  

# 以降のシミュレーションで得られたCCPを要素とするリストを作成
CCP_list <- list()

## シナリオ1) 推定した際と全く同じ状況 ----
CCP_list[['Baseline']] <- policy_operator_nfxp(theta_nfxp, beta, G, state_df)

# シナリオ1における、各価格の購買確率を計算
result_df_edlp <- data_gen %>% 
  # データよりそれぞれの状態変数が観察された数を数える
  dplyr::group_by(state_id, price_id, price) %>% 
  dplyr::summarize(num_obs = n(),
                   .groups = 'drop') %>% 
  dplyr::arrange(state_id) %>% 
  # シナリオ1で得られたCCPを新たな列に加える
  dplyr::mutate(prob_buy_baseline = CCP_list[['Baseline']][, 2]) %>% 
  # 走行距離について、購買確率の加重平均を取る（走行距離の頻度について積分する）
  dplyr::group_by(price_id, price) %>% 
  dplyr::summarize(prob_buy_baseline = weighted.mean(prob_buy_baseline, num_obs), 
                   .groups = 'drop')

## シナリオ2) 価格が最も低い価格に永続的に固定される状況 ----

# このシナリオでは価格の変動がないため、走行距離についてのみの遷移行列がこのモデルの遷移行列となる
G_fixed_price <- gen_mileage_trans(kappa_est)

# 価格を固定した場合のCCPを計算
# 後ほどの需要と収入の計算で使うため、価格を固定した場合のCCPをすべての価格に対して計算
for (fixed_price in price_states) {
  # 価格が固定された場合の、状態変数を示すデータフレームを作成
  state_df_fixed_price <- state_df %>% 
    dplyr::filter(price == fixed_price) %>% 
    dplyr::arrange(mileage_id)
  
  # 遷移行列と状態変数を示すデータフレームが変わっていることに注意してCCPを計算
  CCP_list[[stringr::str_c('edlp', fixed_price * 100)]] <- 
    policy_operator_nfxp(theta = theta_nfxp, 
                         beta = beta, 
                         G = G_fixed_price, 
                         state_df = state_df_fixed_price)
}

# 200万円に価格を固定した場合の購買確率を計算
prob_buy_edlp200 <- data_gen %>% 
  # データよりそれぞれの状態変数が観察された数を数える
  dplyr::group_by(mileage_id, mileage) %>% 
  dplyr::summarize(num_obs = n(),
                   .groups = 'drop') %>% 
  dplyr::arrange(mileage_id) %>% 
  # 価格が固定された場合のCCPを新たな列に加える
  dplyr::mutate(prob_buy = CCP_list[['edlp200']][, 2]) %>% 
  # 走行距離について、購買確率の加重平均を取る（走行距離の頻度について積分する）
  dplyr::summarize(prob_buy = weighted.mean(prob_buy, num_obs)) %>% 
  as.matrix()


# 価格ごとの購買確率を図示する
result_df_edlp %>% 
  ggplot2::ggplot(aes(x = price * 100, y = prob_buy_baseline)) + 
  ggplot2::geom_bar(stat = 'identity') +
  ggplot2::geom_hline(yintercept = prob_buy_edlp200,
                      color = 'black', linewidth = 1.0,
                      linetype = 'dashed') +
  ggplot2::annotate('text', 
                    x = 250,
                    y = prob_buy_edlp200 + 0.005, 
                    color= 'black', 
                    label = str_c('EDLP: ', round(prob_buy_edlp200, 5))) +
  ggplot2::scale_x_continuous(labels = scales::label_number(scale = 1),
                              breaks = seq(200, 250, 10)) +
  ggplot2::scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.title.y = element_text(angle = 0, vjust = 0.5),
    text = element_text(family = 'HiraKakuPro-W3'),
    legend.position = 'none',
    panel.grid.major.y = element_line(color = 'grey',
                                      linewidth = 0.3, linetype = 'dotted'),
    panel.grid.minor.y = element_line(color = 'grey',
                                      linewidth = 0.3, linetype = 'dotted'),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = 'black', linewidth = 0.3)
  ) +
  ggplot2::labs(y = str_wrap('購 買 確 率', width = 1), x = '価格（万円）')

ggsave(here("03_Dynamic_Single_Agent_Ch06_07/output/fig7_1_Counter1_EDLP.pdf"), device = cairo_pdf)


# 反実仮想分析2: 永続的・一時的な値下げ ----

#以下の３つのシナリオを比較することで、値下げの効果を分析する。
# シナリオ1: ベースライン
# シナリオ2: 10万円の永続的な値下げ
# シナリオ3: 10万円の一時的な値下げ

## 値下げが起きた時のCCPの計算

# まず、それぞれのシナリオにおけるCCPを計算する。

##シナリオ2) 10万円の永続的な値下げ ----
# priceを0.1 (10万円)下げた、state_dfを用意し、CCPを再計算する
state_df_discount <- state_df %>% 
  dplyr::mutate(price = price - 0.1)

CCP_list[['Permanent']] <- policy_operator_nfxp(theta_nfxp, beta, G, state_df_discount)

## シナリオ3) 10万円の一時的な値下げ ----

# 値下げが起きたときの今期の効用を計算
U_discount <- flow_utility(theta_nfxp, state_df_discount)
# 元の価格設定における事前の価値関数を計算
V <- contraction_nfxp(theta_nfxp, beta, G, state_df)
# 選択肢ごとの価値関数を計算
CV_temporary <- U_discount + beta * cbind(G[,, 1] %*% V, G[,, 2] %*% V)
# CCPを計算
CCP_list[['Temporary']] <-exp(CV_temporary) / rowSums(exp(CV_temporary))


# 三つのシナリオについて、元々の価格が220万円の場合の、走行距離ごとの購入確率を図示
state_df %>% 
  # 得られたCCPを一つのdataframeにまとめる
  dplyr::mutate(ProbBaseline = CCP_list[['Baseline']][, 2],
                ProbPermanent = CCP_list[['Permanent']][, 2],
                ProbTemporary = CCP_list[['Temporary']][, 2]) %>% 
  # 220万円での購買確率に絞る
  dplyr::filter(price == 2.2) %>% 
  dplyr::select(mileage, ProbBaseline, 
                ProbPermanent, ProbTemporary) %>% 
  tidyr::pivot_longer(cols = - mileage) %>% 
  dplyr::rename(scenario = name, prob_buy = value) %>% 
  ggplot2::ggplot(aes(x = mileage * 100, y = prob_buy, color = scenario)) +
  ggplot2::geom_line() +
  ggplot2::geom_point() +
  ggplot2::theme_minimal()+
  ggplot2::scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  ggplot2::geom_hline(yintercept = 0,
                      linetype = 'dashed',
                      color = 'black', size = 0.3)+
  ggplot2::theme(
    axis.title.y = element_text(angle = 0, vjust = 0.5),
    text = element_text(family = 'HiraKakuPro-W3'),
    panel.grid.major.y = element_line(color = 'grey',
                                      linewidth = 0.3, linetype = 'dotted'),
    panel.grid.minor.y = element_line(color = 'grey',
                                      linewidth = 0.3, linetype = 'dotted'),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = 'black', linewidth = 0.2)
  ) +
  scale_x_continuous(breaks = seq(0, 10, by = 2)) +
  ggplot2::labs(y = str_wrap('購 買 確 率', width = 1), x = '走行距離（万km）') +
  ggplot2::scale_color_grey(name = 'シナリオ', 
                            labels = c('ベースライン','永続的', '一時的'))

ggsave(here("03_Dynamic_Single_Agent_Ch06_07/output/Counter2_CCP_mile.pdf"), device = cairo_pdf)


## 消費者の分布のシミュレーション（値下げ）

# データから状態変数毎の消費者の分布を計算
# state_idの順番に並んだ、(1 * num_states) の行ベクトル
consumer_dist_obs <- data_gen %>% 
  dplyr::group_by(state_id, 
                  price_id, 
                  mileage_id) %>% 
  dplyr::summarize(num_obs = n(),
                   .groups = 'drop') %>% 
  dplyr::mutate(consumer_dist_obs = num_obs / sum(num_obs)) %>% 
  dplyr::arrange(state_id) %>% 
  dplyr::select(consumer_dist_obs) %>% 
  as.matrix() %>% 
  t()

# 購買を考慮した遷移行列を作成
discount_scenerio_vec <- c('Baseline','Permanent','Temporary')

# CCPを考慮した遷移行列を要素とするリストを作成
G_CCP_list <- list()

# それぞれの値下げシナリオに応じた、購買を考慮した遷移行列を作成
for (scenario in discount_scenerio_vec) {
  G_CCP_list[[scenario]] <- 
    (CCP_list[[scenario]][,1] %*% matrix(1, 1, num_states)) * G[,,1] +
    (CCP_list[[scenario]][,2] %*% matrix(1, 1, num_states)) * G[,,2]
}

# 需要と収入の計算に使うため、EDLPを行ったときの購買を考慮した遷移行列も作成
for (scenario in stringr::str_c('edlp', price_states * 100)) {
  G_CCP_list[[scenario]] <- 
    (CCP_list[[scenario]][,1] %*% matrix(1, 1, nrow(G_fixed_price[,,1]))) * G_fixed_price[,,1] +
    (CCP_list[[scenario]][,2] %*% matrix(1, 1, nrow(G_fixed_price[,,1]))) * G_fixed_price[,,2]
}

# シミュレーションする消費者数と期間を定義
num_consumer_sim <- 1000
num_period_sim <- 10

# 値下げによるシミュレーション結果を表すデータフレームを前もって作成
discount_sim_df <- dplyr::tibble(period = 1:num_period_sim)


# それぞれのシナリオについてシミュレーションを行う
for (scenario in discount_scenerio_vec) {
  
  # 消費者の分布の推移を表す行列を前もって作成
  consumer_dist_sim <- matrix(0, num_period_sim, num_states)
  # 初期の分布をデータ上の分布とする
  consumer_dist_sim[1,] <- consumer_dist_obs
  
  # 購買割合、需要、収入の空の列を作成
  discount_sim_df <- discount_sim_df %>% 
    dplyr::mutate(
      prob_buy_sim = NA,
      demand = NA,
      revenue = NA
    )
  
  # t期における消費者の分布、購買割合、需要、収入を計算
  for (t in 1:num_period_sim) {
    # 初期、かつ、一時的な値下げの場合にのみ、一時的な値下げのCCPを参照する
    if (t == 1 & scenario == 'Temporary') {
      # 分布の推移
      consumer_dist_sim[t+1,] <- consumer_dist_sim[t,] %*% G_CCP_list[['Temporary']]
      # 購買割合の計算
      discount_sim_df$prob_buy_sim[t] <- consumer_dist_sim[t,] %*% CCP_list[['Temporary']][, 2]
      # 需要量の計算
      discount_sim_df$demand[t] <-discount_sim_df$prob_buy_sim[t] * num_consumer_sim
      # 収入の計算 
      discount_sim_df$revenue[t] <- state_df %>% 
        dplyr::mutate(consumer_dist = consumer_dist_sim[t,],
                      CCP = CCP_list[['Temporary']][, 2],
                      revenue = (price - 0.1) * consumer_dist * CCP * num_consumer_sim) %>% 
        dplyr::summarize(revenue = sum(revenue)) %>% 
        as.matrix()
    } else {
      # 実際に消費者が従うCCPのシナリオ名を得る
      # 注意: 一時的な値下げの場合、t >= 2においてはベースラインのCCPを参照する
      scenario_current <- 
        case_when(scenario %in% c('Baseline', 'Permanent') ~ scenario,
                  scenario == 'Temporary' ~ 'Baseline')
      # 分布の推移
      # 注意: 最終期は計算する必要がない
      if (t != num_period_sim) {
        consumer_dist_sim[t + 1,] <- consumer_dist_sim[t,] %*% G_CCP_list[[scenario_current]]
      }
      # 購買割合の計算
      discount_sim_df$prob_buy_sim[t] <- consumer_dist_sim[t,] %*% CCP_list[[scenario_current]][, 2]
      # 需要量の計算
      discount_sim_df$demand[t] <- discount_sim_df$prob_buy_sim[t] * num_consumer_sim
      # 収入の計算
      # 値下げをしている場合のみ、値下げ額を取得
      discount <- 
        case_when(scenario %in% c('Baseline', 'Temporary') ~ 0,
                  scenario == 'Permanent' ~ 0.1)
      discount_sim_df$revenue[t] <- state_df %>% 
        dplyr::mutate(consumer_dist = consumer_dist_sim[t,],
                      CCP = CCP_list[[scenario_current]][, 2],
                      revenue = (price - discount) * consumer_dist * CCP * num_consumer_sim) %>% 
        dplyr::summarize(revenue = sum(revenue)) %>% 
        as.matrix()
    }
  }
  # 購買割合、需要、収入の変数名をシナリオに応じて変更
  discount_sim_df <- discount_sim_df %>% 
    dplyr::rename(!!str_c('prob_buy_sim_', scenario) := 'prob_buy_sim',
                  !!str_c('demand_', scenario) := 'demand',
                  !!str_c('revenue_', scenario) := 'revenue')
}


# 値下げの効果

# ベースラインを基準とした、購買割合の変化を計算
discount_sim_df %>%
  dplyr::mutate(PermanentChange = prob_buy_sim_Permanent - prob_buy_sim_Baseline,
                TemporaryChange = prob_buy_sim_Temporary - prob_buy_sim_Baseline) %>%
  dplyr::select(period, PermanentChange, TemporaryChange) %>%
  tidyr::pivot_longer(cols = -period) %>%
  dplyr::rename(scenario = name, prob_buy = value) %>%
  ggplot2::ggplot(aes(x = period, y = prob_buy, color = scenario)) +
  geom_line() +
  geom_point() +
  ggplot2::theme_minimal() +
  ggplot2::scale_x_continuous(breaks = seq(0, num_period_sim, by = 5)) +
  ggplot2::geom_hline(yintercept = 0,
                      linetype = 'dashed',
                      color = 'black', size = 0.3) +
  ggplot2::theme(
    axis.title.y = element_text(angle = 0, vjust = 0.5),
    text = element_text(family = 'HiraKakuPro-W3'),
    panel.grid.major.y = element_line(color = 'grey',
                                      linewidth = 0.3, linetype = 'dotted'),
    panel.grid.minor.y = element_line(color = 'grey',
                                      linewidth = 0.3, linetype = 'dotted'),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = 'black', linewidth = 0.2)
  ) +
  scale_x_continuous(breaks = seq(0, 10, by = 2)) +
  ggplot2::labs(y = str_wrap('購 買 確 率', width = 1), x = '期') +
  ggplot2::scale_color_grey(name = 'シナリオ', labels = c('永続的', '一時的'))

ggsave(here("03_Dynamic_Single_Agent_Ch06_07/output/fig7_2_Counter2_sim_choice_change.pdf"), device = cairo_pdf)

# 需要と収入の計算 ----

# 今回の反実仮想で扱った全てのシナリオ（EDLPと値下げ）の需要と収入を計算する。
# そのために、まずEDLPの場合の消費者の分布をシミュレートする。

## 消費者の分布のシミュレーション（EDLP）----

# データから状態変数毎の消費者の分布を計算。
# mileage_idの順番に並んだ、(1 * 21) の行ベクトル。
consumer_dist_obs_edlp <- data_gen %>% 
  dplyr::group_by(mileage_id) %>% 
  dplyr::summarize(num_obs = n(),
                   .groups = 'drop') %>% 
  dplyr::mutate(consumer_dist_obs_edlp = num_obs / sum(num_obs)) %>% 
  dplyr::arrange(mileage_id) %>% 
  dplyr::select(consumer_dist_obs_edlp) %>% 
  as.matrix() %>% 
  t()

# EDLPのシミュレーション結果を表すデータフレームを前もって作成
edlp_sim_df <- dplyr::tibble(period = 1:num_period_sim)

# すべての価格についてシミュレーションを行う
for (fixed_price in price_states) {
  scenario_edlp <- stringr::str_c('edlp', fixed_price * 100)
  
  # 消費者の分布の推移を表す行列を前もって作成
  consumer_dist_sim <- matrix(0, num_period_sim, length(consumer_dist_obs_edlp))
  # 初期の分布をデータ上の分布とする
  consumer_dist_sim[1,] <- consumer_dist_obs_edlp
  
  # 購買割合、需要、収入の空の列を作成
  edlp_sim_df <- edlp_sim_df %>% 
    dplyr::mutate(
      prob_buy_sim = NA,
      demand = NA,
      revenue = NA
    )
  
  # t期における消費者の分布、購買割合、需要、収入を計算
  for (t in 1:num_period_sim) {
    # 分布の推移
    # 注意: 最終期は計算する必要がない
    if (t != num_period_sim){
      consumer_dist_sim[t + 1,] <- consumer_dist_sim[t,] %*% G_CCP_list[[scenario_edlp]]
    }
    # 購買割合の計算
    edlp_sim_df$prob_buy_sim[t] <- consumer_dist_sim[t,] %*% CCP_list[[scenario_edlp]][,2]
    # 需要量の計算
    edlp_sim_df$demand[t] <- edlp_sim_df$prob_buy_sim[t] * num_consumer_sim
    # 収入の計算
    edlp_sim_df$revenue[t] <- state_df %>% 
      dplyr::filter(price == fixed_price) %>% 
      dplyr::arrange(mileage_id) %>% 
      dplyr::mutate(consumer_dist = consumer_dist_sim[t,],
                    CCP = CCP_list[[scenario_edlp]][, 2],
                    revenue = fixed_price * consumer_dist * CCP * num_consumer_sim) %>% 
      dplyr::summarize(revenue = sum(revenue)) %>% 
      as.matrix()
  }
  # 購買割合、需要、収入の変数名をシナリオに応じて変更
  edlp_sim_df <- edlp_sim_df %>% 
    dplyr::rename(!!str_c('prob_buy_sim_', scenario_edlp) := 'prob_buy_sim',
                  !!str_c('demand_', scenario_edlp) := 'demand',
                  !!str_c('revenue_', scenario_edlp) := 'revenue')
}

# 結果を図示するためにすべてのシミュレーション結果を統合
sim_df <- discount_sim_df %>% 
  dplyr::left_join(edlp_sim_df, by = 'period')

## シナリオ毎の需要と収入の比較 ----

### EDLPの需要と収入の推移 ----

# 需要の推移を図示
sim_df %>% 
  dplyr::select(period, stringr::str_c('demand_edlp', price_states * 100)) %>% 
  tidyr::pivot_longer(cols = - period) %>% 
  dplyr::rename(scenario = name, demand = value) %>% 
  ggplot2::ggplot(aes(x = period, y = demand, color = scenario)) +
  ggplot2::geom_line() +
  ggplot2::geom_point() +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    axis.title.y = element_text(angle = 0, vjust = 0.5),
    text = element_text(family = 'HiraKakuPro-W3'),
    panel.grid.major.y = element_line(color = 'grey',
                                      linewidth = 0.3, linetype = 'dotted'),
    panel.grid.minor.y = element_line(color = 'grey',
                                      linewidth = 0.3, linetype = 'dotted'),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = 'black', linewidth = 0.2)
  ) +
  scale_x_continuous(breaks = seq(0, 10, by = 2)) +
  ggplot2::labs(y = str_wrap('需 要 量', width = 1), x = '期') +
  ggplot2::scale_color_grey(name = 'シナリオ',
                            labels = c('EDLP200万円', 'EDLP210万円',
                                       'EDLP220万円', 'EDLP230万円',
                                       'EDLP240万円', 'EDLP250万円'))

ggsave(here("03_Dynamic_Single_Agent_Ch06_07/output/Counter2_sim_demand_EDLP.pdf"), device = cairo_pdf)

# 収入(万円)の推移
sim_df %>% 
  dplyr::select(period, stringr::str_c('revenue_edlp', price_states * 100)) %>% 
  tidyr::pivot_longer(cols = - period) %>% 
  dplyr::rename(scenario = name, revenue = value) %>% 
  ggplot2::ggplot(aes(x = period, y = revenue * 100, color = scenario)) +
  ggplot2::geom_line() +
  ggplot2::geom_point() +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    axis.title.y = element_text(angle = 0, vjust = 0.5),
    text = element_text(family = 'HiraKakuPro-W3'),
    panel.grid.major.y = element_line(color = 'grey',
                                      linewidth = 0.3, linetype = 'dotted'),
    panel.grid.minor.y = element_line(color = 'grey',
                                      linewidth = 0.3, linetype = 'dotted'),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = 'black', linewidth = 0.2)
  ) +
  scale_x_continuous(breaks = seq(0, 10, by = 2)) +
  ggplot2::labs(y = str_wrap('収 入', width = 1), x = '期') +
  ggplot2::scale_color_grey(name = 'シナリオ',
                            labels = c('EDLP200万円', 'EDLP210万円',
                                       'EDLP220万円', 'EDLP230万円',
                                       'EDLP240万円', 'EDLP250万円'))

ggsave(here("03_Dynamic_Single_Agent_Ch06_07/output/Counter2_sim_revenue_EDLP.pdf"), device = cairo_pdf)

### 値下げをした場合の需要と収入の推移 ----

# 需要の推移を図示
sim_df %>% 
  dplyr::select(period, stringr::str_c('demand_', discount_scenerio_vec)) %>% 
  tidyr::pivot_longer(cols = - period) %>% 
  dplyr::rename(scenario = name, demand = value) %>% 
  ggplot2::ggplot(aes(x = period, y = demand, color = scenario)) +
  ggplot2::geom_line() +
  ggplot2::geom_point() + 
  ggplot2::theme_minimal() +
  ggplot2::theme(
    axis.title.y = element_text(angle = 0, vjust = 0.5),
    text = element_text(family = 'HiraKakuPro-W3'),
    panel.grid.major.y = element_line(color = 'grey',
                                      linewidth = 0.3, linetype = 'dotted'),
    panel.grid.minor.y = element_line(color = 'grey',
                                      linewidth = 0.3, linetype = 'dotted'),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = 'black', linewidth = 0.2)
  ) +
  scale_x_continuous(breaks = seq(0, 10, by = 2)) +
  ggplot2::labs(y = str_wrap('需 要 量', width = 1), x = '期') +
  ggplot2::scale_color_grey(name = 'シナリオ', 
                            labels = c('ベースライン', '永続的', '一時的'))

ggsave(here("03_Dynamic_Single_Agent_Ch06_07/output/Counter2_sim_demand_down.pdf"), device = cairo_pdf)


# 需要の累積和を図示
sim_df %>% 
  dplyr::filter(period <= 10) %>% 
  dplyr::select(period, stringr::str_c('demand_', discount_scenerio_vec)) %>% 
  dplyr::mutate(across(starts_with('demand'), cumsum)) %>% 
  tidyr::pivot_longer(cols = - period) %>% 
  dplyr::rename(scenario = name, demand = value) %>% 
  ggplot2::ggplot(aes(x = period, y = demand, color = scenario)) +
  ggplot2::geom_line() +
  ggplot2::geom_point() + 
  ggplot2::theme_minimal() +
  ggplot2::theme(
    axis.title.y = element_text(angle = 0, vjust = 0.5),
    text = element_text(family = 'HiraKakuPro-W3'),
    panel.grid.major.y = element_line(color = 'grey',
                                      linewidth = 0.3, linetype = 'dotted'),
    panel.grid.minor.y = element_line(color = 'grey',
                                      linewidth = 0.3, linetype = 'dotted'),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = 'black', linewidth = 0.2)
  ) +
  scale_x_continuous(breaks = seq(0, 10, by = 2)) +
  ggplot2::labs(y = str_wrap('累 積 需 要 量', width = 1), x = '期') +
  ggplot2::scale_color_grey(name = 'シナリオ', 
                            labels = c('ベースライン', '永続的', '一時的'))

ggsave(here("03_Dynamic_Single_Agent_Ch06_07/output/Counter2_sim_cumdemand_down.pdf"), device = cairo_pdf)


# ベースラインを100とした時の、需要の累積和を図示
sim_df %>% 
  dplyr::filter(period <= 10) %>% 
  dplyr::select(period, stringr::str_c('demand_', discount_scenerio_vec)) %>% 
  dplyr::mutate(across(starts_with('demand'), cumsum)) %>% 
  dplyr::mutate(Permanent = demand_Permanent / demand_Baseline * 100,
                Temporary = demand_Temporary / demand_Baseline * 100) %>%
  dplyr::select(period, Permanent, Temporary) %>%
  tidyr::pivot_longer(cols = - period) %>% 
  dplyr::rename(scenario = name, demand = value) %>% 
  ggplot2::ggplot(aes(x = period, y = demand, color = scenario)) +
  ggplot2::geom_line() +
  ggplot2::geom_point() + 
  ggplot2::theme_minimal() + 
  ggplot2::theme(
    axis.title.y = element_text(angle = 0, vjust = 0.5),
    text = element_text(family = 'HiraKakuPro-W3'),
    panel.grid.major.y = element_line(color = 'grey',
                                      linewidth = 0.3, linetype = 'dotted'),
    panel.grid.minor.y = element_line(color = 'grey',
                                      linewidth = 0.3, linetype = 'dotted'),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = 'black', linewidth = 0.2)
  ) +
  scale_x_continuous(breaks = seq(0, 10, by = 2)) +
  ggplot2::labs(y = str_wrap('累 積 需 要 量', width = 1), x = '期') +
  ggplot2::scale_color_grey(name = 'シナリオ', 
                            labels = c('永続的', '一時的'))

ggsave(here("03_Dynamic_Single_Agent_Ch06_07/output/Counter2_sim_cumdemand_down_100.pdf"), device = cairo_pdf)


# 収入(万円)の推移を図示
sim_df %>% 
  dplyr::select(period, stringr::str_c('revenue_', discount_scenerio_vec)) %>% 
  tidyr::pivot_longer(cols = - period) %>% 
  dplyr::rename(scenario = name, revenue = value) %>% 
  ggplot2::ggplot(aes(x = period, y = revenue * 100, color = scenario)) +
  ggplot2::geom_line() +
  ggplot2::geom_point() +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    axis.title.y = element_text(angle = 0, vjust = 0.5),
    text = element_text(family = 'HiraKakuPro-W3'),
    panel.grid.major.y = element_line(color = 'grey',
                                      linewidth = 0.3, linetype = 'dotted'),
    panel.grid.minor.y = element_line(color = 'grey',
                                      linewidth = 0.3, linetype = 'dotted'),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = 'black', linewidth = 0.2)
  ) +
  scale_x_continuous(breaks = seq(0, 10, by = 2)) +
  ggplot2::labs(y=str_wrap('収 入', width = 1), x = '期') +
  ggplot2::scale_color_grey(name = 'シナリオ', 
                            labels = c('ベースライン', '永続的', '一時的'))

ggsave(here("03_Dynamic_Single_Agent_Ch06_07/output/Counter2_sim_revenue_down.pdf"), device = cairo_pdf)


### 合計の需要と収入 ----

demand_df <- sim_df %>% 
  dplyr::summarize(across(starts_with('demand'), sum)) %>% 
  dplyr::rename_with(~str_replace(., 'demand_', '')) %>% 
  tidyr::pivot_longer(cols = everything()) %>% 
  dplyr::rename(scenario = name, demand = value)

revenue_df <- sim_df %>% 
  # 割引現在価値を計算するため、割引因子の列を追加
  dplyr::mutate(beta = beta^(period - 1)) %>% 
  dplyr::summarize(across(starts_with('revenue'), ~sum(.x * beta))) %>% 
  dplyr::rename_with(~str_replace(., 'revenue_', '')) %>% 
  tidyr::pivot_longer(cols = everything()) %>% 
  dplyr::rename(scenario = name, revenue = value)

demand_df %>% 
  dplyr::left_join(revenue_df, by = 'scenario') %>% 
  dplyr::mutate(revenue = revenue / 100) %>%  # 単位は億円
  dplyr::mutate(per_rev = revenue / demand * 10000) %>% { # 単位は万円 
    txt <- capture.output(knitr::kable(., format = "simple"))
    writeLines(txt, here("03_Dynamic_Single_Agent_Ch06_07/output/tab7_2_compare_scenario.txt"), useBytes = TRUE)
  }    

