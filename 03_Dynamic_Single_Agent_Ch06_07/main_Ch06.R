
# パッケージを読み込む
library(tidyverse)
library(dplyr)
library(knitr)
library(skimr)
library(numDeriv)
library(evd)
library(here)
library(ggplot2)

# plot3Dパッケージは内部で tcltk / X11 を呼び出している
# macOSではX11が標準で入っていないため、XQuartzをインストールする必要がある
library(plot3D)

source(here("03_Dynamic_Single_Agent_Ch06_07/function.R"))


# パラメータの設定 ----

# 走行距離にかかるtheta_cと車の価格にかかるtheta_p
theta_true <- c(theta_c = 15, theta_p = 6)

# 時間割引率
beta <- 0.95

# オイラー定数
Euler_const <- digamma(1)

# 消費者は車を購入するかどうか決めるため選択肢の数は2つ
num_choice <- 2


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
# したがって、状態変数の数は価格の状態変数の数と走行距離の状態変数の数の積となる
num_states <- num_price_states * num_mileage_states

# 価格、走行距離の状態変数の組み合わせ(p,m)を1つのデータフレームで表す
state_df <- 
  dplyr::tibble(
    # stateにidを番号付ける
    state_id = 1:num_states,
    # 順番は (p,m) = (2000,0), (2100,0), ...,(2500,0), (2000,5), (2100,5), ...
    price_id = rep(1:num_price_states,
                   times = num_mileage_states),
    mileage_id = rep(1:num_mileage_states,
                     each = num_price_states),
    price = rep(price_states,
                times = num_mileage_states),
    mileage = rep(mileage_states,
                  each = num_price_states)
  )

# 下3行を表示
state_df %>% tail(3)

# 走行距離の遷移行列のパラメタを設定し、遷移行列を作成する
kappa_true <- c(0.4, 0.1)

mileage_trans_mat_true <- gen_mileage_trans(kappa_true)

# 走行距離の遷移行列の4行4列までを表示
mileage_trans_mat_true[1:4, 1:4, 1]

# 価格の遷移行列のパラメタを設定し、遷移行列を作成する
lambda_true <- c(0.1, 0.2, 0.2, 0.2, 0.2,
                 0.1, 0.2, 0.2, 0.2, 0.2,
                 0.1, 0.1, 0.2, 0.2, 0.1,
                 0.1, 0.1, 0.2, 0.2, 0.1,
                 0.05, 0.05, 0.1, 0.1, 0.2,
                 0.05, 0.05, 0.1, 0.1, 0.2)

price_trans_mat_true <- gen_price_trans(lambda_true)

# 価格の遷移行列を表示
price_trans_mat_true

# コントロール変数毎の遷移行列を作成
trans_mat_true <- list()

# 車を購入しない場合の遷移行列
trans_mat_true$not_buy <- mileage_trans_mat_true[,,1] %x% price_trans_mat_true

# 車を購入する場合の遷移行列
trans_mat_true$buy <- mileage_trans_mat_true[,,2] %x% price_trans_mat_true


# 定常状態での価格の分布を計算
# 以下を満たすような price_dist_steady を求める
# price_dist_steady %*% price_trans_mat == price_dist_steady

# 固有値/固有ベクトルを求める
# 固有値が1となる固有ベクトルは1つだけ（1つめ）
price_trans_eigen <- eigen(t(price_trans_mat_true))

# 価格の定常分布を求める
price_dist_steady <- price_trans_eigen$vectors[, 1] / sum(price_trans_eigen$vectors[, 1])

price_dist_steady


# Vを求める
start_time <- proc.time()

V_true <- contraction(theta_true, beta, trans_mat_true, state_df)

end_time <- proc.time()

cat('Runtime:\n')
print((end_time - start_time)[[3]])

# 選択毎の価値関数を定義する
U_true <- flow_utility(theta_true, state_df)

V_CS_true <- U_true + beta * cbind(trans_mat_true$not_buy %*% V_true,
                                   trans_mat_true$buy %*% V_true)

colnames(V_CS_true) <- c('V_not_buy', 'V_buy')

# state(p,m)ごとに、logitで計算される理論上の条件付き購入確率を計算
prob_buy_true_mat <- matrix(exp(V_CS_true[, 'V_buy']) / rowSums(exp(V_CS_true)), 
                            nrow = num_price_states,
                            ncol = num_mileage_states)
prob_buy_true_mat

# logitで計算される理論上の条件付き購入確率を図示
# 3次元プロット
cairo_pdf(here('03_Dynamic_Single_Agent_Ch06_07/output/CCP_true_3D.pdf'))

par(family = 'HiraKakuProN-W3')

hist3D(x = mileage_states * 100,
       y = price_states * 100,
       z = t(prob_buy_true_mat),
       zlim = c(0, 0.8),
       bty = 'g',
       phi = 10,
       theta = -60,
       axes = TRUE,
       label = TRUE,
       xlab = '走行距離（万km）',
       ylab = '価格（万円）',
       zlab = '購入確率',
       col = 'grey',
       border = 'black',
       shade = 0.4,
       ticktype = 'detailed',
       space = 0.05,
       d = 2,
       cex.axis = 0.8)

dev.off()

# logitで計算される理論上の条件付き購入確率を図示
# 2次元プロット:価格ごとに図示

data.frame(price = state_df$price,
           mileage = state_df$mileage,
           prob_buy = as.vector(prob_buy_true_mat)) %>%
  ggplot2::ggplot(aes(x = mileage * 100,
                      y = prob_buy,
                      color = factor(price))) +
  ggplot2::geom_line() +
  ggplot2::geom_point() +
  ggplot2::theme_minimal() +
  ggplot2::scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  ggplot2::theme(
    axis.title.y = element_text(angle = 0, vjust = 0.5),
    text = element_text(family = 'HiraKakuPro-W3'),
    panel.grid.major.y = element_line(color = 'grey',
                                      linewidth = 0.3,
                                      linetype = 'dotted'),
    panel.grid.minor.y = element_line(color = 'grey',
                                      linewidth = 0.3,
                                      linetype = 'dotted'),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = 'black', linewidth = 0.2)
  ) +
  scale_x_continuous(breaks = seq(0, 10, by = 2)) +
  ggplot2::labs(y = str_wrap('購 買 確 率', width = 1),
                x = '走行距離（万km）') +
  ggplot2::scale_color_grey(name = '価格', 
                            labels = c('200万円', '210万円',
                                       '220万円', '230万円',
                                       '240万円', '250万円')) 

ggsave(here('03_Dynamic_Single_Agent_Ch06_07/output/fig6_3_CCP_true_2D.pdf'), device = cairo_pdf)


# シミュレーション ----

# サンプルサイズを決める

# 1000人の消費者が存在
num_consumer <- 1000

# 50年分の年次データを生成した後、最後の10年のみが観察できるとする
num_period <- 50
num_period_obs <- 10

# 総観察数
num_obs <- num_consumer * num_period

# 累積分布確率を持つように遷移行列を変換（行方向に足し上げる）
trans_mat_cum <- list()
trans_mat_cum$not_buy <- t(apply(trans_mat_true$not_buy, 1, cumsum))
trans_mat_cum$buy <- t(apply(trans_mat_true$buy, 1, cumsum))

# 乱数を固定
set.seed(1)

# 生成するデータの元となるdata.frameを作成
data_gen <- 
  dplyr::tibble(
    consumer = rep(1:num_consumer, each = num_period),
    period = rep(1:num_period, times = num_consumer),
    eps_type1_not_buy = evd::rgev(num_obs),
    eps_type1_buy = evd::rgev(num_obs),
    eps_unif = runif(num_obs),
    eps_price_state_unif = runif(num_obs),
    state_id = 0,
    action = 0
  )


data_gen <- data_gen %>%
  # 消費者ごとにデータを分割
  dplyr::group_split(consumer) %>%
  # 上記の関数で定義したデータの生成過程をすべての消費者に対して行う
  purrr::map_dfr(generate_data,
                 V_CS = V_CS_true,
                 state_df = state_df, 
                 price_dist_steady = price_dist_steady) %>% 
  # 最後の10年のみが観察できるとする
  dplyr::filter(period > (num_period - num_period_obs)) %>% 
  # 状態変数を表した列を追加
  dplyr::left_join(state_df, by = 'state_id')

data_gen %>% tail(3)

# 擬似データの出力・保存
data_gen %>% write.csv(here("03_Dynamic_Single_Agent_Ch06_07/intermediate/Chap6_data.csv")) 

# 不要なオブジェクトの削除
rm(V_CS_true, trans_mat_cum)


# 生成したデータの要約統計
data_gen %>% 
  dplyr::select(price, mileage, action) %>%
  skimr::skim() %>% 
  skimr::yank('numeric') %>% 
  dplyr::select(skim_variable, mean, sd, p0, p100) %>% { 
    txt <- capture.output(knitr::kable(., format = "simple"))
    writeLines(txt, here("03_Dynamic_Single_Agent_Ch06_07/output/tab6_2_descr.txt"), useBytes = TRUE)
    }

# 走行距離の分布
data_gen %>%
  ggplot2::ggplot(aes(x = mileage * 100)) + 
  ggplot2::geom_histogram(bins = 21) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.title.y = element_text(angle = 0, vjust = 0.5),
    text = element_text(family = 'HiraKakuPro-W3'),
    legend.position = 'none',
    panel.grid.major.y = element_line(color = 'grey', 
                                      linewidth = 0.3,
                                      linetype = 'dotted'),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = 'black',
                             linewidth = 0.2)
  ) +
  scale_x_continuous(breaks = seq(0, 10, 2)) + 
  ggplot2::labs(y = str_wrap('頻 度', width = 1),
                x = '走行距離（万km）')

ggsave(here('03_Dynamic_Single_Agent_Ch06_07/output/dist_mile.pdf'), device = cairo_pdf)

# 価格の分布
data_gen %>%
  ggplot2::ggplot(aes(x = price * 100)) + 
  ggplot2::geom_histogram(bins = 6) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.title.y = element_text(angle = 0, vjust = 0.5),
    text = element_text(family = 'HiraKakuPro-W3'),
    legend.position = 'none',
    panel.grid.major.y = element_line(color = 'grey',
                                      linewidth = 0.3,
                                      linetype = 'dotted'),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = 'black',
                             linewidth = 0.2)
  ) +
  scale_x_continuous(breaks = seq(200, 250, by = 10)) +
  ggplot2::labs(y = str_wrap('頻 度', width = 1), x = '価格（万円）')

ggsave(here('03_Dynamic_Single_Agent_Ch06_07/output/dist_price.pdf'), device = cairo_pdf)



# それぞれの走行距離において購入した割合を観察
data_gen %>% 
  dplyr::group_by(mileage) %>% 
  dplyr::summarize(num_state = n(),
                   sum_action = sum(action),
                   .groups = 'drop') %>% 
  dplyr::mutate(prob_buy = sum_action / num_state) %>% 
  ggplot2::ggplot(aes(x = mileage * 100, y = prob_buy)) + 
  ggplot2::geom_bar(stat = 'identity') +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.title.y = element_text(angle = 0, vjust = 0.5),
    text = element_text(family = 'HiraKakuPro-W3'),
    legend.position = 'none',
    panel.grid.major.y = element_line(color = 'grey',
                                      linewidth = 0.3,
                                      linetype = 'dotted'),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = 'black', linewidth = 0.2)
  ) +
  scale_x_continuous(breaks = seq(0, 10, 2)) + 
  ggplot2::labs(y = str_wrap('購 入 確 率', width = 1),
                x = '走行距離（万km）')

ggsave(here('03_Dynamic_Single_Agent_Ch06_07/output/CCP_mile.pdf'), device = cairo_pdf)

# それぞれの価格において購入した割合を観察
data_gen %>% 
  dplyr::group_by(price) %>% 
  dplyr::summarize(num_state = n(),
                   sum_action = sum(action),
                   .groups = 'drop') %>% 
  dplyr::mutate(prob_buy = sum_action / num_state) %>% 
  ggplot2::ggplot(aes(x = price * 100, y = prob_buy)) + 
  ggplot2::geom_bar(stat = 'identity') +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.title.y = element_text(angle = 0, vjust = 0.5),
    text = element_text(family = 'HiraKakuPro-W3'),
    legend.position = 'none',
    panel.grid.major.y = element_line(color = 'grey',
                                      linewidth = 0.3,
                                      linetype = 'dotted'),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = 'black', linewidth = 0.2)
  ) +
  scale_x_continuous(breaks = seq(200, 250, by = 10)) +
  ggplot2::labs(y=str_wrap('購 入 確 率', width = 1),
                x = '価格（万円）')

ggsave(here('03_Dynamic_Single_Agent_Ch06_07/output/CCP_price.pdf'), device = cairo_pdf)


# state(p,m)ごとに、観測された条件付き購入確率を計算
prob_buy_obs_mat <- 
  data_gen %>%
  dplyr::group_by(mileage,price) %>%
  dplyr::summarize(num_state = n(),
                   sum_action = sum(action),
                   .groups = 'drop') %>%
  dplyr::mutate(prob_buy = sum_action / num_state) %>% 
  dplyr::select(prob_buy) %>% 
  as.matrix() %>% 
  matrix(nrow = num_price_states, ncol = num_mileage_states)

prob_buy_obs_mat


cairo_pdf(here('03_Dynamic_Single_Agent_Ch06_07/output/CCP_3D.pdf'))

par(family = 'HiraKakuProN-W3')

hist3D(x = mileage_states * 100,
       y = price_states * 100,
       z = t(prob_buy_obs_mat),
       zlim = c(0, 0.8),
       bty = 'g',
       phi = 10,
       theta = -60,
       axes = TRUE,
       label = TRUE,
       xlab = '走行距離（万km）',
       ylab = '価格（万円）',
       zlab = '購入確率',
       col = 'grey',
       border = 'black',
       shade = 0.4,
       ticktype = 'detailed',
       space = 0.05,
       d = 2,
       cex.axis = 0.8)

dev.off()



# 遷移行列の推定 ----

# 走行距離の遷移行列の推定

# 遷移行列の推定で使うため、ラグ変数を追加
data_gen <- data_gen %>% 
  dplyr::group_by(consumer) %>% 
  dplyr::mutate(lag_price_id = lag(price_id),
                lag_mileage_id = lag(mileage_id),
                lag_action = lag(action)) %>% 
  dplyr::ungroup() 

# それぞれの確率が実現した観察の数を数える
num_cond_obs_mileage <- data_gen %>% 
  # 1期目は推定に使えないため落とす
  dplyr::filter(period != (num_period - num_period_obs + 1)) %>% 
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

kappa_est[1] <- (num_cond_obs_mileage[2] *
                   (num_cond_obs_mileage[2] + num_cond_obs_mileage[3] + num_cond_obs_mileage[4])) / 
  ((num_cond_obs_mileage[2] + num_cond_obs_mileage[3]) * 
     (num_cond_obs_mileage[1] + num_cond_obs_mileage[2] + num_cond_obs_mileage[3] + num_cond_obs_mileage[4]))

kappa_est[2] <- (num_cond_obs_mileage[3] * 
                   (num_cond_obs_mileage[2] + num_cond_obs_mileage[3] + num_cond_obs_mileage[4])) /
  ((num_cond_obs_mileage[2] + num_cond_obs_mileage[3]) * 
     (num_cond_obs_mileage[1] + num_cond_obs_mileage[2] + num_cond_obs_mileage[3] + num_cond_obs_mileage[4]))


# 最尤法の解析解により標準誤差を求める
Infomat_mileage_est <- matrix(0, nrow = 2, ncol = 2)

# 最尤法のフィッシャー情報量を求める
Infomat_mileage_est[1,1] <- 
  (num_cond_obs_mileage[1] / (1 - kappa_est[1] - kappa_est[2])^2) +
  (num_cond_obs_mileage[2] / kappa_est[1]^2) +
  (num_cond_obs_mileage[4] / (kappa_est[1]+kappa_est[2])^2)

Infomat_mileage_est[1,2] <- 
  (num_cond_obs_mileage[1] / (1 - kappa_est[1] - kappa_est[2])^2) +
  (num_cond_obs_mileage[4] / (kappa_est[1]+kappa_est[2])^2)

Infomat_mileage_est[2,1] <- Infomat_mileage_est[1,2]

Infomat_mileage_est[2,2] <- 
  (num_cond_obs_mileage[1] / (1 - kappa_est[1] - kappa_est[2])^2) +
  (num_cond_obs_mileage[3] / kappa_est[2]^2) +
  (num_cond_obs_mileage[4] / (kappa_est[1]+kappa_est[2])^2)

# 逆行列の対角要素の平方根が標準誤差になる
kappa_se <- sqrt(diag(solve(Infomat_mileage_est)))

dplyr::tibble(kappa_est, kappa_se) %>% {
  txt <- capture.output(knitr::kable(., format = "simple"))
  writeLines(txt, here("03_Dynamic_Single_Agent_Ch06_07/output/tbl_kappa_est.txt"), useBytes = TRUE)
  }


## 価格の遷移行列の推定

# それぞれの確率が実現した観察の数を数える
num_cond_obs_price <- data_gen %>% 
  # 1期目は推定に使えないため落とす
  dplyr::filter(period != (num_period - num_period_obs + 1)) %>% 
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

# 最尤法の解析解により標準誤差を求める
lambda_se <- c()

for (i in 1:num_price_states) {
  # 最尤法のフィッシャー情報量を求める
  Infomat_price_est <- 
    diag(num_cond_obs_price[i,],
         num_price_states)[-i, -i] / 
    (lambda_est_mat[-i, -i]^2) + 
    (num_cond_obs_price[i, i] / 
       lambda_est_mat[i, i]^2) *
    matrix(1, num_price_states, num_price_states)[-i, -i]
  lambda_se <- c(
    lambda_se,
    # 逆行列の対角要素の平方根が標準誤差になる
    sqrt(diag(solve(Infomat_price_est)))
  )
}


lambda_se_mat <- 
  c(0, lambda_se[1], lambda_se[2], lambda_se[3], lambda_se[4], lambda_se[5],
    lambda_se[6], 0, lambda_se[7], lambda_se[8], lambda_se[9], lambda_se[10],
    lambda_se[11], lambda_se[12], 0, lambda_se[13], lambda_se[14], lambda_se[15],
    lambda_se[16], lambda_se[17], lambda_se[18], 0, lambda_se[19], lambda_se[20],
    lambda_se[21], lambda_se[22], lambda_se[23], lambda_se[24], 0, lambda_se[25],
    lambda_se[26], lambda_se[27], lambda_se[28], lambda_se[29], lambda_se[30], 0) %>% 
  matrix(ncol = num_price_states, nrow = num_price_states, byrow = T)

lambda_se_mat

lambda_est <- as.vector(t(lambda_est_mat))[c(-1, -8, -15, -22, -29, -36)]

dplyr::tibble(lambda_est, lambda_se) %>% { 
  txt <- capture.output(knitr::kable(., format = "simple"))
  writeLines(txt, here("03_Dynamic_Single_Agent_Ch06_07/output/tbl_lambda_est.txt"), useBytes = TRUE)
  }


# パラメータの推定 ----

## 静学的なロジットによる推定 ----

#ここでは消費者が動学的にではなく静学的に意思決定を行っていると仮定し、単純なロジットモデルを用いてtheta_cとtheta_pを推定する。

# パラメータtheta_cとtheta_pを推定する
start_time <- proc.time()

# 最適化
logit_stat_opt <- optim(theta_true,
                        logLH_stat,
                        state_df = state_df,
                        df = data_gen, 
                        control = list(fnscale = -1), 
                        method = 'Nelder-Mead')

end_time <- proc.time()

cat('Runtime:\n')
print((end_time - start_time)[[3]])

theta_est_stat <- logit_stat_opt$par
theta_est_stat

# 標準誤差を推定する
hessian_stat <- numDeriv::hessian(func = logLH_stat,
                                  x = theta_est_stat, 
                                  state_df = state_df,
                                  df = data_gen)

theta_se_stat <- sqrt(diag(solve(-hessian_stat)))

dplyr::tibble(theta_est_stat, theta_se_stat) %>% { 
  txt <- capture.output(knitr::kable(., format = "simple"))
  writeLines(txt, here("03_Dynamic_Single_Agent_Ch06_07/output/tbl_theta_est_static.txt"), useBytes = TRUE)
  }

## 入れ子不動点アルゴリズムによる推定 ----

# 推定された遷移行列を取得
trans_mat_hat <- list()

trans_mat_hat$not_buy <- gen_mileage_trans(kappa_est)[,, 1] %x% gen_price_trans(lambda_est)

trans_mat_hat$buy <- gen_mileage_trans(kappa_est)[,, 2] %x% gen_price_trans(lambda_est)

start_time <- proc.time()

# 最適化
NFXP_opt <- optim(theta_true,
                  logLH_nfxp,
                  beta = beta,
                  trans_mat = trans_mat_hat,
                  state_df = state_df,
                  df = data_gen, 
                  control = list(fnscale = -1), 
                  method = 'Nelder-Mead')

end_time <- proc.time()

cat('Runtime:\n')
print((end_time - start_time)[[3]])

theta_est <- NFXP_opt$par
theta_est

# 標準誤差を推定する
hessian <- numDeriv::hessian(func = logLH_nfxp,
                             x = theta_est, 
                             beta = beta,
                             trans_mat = trans_mat_hat,
                             state_df = state_df,
                             df = data_gen)

theta_se <- sqrt(diag(solve(-hessian)))

dplyr::tibble(theta_est, theta_se) %>% { 
  txt <- capture.output(knitr::kable(., format = "simple"))
  writeLines(txt, here("03_Dynamic_Single_Agent_Ch06_07/output/tbl_theta_est_nfxp_ch6.txt"), useBytes = TRUE)
  }

# 推定結果の比較

dplyr::tibble(
  algorithm = c('Static', 'NFXP', 'True'),
  theta_c = c(theta_est_stat[1], theta_est[1], theta_true[1]),
  theta_se_c = c(theta_se_stat[1], theta_se[1], NA),
  theta_p = c(theta_est_stat[2], theta_est[2], theta_true[2]),
  theta_se_p = c(theta_se_stat[2], theta_se[2], NA)
  ) %>% { 
  txt <- capture.output(knitr::kable(., format = "simple"))
  writeLines(txt, here("03_Dynamic_Single_Agent_Ch06_07/output/tab6_3_compare_est.txt"), useBytes = TRUE)
  }



