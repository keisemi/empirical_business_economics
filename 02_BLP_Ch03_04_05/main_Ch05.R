
# パッケージの読み込み
library(tidyverse)
library(sjmisc)
library(knitr)
library(tictoc) #時間計測に用いるパッケージ
library(here)

source(here("02_BLP_CH03_04_05/function_Ch05.R"))

# 時間計測開始
tictoc::tic()

# データの準備 ----

data <- readr::read_csv(file = here("02_BLP_CH03_04_05/intermediate/data_for_estimation.csv"))

# ダミー変数の作成
# 外国車ダミー
data <- data %>%
  dplyr::mutate(Foreign_d = if_else(Type=="Foreign", 1, 0))

# FuelType=レギュラー のダミー
data <- data %>%
  dplyr::mutate(FuelRegular_d = if_else(FuelType=="レギュラー", 1, 0))

# 年ダミー
data <- data %>%
  dplyr::mutate(year = as.factor(year)) %>%
  sjmisc::to_dummy(year, suffix = "label") %>%
  dplyr::bind_cols(data) %>%
  dplyr::mutate(year = as.numeric(year))

# year 2006を基準に
data <- data %>%
  dplyr::select(-year_2006)

# capacity ダミー(4以下と5以上で分ける)
data <- data %>%
  dplyr::mutate(capacity_d = if_else(capacity > 4, 1, 0))

# ダミー変数を最後の列に移動
data <- data %>%
  dplyr::relocate(ends_with("_d"), starts_with("year_"), .after = last_col())

# ランダム係数ロジットモデルの推定 ----
# コード自体はCh04とほぼ同様

## データの下準備 ----

# まずデータをソートする。マーケット順、メーカー順かつ価格順
data <- data %>% 
  dplyr::arrange(year, Maker, price)

# データセットより抽出した行列等をdatalistに格納する
datalist <- list()

# 平均効用に入ってくる部分。内生性がある価格を含んでいる。
datalist$X1 <- data %>%
  dplyr::mutate(cons = 1) %>%
  # 価格を必ず2つ目にする
  dplyr::select(
    cons, price, FuelEfficiency, hppw, size, 
    capacity_d, FuelRegular_d, Foreign_d,
    starts_with("year_")
  ) %>% 
  as.matrix()

# ランダム係数とInteractする部分。
datalist$X2 <- data %>%
  dplyr::mutate(cons = 1) %>%
  # 価格を必ず1つ目にする
  dplyr::select(price) %>% 
  as.matrix()

# 操作変数行列。ここは外生変数と追加的な操作変数を含む
datalist$Z <- data %>%
  dplyr::mutate(cons = 1) %>%
  dplyr::select(
    cons, FuelEfficiency, hppw, size,
    capacity_d, FuelRegular_d, Foreign_d,
    starts_with("year_"),
    starts_with("iv_GH"),
    -ends_with("nest")
  ) %>% 
  as.matrix()

# 市場シェア
datalist$ShareVec <- data %>%
  dplyr::select(share) %>%
  as.matrix()

datalist$marketindex <- data$year
datalist$logitshare <- data$logit_share

# 観察数、年数を取得
N <- length(datalist$marketindex)
T <- length(unique(datalist$marketindex))                        

# 乱数を設定

set.seed(111)
Nsim <- 1000
datalist$draw_vec <- rnorm(Nsim * ncol(datalist$X2))

parameter <- list()
parameter$Nsim <- Nsim
parameter$T <- T 
parameter$N <- N
marketindex <- datalist$marketindex
uniquemarketindex <- sort(unique(marketindex))
temp1 <- matrix(rep(uniquemarketindex, N), T, N) %>% t()
temp2 <- matrix(rep(marketindex, T), N, T)
mkt_denom_d <- (temp1 == temp2) * 1
datalist$mkt_denom_d <- mkt_denom_d


## パラメタ推定及び標準誤差の計算 ----
 
# GMMの荷重行列のオプション
datalist$weight_mat_option <- "2SLS"

# グローバル変数としてアップデートする平均効用
delta_global <- data$logit_share

# 最適化
GMM_nonlinear <- stats::optim(
  par = c(0.7), 
  fn = GMM_obj,
  method = "L-BFGS-B",
  lower = c(0), 
  parameter = parameter,
  datalist = datalist,
  option = 0
)

# nonlinear parameter estimates
theta2_hat <- GMM_nonlinear$par

# 標準誤差の計算

GMM_linear <- GMM_obj(theta2_hat, parameter, datalist, option = 1)

# deltaと観測されないXiをdataに保存
delta <- GMM_linear$delta
Xi <- GMM_linear$Xi 

data <- data %>%
  dplyr::mutate(Xi = Xi, delta = delta)

# linear parameter estimates
theta1_hatmat <- GMM_linear$theta1
theta1_hat <- as.vector(theta1_hatmat)

se <- calculate_standard_error(theta2 = theta2_hat, parameter = parameter, datalist = datalist)


# 推定結果

theta2_hatmat <- as.matrix(theta2_hat)
rownames(theta2_hatmat)[1] <- "random_price"
beta_hat <- rbind(theta1_hatmat, theta2_hatmat) 
colnames(beta_hat) <- "coeff"
names(se)[length(GMM_linear$theta1) + 1] <- "random_price"

est_table <- as.data.frame(cbind(beta_hat, se))

est_table %>%
  knitr::kable(
    align = c("l", rep("c", 2)),
    digits = 2, 
    booktabs = FALSE, 
    format = "html",
    caption = "Estimates of The Parameters"
  ) %>%
  kableExtra::kable_styling(full_width = FALSE)

est_table %>%  {
  txt <- capture.output(knitr::kable(., format = "simple"))
  writeLines(txt, con = here("02_BLP_CH03_04_05/output/est_table.txt"), useBytes = TRUE)
}

# 限界費用の推定 ----

# 価格弾力性行列の計算。限界費用推定において用いる。
elasticity <- calculate_elasticity(
  datalist = datalist, 
  parameter = parameter,
  theta1 = theta1_hat, 
  theta2 = theta2_hat, 
  delta = delta
)

# 作成した関数を実行する
marketindex <- datalist$marketindex
Marginal_Cost_Est <- estimate_merginal_cost(data, marketindex, elasticity)

# 2016年の3つの企業のそれぞれの車種の限界費用を取り出す
# まず適当に3つの企業の名前を選ぶ
Maker <- unique(data$Maker)
Choice_Maker <- Maker[c(8, 10, 11)]

#2016年の限界費用を抽出
Marginal_Cost_2016_Est <- Marginal_Cost_Est$Marginal_Cost_2016

# Nevo (2000) Table 4
Table_4dt <-  Marginal_Cost_2016_Est %>%
  dplyr::select(Maker, Name, price, mc, margin) %>%
  # この後の表記に合わせて企業名を、Honda=Nippyo, Nissan=Brand A, Subaru=Brand B, Toyota=Brand Cとする
  dplyr::mutate(
    Maker = recode(Maker, Honda = "Nippyo", Nissan = "Brand_A", Subaru = "Brand_B", Toyota = "Brand_C")
  )

colnames(Table_4dt) <- c("Maker", "Name", "Price", "Marginal Cost", "Margin (p-mc)/p")

Table_4dt %>% 
  knitr::kable(
    align = rep("c", 5),
    digits = 3,
    booktabs = FALSE,
    format = "html",
    caption = "Predicted Marginal Costs"
  ) %>%
  kableExtra::kable_styling(full_width = FALSE)

Table_4dt %>%  {
  txt <- capture.output(knitr::kable(., format = "simple"))
  writeLines(txt, con = here("02_BLP_CH03_04_05/output/Table_4dt.txt"), useBytes = TRUE)
}

# distribution of margin
g <- Marginal_Cost_2016_Est %>%
  ggplot2::ggplot(aes(x = margin)) + 
  ggplot2::geom_histogram(binwidth = 0.01) +
  ggplot2::theme_bw()

plot(g)

ggsave(g, file = here("02_BLP_CH03_04_05/output/chap5_margin.pdf"), width = 9, height = 6, device = cairo_pdf)

# 合併シミュレーション ----

## 準備 ----

# 最新の2016年のみ使用
data_2016 <- data %>%
  dplyr::filter(year == 2016) %>%
  dplyr::left_join(Marginal_Cost_2016_Est, by = c("NameID", "Maker", "Name", "price")) %>%
  # 企業名を変更する
  dplyr::mutate(
    Maker = recode(Maker, 
                   Honda = "Nippyo", Nissan = "Brand_A", Subaru = "Brand_B", Toyota = "Brand_C")
  ) %>%
  # 合併後の企業
  dplyr::mutate(
    MakerNippyoA = recode(Maker, Nippyo = "Nippyo_A", Brand_A = "Nippyo_A"), 
    MakerNippyoB = recode(Maker, Nippyo = "Nippyo_B", Brand_B = "Nippyo_B")
  )

### 所有構造行列を作成する ----

# ownership matrix
J <- sum(marketindex == 2016)
Ownership_NippyoA <- matrix(0, J, J)
Ownership_NippyoB <- matrix(0, J, J)
Ownership_true <- matrix(0, J, J)

for (j in 1:J) {
  for (k in 1:J) {
    if (data_2016$MakerNippyoA[j] == data_2016$MakerNippyoA[k]) {
      Ownership_NippyoA[j, k] <- 1
    }
    if (data_2016$MakerNippyoB[j] == data_2016$MakerNippyoB[k]) {
      Ownership_NippyoB[j, k] <- 1
    }
    if (data_2016$Maker[j] == data_2016$Maker[k]) {
      Ownership_true[j, k] <- 1
    }
  }
}

### BLP推定同様の下準備 ----

# 2016年に絞った状態で推定に必要な行列・リストを作成
datalist_2016 <- list()

# マーケットとモデルの情報をGET
N <- length(data_2016$year)
T <- length(unique(data_2016$year))

# 平均効用に入ってくる部分。内生性がある価格を含んでいる。
datalist_2016$X1 <- as.data.frame(datalist$X1) %>%
  dplyr::mutate(year = datalist$marketindex) %>%
  dplyr::filter(year == 2016) %>%
  dplyr::select(-year) %>%
  as.matrix()

# ランダム係数とInteractする部分。
datalist_2016$X2 <- as.data.frame(datalist$X2) %>%
  dplyr::mutate(year = datalist$marketindex) %>%
  dplyr::filter(year == 2016) %>%
  dplyr::select(-year) %>%
  as.matrix()

# 操作変数行列。ここは外生変数と追加的な操作変数を含む
datalist_2016$Z <-  as.data.frame(datalist$Z) %>%
  dplyr::mutate(year = datalist$marketindex) %>%
  dplyr::filter(year == 2016) %>%
  dplyr::select(-year) %>%
  as.matrix()

# 市場シェア
datalist_2016$ShareVec <- data_2016 %>%
  dplyr::select(share) %>%
  as.matrix()

datalist_2016$marketindex <- data_2016$year
datalist_2016$logitshare <- data_2016$logit_share

# 乱数を用意
set.seed(111)

Nsim <- 1000

datalist_2016$draw_vec <- rnorm(Nsim * ncol(datalist_2016$X2))

parameter <- list()
parameter$Nsim <- Nsim

marketindex <- datalist_2016$marketindex
uniquemarketindex <- sort(unique(marketindex))
temp1 <- matrix(rep(uniquemarketindex, N), T, N) %>% t()
temp2 <- matrix(rep(marketindex, T), N, T)
mkt_denom_d <- (temp1 == temp2) * 1
datalist_2016$mkt_denom_d <- mkt_denom_d

### 均衡計算 ----

# 合併後の価格を計算

# 推定した費用を用いる
mc <- as.matrix(data_2016$mc)
# 残差として得られた観測されないXiは固定されているものとする
Xi <- as.matrix(data_2016$Xi)

## 日評自動車とA社が合併したケースのシミュレーション ----

# Nippyo-Brand A Simulation
# set the initial price
p_ini <- data_2016$price

p_NippyoA <- solve_equilibrium_price(
  datalist_2016, 
  p_ini, 
  Ownership_NippyoA, 
  parameter, 
  theta1_hat, 
  theta2_hat, 
  mc, 
  Xi
)

data_2016 <- data_2016 %>%
  dplyr::mutate(p_NippyoA = p_NippyoA[, 1])

## 日評自動車とB社が合併したケースのシミュレーション ----

# Nippyo-Brand B Simulation
# set the initial price
p_ini <- data_2016$price

p_NippyoB <- solve_equilibrium_price(
  datalist_2016, 
  p_ini, 
  Ownership_NippyoB, 
  parameter, 
  theta1_hat, 
  theta2_hat, 
  mc, 
  Xi
)

data_2016 <- data_2016 %>%
  dplyr::mutate(p_NippyoB = p_NippyoB[, 1])

# 合併時のシェアを計算
share_NippyoA <- calculate_mktshare_sim(
  datalist_2016, 
  data_2016$p_NippyoA, 
  parameter, 
  theta1_hat, 
  theta2_hat, 
  Xi
)

share_NippyoB <- calculate_mktshare_sim(
  datalist_2016, 
  data_2016$p_NippyoB, 
  parameter, 
  theta1_hat, 
  theta2_hat, 
  Xi
)

data_2016 <- data_2016 %>%
  dplyr::mutate(
    share_NippyoA = share_NippyoA,
    share_NippyoB = share_NippyoB
  )

## シミュレーション結果のまとめ ----

## 合併シミュレーションによる価格・販売台数変化

Table_5dt <- data_2016 %>%
  dplyr::mutate(
    p_NippyoA_PPC = (p_NippyoA - price) / price * 100,
    share_NippyoA_PPC = (share_NippyoA - share) / share * 100,
    p_NippyoB_PPC = (p_NippyoB - price) / price * 100,
    share_NippyoB_PPC = (share_NippyoB - share) / share * 100
  ) %>%
  dplyr::select(
    Maker, Name, 
    p_NippyoA_PPC, share_NippyoA_PPC, p_NippyoB_PPC, share_NippyoB_PPC
  ) %>%
  dplyr::filter(Maker %in% c("Nippyo","Brand_A","Brand_B","Brand_C"))

# 名前を付けなおす
colnames(Table_5dt) <- c("Maker", "Name", "p", "q", "p", "q")

# 表を作る
Table_5dt %>%
  knitr::kable(
    align = rep("c", 6),
    digits = 3,
    booktabs = FALSE,
    format = "html", 
    caption = "Predicted Percent Change in Prices and Quantities as a Result of Mergers"
  ) %>%
  kableExtra::add_header_above(c(" ", " ", "Nippyo and Brand A" = 2, "Nippyo and Brand B" = 2)) %>%
  kableExtra::kable_styling(full_width = FALSE)


Table_5dt %>%  {
  txt <- capture.output(knitr::kable(., format = "simple"))
  writeLines(txt, con = here("02_BLP_CH03_04_05/output/Table_5dt.txt"), useBytes = TRUE)
}


## 日評自動車とブランドAの合併後に価格が変化しないような限界費用を計算する

# 合併前の弾力性行列
elasmat_2016 <- list(elasticity$elasmat_2016)

# 日評自動車とブランドA社の名前を同じにする（所有構造を変更）
data_2016_NippyoA <- data_2016 %>%
  dplyr::mutate(Maker = MakerNippyoA) %>%
  dplyr::select(year, Maker, Name, NameID, price, share)

# 限界費用を計算（シェアと価格、弾力性は合併以前のものを使用）
mc_NippyoA_pfix <- estimate_merginal_cost(data_2016_NippyoA, marketindex, elasmat_2016)

mc_NippyoA_pfix_2016 <- mc_NippyoA_pfix$Marginal_Cost_2016

data_2016 <- data_2016 %>%
  dplyr::mutate(mc_NippyoA_pfix = mc_NippyoA_pfix_2016$mc)

## 日評自動車とブランドBの合併後に価格が変化しないような限界費用を計算する

# 日評自動車とブランドB社の名前を同じにする（所有構造を変更）
data_2016_NippyoB <- data_2016 %>%
  dplyr::mutate(Maker = MakerNippyoB) %>%
  dplyr::select(year, Maker, Name, NameID, price, share)

# 限界費用を計算（シェアと価格、弾力性は合併以前のものを使用）
mc_NippyoB_pfix <- estimate_merginal_cost(data_2016_NippyoB, marketindex, elasmat_2016)

mc_NippyoB_pfix_2016 <- mc_NippyoB_pfix$Marginal_Cost_2016

data_2016 <- data_2016 %>%
  dplyr::mutate(mc_NippyoB_pfix = mc_NippyoB_pfix_2016$mc)

# 合併による価格上昇を打ち消すのに必要な限界費用変化を求める

Table_6dt <- data_2016 %>%
  dplyr::mutate(
    mc_NippyoA_PPC = (mc_NippyoA_pfix - mc) / mc * 100,
    mc_NippyoB_PPC = (mc_NippyoB_pfix - mc) / mc * 100
  ) %>%
  dplyr::select(Maker, Name, mc_NippyoA_PPC, mc_NippyoB_PPC) %>%
  dplyr::filter(Maker %in% c("Nippyo","Brand_A","Brand_B","Brand_C"))

# 名前を付けなおす
colnames(Table_6dt) <- c("Origin Maker", "Name", "Nippyo and Brand A", "Nippyo and Brand B")

# 表を作る
Table_6dt %>%
  knitr::kable(
    align = rep("c", 4),
    digits = 3,
    booktabs = FALSE,
    format = "html", 
    caption = "Percent Reduction in Marginal Costs Required for No Change in Predicted Postmerger Prices"
  ) %>%
  kableExtra::kable_styling(full_width = FALSE)

Table_6dt %>%  {
  txt <- capture.output(knitr::kable(., format = "simple"))
  writeLines(txt, con = here("02_BLP_CH03_04_05/output/Table_6dt.txt"), useBytes = TRUE)
}

## 合併シミュレーションの厚生分析 ----

HH_2016 <- unique(data_2016$HH)

CS_2016 <- calculate_consumer_surplus(
  datalist_2016, 
  data_2016$price, 
  parameter, 
  theta1_hat, 
  theta2_hat, 
  Xi, 
  HH_2016
)

CS_NippyoA <- calculate_consumer_surplus(
  datalist_2016, 
  data_2016$p_NippyoA, 
  parameter, 
  theta1_hat, 
  theta2_hat, 
  Xi, 
  HH_2016
)

CS_NippyoB <- calculate_consumer_surplus(
  datalist_2016, 
  data_2016$p_NippyoB, 
  parameter, 
  theta1_hat, 
  theta2_hat, 
  Xi, 
  HH_2016
)

CV_NippyoA <- CS_NippyoA - CS_2016
CV_NippyoB <- CS_NippyoB - CS_2016

# 総余剰の変化を導

pro_rev_2016 <- calculate_profit(
  data_2016$Maker, 
  data_2016$price, 
  data_2016$mc, 
  data_2016$share, 
  HH_2016
)  

pro_rev_NippyoA <- calculate_profit(
  data_2016$Maker, 
  data_2016$p_NippyoA, 
  data_2016$mc, 
  data_2016$share_NippyoA, 
  HH_2016
) 

pro_rev_NippyoB <- calculate_profit(
  data_2016$Maker, 
  data_2016$p_NippyoB, 
  data_2016$mc, 
  data_2016$share_NippyoB, 
  HH_2016
) 

# Total Surplus Change
TS_Chage_NippyoA <- CV_NippyoA + sum(pro_rev_NippyoA$profit - pro_rev_2016$profit)
TS_Chage_NippyoB <- CV_NippyoB + sum(pro_rev_NippyoB$profit - pro_rev_2016$profit)

# 消費者余剰と総余剰の変化の表
Table_7dt1 <- dplyr::tibble(
  clo1 = c("Consumer surplus", "Total Welfare"),
  col2 = c(CV_NippyoA, TS_Chage_NippyoA),
  col3 = c(CV_NippyoB, TS_Chage_NippyoB)
)

colnames(Table_7dt1) <- c("Measure", "Nippyo and Brand A", "Nippyo and Brand B")

Table_7dt1 %>%
  knitr::kable(
    align = rep("c", 3), 
    digits = 1,
    booktabs = FALSE, 
    format = "html",
    caption = "Consumer Surplus and Welfare"
  ) %>%
  kableExtra::kable_styling(full_width = FALSE)

Table_7dt1 %>%  {
  txt <- capture.output(knitr::kable(., format = "simple"))
  writeLines(txt, con = here("02_BLP_CH03_04_05/output/Table_7dt1.txt"), useBytes = TRUE)
}


Table_7dt2 <- generate_profit_revenue_table(pro_rev_2016, pro_rev_NippyoA, pro_rev_NippyoB)

# 表を作る
Table_7dt2 %>%
  knitr::kable(
    align = rep("c", 5),
    digits = 1,
    booktabs = FALSE, 
    format = "html",
    caption = "Change in Profits and Revenues"
  ) %>%
  kableExtra::add_header_above(c(" ", "Nippyo and Brand A" = 2, "Nippyo and Brand B" = 2)) %>%
  kableExtra::kable_styling(full_width = FALSE)

Table_7dt2 %>%  {
  txt <- capture.output(knitr::kable(., format = "simple"))
  writeLines(txt, con = here("02_BLP_CH03_04_05/output/Table_7dt2.txt"), useBytes = TRUE)
}

# シミュレーションに関する計算時間
# 時間計測終了
tictoc::toc()

# Appendix 1: 追加的なシミュレーション ----
# 「総余剰が変化しない」ために要求される合併企業の限界費用削減について計算

# 日評自動車とブランドAについて限界費用削減を計算 ----

# 時間計測開始
tictoc::tic()

# Nippyo-Brand A 限界費用削減
# cost_redについて総余剰変化は単調減少するため、二分法を適用

# 両企業が限界費用削減 （日評自動車とbrand Aを取り出す）
cost_red_firm <- c("Nippyo", "Brand_A")

distance <- 100
lambda <- 10^(-6)
max_cost_red <- 1
min_cost_red <- 0

# ここのループにおいて、総余剰を合併前の水準に保つために必要な限界費用削減率を計算している。
while (distance > lambda) {
  mid_cost_red <- (max_cost_red + min_cost_red) / 2
  
  mid_eval <- calculate_surplus_change(
    cost_red = mid_cost_red, 
    cost_red_firm = cost_red_firm, 
    Ownership = Ownership_NippyoA,
    data = data_2016, 
    datalist = datalist_2016,  
    parameter = parameter, 
    theta1 = theta1_hat,
    theta2 = theta2_hat,
    HH = HH_2016, 
    p_pre = data_2016$p_NippyoA,
    pro_rev_pre = pro_rev_2016,
    CS_pre = CS_2016
  )
  
  if (mid_eval > 0) {
    min_cost_red <- mid_cost_red
  } else {
    max_cost_red <- mid_cost_red
  }
  
  distance <- abs(mid_eval - 0)
  print(distance)
}

cost_red_NippyoA <- mid_cost_red

data_2016 <- data_2016 %>%
  dplyr::mutate(mc_NippyoA_TSfix = if_else(Maker %in% cost_red_firm, mc * cost_red_NippyoA, mc))

p_NippyoA_TSfix <- solve_equilibrium_price(
  datalist_2016,
  data_2016$p_NippyoA,
  Ownership_NippyoA,
  parameter, 
  theta1_hat, 
  theta2_hat, 
  data_2016$mc_NippyoA_TSfix, 
  Xi
)

share_NippyoA_TSfix <- calculate_mktshare_sim(
  datalist_2016, 
  p_NippyoA_TSfix, 
  parameter, 
  theta1_hat, 
  theta2_hat, 
  Xi
)

data_2016 <- data_2016 %>%
  dplyr::mutate(p_NippyoA_TSfix = p_NippyoA_TSfix, 
                share_NippyoA_TSfix = share_NippyoA_TSfix)

## 日評自動車とブランドBについて限界費用削減を計算 ----

# Nippyo-Brand B 限界費用削減
# cost_redについて総余剰変化は単調減少するため、二分法を適用

# 両企業が限界費用削減 （日評自動車とbrand Bを取り出す）
cost_red_firm <- c("Nippyo","Brand_B")

distance <- 100
lambda <- 10^(-6)
max_cost_red <- 1
min_cost_red <- 0

# ここのループにおいて、総余剰を合併前の水準に保つために必要な限界費用削減率を計算している。
while (distance > lambda){
  mid_cost_red <- (max_cost_red + min_cost_red) / 2
  mid_eval <- calculate_surplus_change(
    cost_red = mid_cost_red, 
    cost_red_firm = cost_red_firm,
    Ownership = Ownership_NippyoB,
    data = data_2016, 
    datalist = datalist_2016,  
    parameter = parameter,
    theta1 = theta1_hat,
    theta2 = theta2_hat,
    HH = HH_2016, 
    p_pre = data_2016$p_NippyoB,
    pro_rev_pre = pro_rev_2016,
    CS_pre = CS_2016
  )
  
  if (mid_eval > 0) {
    min_cost_red <- mid_cost_red
  } else {
    max_cost_red <- mid_cost_red
  }
  
  distance <- abs(mid_eval - 0)
  print(distance)
}

cost_red_NippyoB <- mid_cost_red

data_2016 <- data_2016 %>%
  dplyr::mutate(mc_NippyoB_TSfix = if_else(Maker %in% cost_red_firm, mc * cost_red_NippyoB, mc))

p_NippyoB_TSfix <- solve_equilibrium_price(
  datalist_2016, 
  data_2016$p_NippyoB, 
  Ownership_NippyoB,
  parameter, 
  theta1_hat, 
  theta2_hat, 
  data_2016$mc_NippyoB_TSfix, 
  Xi
)

share_NippyoB_TSfix <- calculate_mktshare_sim(
  datalist_2016, 
  p_NippyoB_TSfix, 
  parameter, 
  theta1_hat, 
  theta2_hat, 
  Xi
)

data_2016 <- data_2016 %>%
  dplyr::mutate(p_NippyoB_TSfix = p_NippyoB_TSfix, 
                share_NippyoB_TSfix = share_NippyoB_TSfix)

print(cost_red_NippyoB)

Table_7dt3 <- dplyr::tibble(
  col1 = c("Cost reduction"),
  col2 = (1 - cost_red_NippyoA) * 100,
  col3 = (1 - cost_red_NippyoB) * 100
)

colnames(Table_7dt3) <- c("Measure", "Nippyo and Brand A", "Nippyo and Brand B")

Table_7dt3 %>%
  knitr::kable(
    align = rep("c", 3),
    digits = 2,
    booktabs = FALSE, 
    format = "html",
    caption = "Cost Reduction (Percent)"
  ) %>% 
  kableExtra::kable_styling(full_width = FALSE)

Table_7dt3 %>% {
  txt <- capture.output(knitr::kable(., format = "simple"))
  writeLines(txt, con = here("02_BLP_CH03_04_05/output/Table_7dt3.txt"), useBytes = TRUE)
}

## 限界費用削減をした場合の利潤と収入の変化 ----

# 日評自動車とブランドA社が合併した場合の利益の導出
pro_rev_NippyoA_rc <- calculate_profit(
  data_2016$Maker, 
  data_2016$p_NippyoA_TSfix, 
  data_2016$mc_NippyoA_TSfix,
  data_2016$share_NippyoA_TSfix, 
  HH_2016
) 

# 日評自動車とブランドB社が合併した場合の利益の導出
pro_rev_NippyoB_rc <- calculate_profit(
  data_2016$Maker, 
  data_2016$p_NippyoB_TSfix, 
  data_2016$mc_NippyoB_TSfix,
  data_2016$share_NippyoB_TSfix, 
  HH_2016
) 

Table_7dt4 <- generate_profit_revenue_table(pro_rev_2016, pro_rev_NippyoA_rc, pro_rev_NippyoB_rc)

# 表を作る
Table_7dt4 %>%
  knitr::kable(
    align = rep("c", 5),
    digits = 2,
    booktabs = FALSE, 
    format = "html",
    caption = "Change in Profits and Revenues"
  ) %>% 
  kableExtra::add_header_above(c(" ", "Nippyo and Brand A" = 2, "Nippyo and Brand B" = 2)) %>%
  kableExtra::kable_styling(full_width = FALSE)

Table_7dt4 %>% {
  txt <- capture.output(knitr::kable(., format = "simple"))
  writeLines(txt, con = here("02_BLP_CH03_04_05/output/Table_7dt4.txt"), useBytes = TRUE)
}

## 限界費用削減をした場合の全体の利潤と収入の変化率
print((sum(pro_rev_NippyoA_rc$profit) - sum(pro_rev_2016$profit)) / sum(pro_rev_2016$profit) * 100)
print((sum(pro_rev_NippyoA_rc$revenue) - sum(pro_rev_2016$revenue)) / sum(pro_rev_2016$revenue) * 100)
print((sum(pro_rev_NippyoB_rc$profit) - sum(pro_rev_2016$profit)) / sum(pro_rev_2016$profit) * 100)
print((sum(pro_rev_NippyoB_rc$revenue) - sum(pro_rev_2016$revenue)) / sum(pro_rev_2016$revenue) * 100)

# 時間計測終了
tictoc::toc()


