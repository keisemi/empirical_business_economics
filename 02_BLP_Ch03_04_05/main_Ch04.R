# 変数消去
rm(list = ls())

# パッケージの読み込み
library(tidyverse)
library(fixest)
library(summarytools)
library(sjmisc)
library(here)
library(showtext)
library(tictoc)

source(here("02_BLP_CH03_04_05/function_Ch03_04.R"))

# 日本語のPDF出力
showtext_auto()

# main_Ch03.Rで作成したクリーニング済みのデータを読み込む
data <- readr::read_csv(file = here("02_BLP_CH03_04_05/intermediate/data_cleaned.csv"))

# Logitモデルにおける価格弾力性
## イントロダクションの日評自動車

# 日評自動車のデータセットを作成する
data_NIPPYO <- data %>% 
  dplyr::filter(Nippyo == 1) %>% 
  dplyr::select(year, share, NameID, Sales, price, hppw, FuelEfficiency, size, Name) %>% 
  dplyr::mutate(
    log_sales = log(Sales),
    log_price = log(price)
  )

# 価格弾力性行列を作成する

data <- data %>%
  dplyr::mutate(logit_share = log(share) - log(share0))

# Differentiation IVを用いた結果
iv_GH <- fixest::feols(
  logit_share ~ hppw + FuelEfficiency + size | 0 |
    price ~ iv_GH_own_hppw + iv_GH_own_FuelEfficiency + iv_GH_own_size +
    iv_GH_other_hppw + iv_GH_other_FuelEfficiency + iv_GH_other_size,
  data = data
)

dt2016 <- data_NIPPYO %>% 
  dplyr::filter(year == 2016) %>% 
  dplyr::select(price, share, NameID, Name) %>% 
  dplyr::arrange(NameID)

price <- dt2016$price
share <- dt2016$share
NameID <- dt2016$NameID

own_elas <- iv_GH$coefficients["fit_price"] * price * (1 - share)

cross_elas <- (-1) * iv_GH$coefficients["fit_price"] * price * share
J <- length(own_elas)

elas_mat <- matrix(rep(cross_elas, J), nrow = J, ncol = J)
diag(elas_mat) <- own_elas


betard <- dt2016 %>%
  dplyr::filter(Name == "アルファード") %>%
  dplyr::pull(NameID)

sedan <- dt2016 %>%
  dplyr::filter(Name == "カローラ") %>%
  dplyr::pull(NameID)

suv <- dt2016 %>%
  dplyr::filter(Name == "ジューク") %>%
  dplyr::pull(NameID)

kei <- dt2016 %>%
  dplyr::filter(Name == "タント") %>%
  dplyr::pull(NameID)

elas_mat_restricted <- elas_mat[NameID %in% c(betard, sedan, suv, kei), NameID %in% c(betard, sedan, suv, kei)]
Name <- c("ベータード", "セダン(A)", "SUV(B)", "軽自動車(C)")
colnames(elas_mat_restricted) <- Name
rownames(elas_mat_restricted) <- Name

elas_mat_restricted %>% {
  txt <- capture.output(knitr::kable(., format = "simple"))
  writeLines(txt, con = here("02_BLP_CH03_04_05/output/tab4_1_elas_mat_restricted.txt"), useBytes = TRUE)
}

# 入れ子型ロジット推定のためのBLP操作変数の作成

# まず、マーケット・企業レベルにおける、各製品属性の和と自乗和を計算する
# ここでacross関数は、最初に文字列ベクトルで指定した変数について、後ろにリスト内で定義した操作を適用している
data <- data %>%
  dplyr::group_by(year, Maker, Type) %>%
  dplyr::mutate(
    dplyr::across(c("hppw", "FuelEfficiency", "size"),
                  list(sum_own = ~ sum(.x, na.rm = TRUE))),
    dplyr::across(c("hppw", "FuelEfficiency", "size"),
                  list(sqr_sum_own = ~ sum(.x^2, na.rm = TRUE))),
    group_n = n()
  ) %>%
  dplyr::ungroup()

# 次に、マーケットレベルでの、各製品属性の和を計算する
data <- data %>% 
  dplyr::group_by(year, Type) %>%
  dplyr::mutate( 
    dplyr::across(c("hppw", "FuelEfficiency", "size"),
                  list(sum_mkt = ~ sum(.x, na.rm = TRUE))),
    dplyr::across(c("hppw", "FuelEfficiency", "size"),
                  list(sqr_sum_mkt = ~ sum(.x^2, na.rm = TRUE))),
    mkt_n = n()
  ) %>%
  dplyr::ungroup() 

# 以上で定義した変数を利用して、まずBLP操作変数を構築する。
data <- data %>% 
  dplyr::mutate(
    iv_BLP_own_hppw_nest = hppw_sum_own - hppw,
    iv_BLP_own_FuelEfficiency_nest = FuelEfficiency_sum_own - FuelEfficiency,
    iv_BLP_own_size_nest = size_sum_own - size,
    iv_BLP_other_hppw_nest = hppw_sum_mkt - hppw_sum_own,
    iv_BLP_other_FuelEfficiency_nest = FuelEfficiency_sum_mkt - FuelEfficiency_sum_own,
    iv_BLP_other_size_nest = size_sum_mkt - size_sum_own, 
    iv_BLP_own_num_nest = group_n - 1, 
    iv_BLP_other_num_nest = mkt_n - group_n
  ) 

# 続いて、Differentiation IVを構築する。
data <- data %>% 
  dplyr::mutate(
    iv_GH_own_hppw_nest = (group_n - 1) * hppw^2 + 
      (hppw_sqr_sum_own - hppw^2) - 2 * hppw * (hppw_sum_own - hppw),
    iv_GH_own_FuelEfficiency_nest = (group_n - 1) * FuelEfficiency^2 + 
      (FuelEfficiency_sqr_sum_own - FuelEfficiency^2) - 
      2 * FuelEfficiency * (FuelEfficiency_sum_own - FuelEfficiency),
    iv_GH_own_size_nest = (group_n - 1) * size^2 + 
      (size_sqr_sum_own - size^2) - 2 * size * (size_sum_own - size),
    iv_GH_other_hppw_nest = (mkt_n - group_n) * hppw^2 + 
      (hppw_sqr_sum_mkt - hppw_sqr_sum_own) - 
      2 * hppw * (hppw_sum_mkt - hppw_sum_own),
    iv_GH_other_FuelEfficiency_nest = (mkt_n - group_n) * FuelEfficiency^2 + 
      (FuelEfficiency_sqr_sum_mkt - FuelEfficiency_sqr_sum_own) - 
      2 * FuelEfficiency * (FuelEfficiency_sum_mkt - FuelEfficiency_sum_own),
    iv_GH_other_size_nest = (mkt_n - group_n) * size^2 + 
      (size_sqr_sum_mkt - size_sqr_sum_own) - 2 * size * (size_sum_mkt - size_sum_own)
  ) %>%
  dplyr::select(
    -starts_with("sum_own"),
    -starts_with("sum_mkt"),
    -starts_with("sqr_sum_own"),
    -starts_with("sqr_sum_mkt"),
    -mkt_n,
    -group_n
  )

# 入れ子型ロジットモデルの推定とその応用 ----
# 推定に際しては、OLS、BLP操作変数による結果、そしてDifferentiation IVによる結果の3通りについて比較する。

# まず被説明変数を定義する。
data <- data %>%
  dplyr::mutate(logit_share = log(share) - log(share0))

# インサイドシェアを定義する。
data <- data %>% 
  dplyr::group_by(year, Type) %>% 
  dplyr::mutate(sum_year_body = sum(Sales)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(
    inside_share = Sales / sum_year_body, 
    log_inside_share = log(Sales / sum_year_body)
  )

# データの保存(第５章で使用)
readr::write_csv(data, here("02_BLP_CH03_04_05/intermediate/data_for_estimation.csv"))

# OLSの結果
ols <- fixest::feols(
  logit_share ~ price + log_inside_share + hppw + FuelEfficiency + size | 0,
  data = data
)

# GH操作変数を用いた結果
iv_GH <-
  fixest::feols(
    logit_share ~ hppw + FuelEfficiency + size | 0 |
      price + log_inside_share ~ iv_GH_own_hppw_nest + iv_GH_own_FuelEfficiency_nest +
      iv_GH_own_size_nest + iv_GH_other_hppw_nest + iv_GH_other_FuelEfficiency_nest +
      iv_GH_other_size_nest,
    data = data
  )

# BLP操作変数を用いた結果
iv_BLP_nest <-
  fixest::feols(
    logit_share ~ hppw + FuelEfficiency + size | 0 |
      price + log_inside_share ~ iv_BLP_own_hppw_nest + 
      iv_BLP_own_FuelEfficiency_nest + iv_BLP_own_size_nest +
      iv_BLP_other_hppw_nest + iv_BLP_other_FuelEfficiency_nest + 
      iv_BLP_other_size_nest + iv_BLP_own_num_nest +  iv_BLP_other_num_nest,
    data = data
  )

# 推定結果をレポートする。
tbl_nest_iv <- fixest::etable(
  list("OLS" = ols, "IV_BLP" = iv_BLP_nest),  
  se = "hetero",
  fitstat = c("r2", "n", "ivf"), 
  signif.code = NA, 
  dict = c(
    price = "自動車価格",
    hppw = "馬力／重量",
    FuelEfficiency = "燃費 (キロメートル／ 1 リットル)",
    size = "サイズ",
    `(Intercept)` = "定数項"
  ),
  digits = 2,
  digits.stats = 2,
  depvar = FALSE
  ) %>% {
  txt <- capture.output(knitr::kable(., format = "simple"))
  writeLines(txt, con = here("02_BLP_CH03_04_05/output/tab_4_2_nest_iv.txt"), useBytes = TRUE)
}

## 入れ子型ロジットモデルにおける自己価格弾力性の計算

alpha1 <- ols$coefficients["price"]
sigma1 <- ols$coefficients["log_inside_share"]

alpha2 <- iv_BLP_nest$coefficients["fit_price"]
sigma2 <- iv_BLP_nest$coefficients["fit_log_inside_share"]

data_elas <- data %>% 
  dplyr::mutate(
    own_elas_ols = alpha1 * price * (1 - sigma1 * inside_share - (1 - sigma1) * share) / (1 - sigma1),
    own_elas_ivblp = alpha2 * price * (1 - sigma2 * inside_share - (1 - sigma2) * share) / (1 - sigma2)
  )

tbl_nest_own_elas <- data_elas %>% 
  dplyr::select(starts_with("own_elas")) %>% 
  summarytools::descr(
    transpose = TRUE, 
    stats = c("mean", "sd", "med", "min", "max"),
    order = "preserve"
  ) %>% {
    txt <- capture.output(knitr::kable(., format = "simple"))
    writeLines(txt, con = here("02_BLP_CH03_04_05/output/tbl_nest_own_elas.txt"), useBytes = TRUE)
  }

data_NIPPYO <- data %>% 
  dplyr::filter(Nippyo == 1) %>%
  dplyr::select(year, share, Type, inside_share, NameID, 
                Sales, price, hppw, FuelEfficiency, size) %>% 
  dplyr::mutate(
    log_sales = log(Sales),
    log_price = log(price)
  ) 

dt2016 <- data_NIPPYO %>% 
  dplyr::filter(year == 2016) %>% 
  dplyr::select(price, Type, share, inside_share, NameID) %>% 
  dplyr::arrange(NameID)

price <- dt2016$price
share <- dt2016$share
NameID <- dt2016$NameID
inside_share <- dt2016$inside_share
group <- dt2016$Type


# 自己価格弾力性
own_elas <- alpha2 * price * (1 - sigma2 * inside_share - (1 - sigma2) * share) / (1 - sigma2)

# グループの外との交差価格弾力性
cross_elas_othergroup <- (-1) * alpha2 * price * share
J <- length(own_elas)
cross_elas_othergroup <- matrix(rep(cross_elas_othergroup, J), nrow = J, ncol = J)

# グループの中での交差価格弾力性（行列表記）
price_l_mat <- matrix(rep(price, J), nrow = J, ncol = J)
share_l_mat <- matrix(rep(share, J), nrow = J, ncol = J)
insideshare_l_mat <- matrix(rep(inside_share, J), nrow = J, ncol = J)
cross_elas_samegroup <- (-1) * alpha2 * price_l_mat * (sigma2 * insideshare_l_mat + (1 - sigma2) * share_l_mat) / (1 - sigma2)

# Indicator of group
temp_mat1 <- matrix(rep(group, J), nrow = J, ncol = J)
temp_mat2 <- t(temp_mat1)
ind_same_group <- (temp_mat1 == temp_mat2)
ind_other_group <- (temp_mat1 != temp_mat2)

# Construct Elasticity Matrix
elas_mat_nl <- cross_elas_samegroup * ind_same_group + cross_elas_othergroup * ind_other_group
diag(elas_mat_nl) <- own_elas

elas_mat_nl_restricted <- elas_mat_nl[NameID %in% c(betard, sedan, suv, kei), NameID %in% c(betard, sedan, suv, kei)]

Name <- c("ベータード", "セダン(A)", "SUV(B)", "軽自動車(C)")

colnames(elas_mat_nl_restricted) <- Name
rownames(elas_mat_nl_restricted) <- Name

elas_mat_nl_restricted %>% {
  txt <- capture.output(knitr::kable(., format = "simple"))
  writeLines(txt, con = here("02_BLP_CH03_04_05/output/tab4_3_elas_mat_nl_restricted.txt"), useBytes = TRUE)
}

# ランダム係数ロジットモデルの推定 ----

# 推定のOverview
# Step 1はデータに関する下準備である。
# Step 2から4では、GMM推定のために必要な関数を順番に用意していく。
# 1. 市場シェアを計算する関数 (Step 2)
# 2. 縮小写像によるBerryインバージョンを行う関数 (Step 3)
# 3. GMMの目的関数 (Step 4)
# 市場シェアを計算する関数はStep 3・4両方で使用する。また、Berryインバージョンを行う関数はGMM目的関数を評価する際に必要となる。
# GMM目的関数を定義した後、数値最適化でGMM推定を行う。
# Step 5では推定したパラメタの標準誤差の計算を行う。最後に推定結果をテーブルでまとめる。



## Step 1: 下準備 ----

# まずデータをソートする。マーケット順かつモデル順
data <- data %>% 
  dplyr::arrange(year, NameID)

# マーケットとモデルの情報をGET
Info <- data %>% 
  dplyr::select(year, NameID)

N <- length(Info$year)
T <- length(unique(Info$year))

# 平均効用に入ってくる部分。内生性がある価格を含んでいる。
X1 <- data %>%
  dplyr::mutate(cons = 1) %>%
  dplyr::select(cons, price, FuelEfficiency, hppw, size) %>%
  as.matrix()

# ランダム係数とInteractする部分。
X2 <- data %>%
  dplyr::mutate(cons = 1) %>%
  dplyr::select(price, cons, size) %>%
  as.matrix()

# 操作変数行列。ここは外生変数と追加的な操作変数を含む
Z <- data %>%
  dplyr::mutate(cons = 1) %>%
  dplyr::select(
    cons, FuelEfficiency, hppw, size,
    starts_with("iv_GH"),
    -ends_with("nest")
  ) %>%
  as.matrix()

# 市場シェア
ShareVec <- data %>%
  dplyr::select(share) %>%
  as.matrix()

datalist <- list()
datalist$X1 <- X1
datalist$X2 <- X2
datalist$Z <- Z
datalist$ShareVec <- ShareVec
datalist$marketindex <- data$year
datalist$logitshare <- data$logit_share


set.seed(111)

Nsim <- 500

draw_vec <- rnorm(Nsim * ncol(X2))

datalist$draw_vec <- draw_vec

parameter <- list()
parameter$Nsim <- Nsim
parameter$T <- T 
parameter$N <- N
parameter$theta2 = c(0.001, 0.001)

marketindex <- datalist$marketindex
uniquemarketindex <- sort(unique(marketindex))
temp1 <- matrix(rep(uniquemarketindex, N), T, N) %>% t()
temp2 <- matrix(rep(marketindex, T), N, T)
tempmat <- (temp1 == temp2) * 1
datalist$tempmat <- tempmat

## Step 2-4: GMM推定 ----

# GMMの荷重行列のオプション
datalist$weight_mat_option <- "2SLS"

# Contraction mappingで用いる初期値
datalist$delta_ini <- data$logit_share

# 最適化
tic()
result <- stats::optim(
  par = c(0.2, 10, 0.1),
  fn = GMM_obj,
  method = "L-BFGS-B",
  lower = c(0, 0, 0),
  parameter = parameter,
  datalist = datalist,
  option = 0
)
toc()

## Step 5: 標準誤差の計算 ----

result2 <-GMM_obj(result$par, parameter, datalist, option = 1)

delta <- result2$delta

se <- calculate_standard_error(theta2 = result$par, parameter = parameter, datalist = datalist)

## ランダム係数ロジットの推定結果 ----
beta_hat <- rbind(result2$beta_hat, result$par[1], result$par[2], result$par[3])

rownames(beta_hat)[rownames(beta_hat) == "cons"] <- "定数項：平均"
rownames(beta_hat)[rownames(beta_hat) == "price"] <- "価格：平均"
rownames(beta_hat)[rownames(beta_hat) == "size"] <- "サイズ：平均"
rownames(beta_hat)[rownames(beta_hat) == "FuelEfficiency"] <- "燃費"
rownames(beta_hat)[rownames(beta_hat) == "hppw"] <- "馬力"
rownames(beta_hat)[length(result2$beta_hat) + 1] <- "価格：標準偏差"
rownames(beta_hat)[length(result2$beta_hat) + 2] <- "定数項：標準偏差"
rownames(beta_hat)[length(result2$beta_hat) + 3] <- "サイズ：標準偏差"

result_table <-
  as.data.frame(beta_hat) %>%
  dplyr::mutate(se = se) %>%
  dplyr::rename(
    推定値 = share,
    標準誤差 = se
  )

result_table <- result_table[c("定数項：平均",
                               "定数項：標準偏差",
                               "価格：平均",
                               "価格：標準偏差",
                               "サイズ：平均",
                               "サイズ：標準偏差",
                               "燃費",
                               "馬力"), ]

result_table %>% {
  txt <- capture.output(knitr::kable(., format = "simple"))
  writeLines(txt, con = here("02_BLP_CH03_04_05/output/tab4_4_rand_coef_logit_result.txt"), useBytes = TRUE)
}

## 価格弾力性行列の計算 ----

# 弾力性の表示の前処理
parameter$theta2 <- result$par

elasticity <- calculate_elasticity(
  datalist = datalist, 
  parameter = parameter, 
  beta_hat = beta_hat, 
  delta = delta
)

NameID2016 <- data %>% 
  dplyr::filter(year == 2016) %>%
  dplyr::select(NameID) %>% 
  dplyr::pull()


elas_mat_blp_restricted <- elasticity$elasmat_2016[NameID2016 %in% c(betard, sedan, suv, kei), NameID2016 %in% c(betard, sedan, suv, kei)]

Name <- c("ベータード", "セダン(A)", "SUV(B)", "軽自動車(C)")
colnames(elas_mat_blp_restricted) <- Name
rownames(elas_mat_blp_restricted) <- Name

elas_mat_blp_restricted  %>% {
  txt <- capture.output(knitr::kable(., format = "simple"))
  writeLines(txt, con = here("02_BLP_CH03_04_05/output/tab4_5_elas_mat_blp_restricted.txt"), useBytes = TRUE)
}

# 応用：プライシング ----
## 需要曲線と収入曲線を描く

NIPPYOautoIDvec <- data_NIPPYO %>%
  dplyr::distinct(NameID) %>%
  dplyr::pull()

NameID_target <- betard

pricevec <- seq(from = 1.73, to = 4, by = 0.01)
pivec <- numeric(length(pricevec))
pivec2 <- numeric(length(pricevec))

for (i in 1:length(pricevec)) {
  # ベータードのみの利潤
  pivec[i] <- calculate_revenue(
    pricevec[i], 
    data, 
    datalist, 
    result2, 
    parameter, 
    option = "ownpi"
  )
  
  # ベータードのみの利潤+日評自動車の他の車種の収入
  pivec2[i] <- calculate_revenue(
    pricevec[i], 
    data, 
    datalist, 
    result2, 
    parameter, 
    option = "totalpi"
  )
}


# ベータードのみ
dt_fig1 <- dplyr::tibble(
  price = pricevec * 100, 
  pi1 = pivec * 100 / 10000) %>% 
  tidyr::pivot_longer(-price, names_to = "name")

fig_pi1 <- ggplot2::ggplot(dt_fig1, aes(price, value)) + 
  ggplot2::geom_point() +
  ggplot2::theme_bw(base_family = "HiraKakuPro-W3") +
  ggplot2::theme(
    axis.title.y = element_text(angle = 0, vjust = 0.5),
    panel.grid.major.y = element_line(
      color = "grey",
      linewidth = 0.3,
      linetype = "dotted"
    ),
    panel.grid.minor.y = element_line(
      color = "grey",
      linewidth = 0.3,
      linetype = "dotted"
    ),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.2)) + 
  ggplot2::scale_x_continuous(breaks = seq(0, 500, 50)) + 
  ggplot2::labs(
    y = "収入(億円)",
    x = "価格(万円)"
  ) 

plot(fig_pi1)

ggsave(fig_pi1, file = here("02_BLP_CH03_04_05/output/fig4_2_revenue_betard.pdf"), width = 9, height = 6, device = cairo_pdf)

# ベータードのみの利潤+日評自動車の他の車種の収入
dt_fig2 <- dplyr::tibble(
  price = pricevec * 100, 
  pi2 = pivec2 * 100 / 10000
) %>% 
  tidyr::pivot_longer(-price, names_to = "name")

fig_pi2 <- ggplot2::ggplot(dt_fig2, aes(price, value)) + 
  ggplot2::geom_point() +
  ggplot2::theme_bw(base_family = "HiraKakuPro-W3") +
  ggplot2::theme(
    axis.title.y = element_text(angle = 0, vjust = 0.5),
    panel.grid.major.y = element_line(
      color = "grey",
      linewidth = 0.3,
      linetype = "dotted"
    ),
    panel.grid.minor.y = element_line(
      color = "grey",
      linewidth = 0.3,
      linetype = "dotted"
    ),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.2)) +
  ggplot2::scale_x_continuous(breaks = seq(0, 500, 50)) + 
  ggplot2::labs(
    y = "収入(億円)",
    x = "価格(万円)"
  ) 

plot(fig_pi2)

ggsave(fig_pi2, file = here("02_BLP_CH03_04_05/output/fig4_2_revenue_all.pdf"), width = 9, height = 6, device = cairo_pdf)

# まとめたテーブル
dt_pricing <- dplyr::tibble(
  price = pricevec * 100, 
  pi1 = pivec * 100 / 10000,
  pi2 = pivec2 * 100 / 10000) 

dt_pricing %>% 
  tidyr::pivot_longer(-price, names_to = "name") -> dt_forfig

write.csv(dt_forfig, file = here("02_BLP_Ch03_04_05/output/dt_prof_Betard.csv"))


# 数値最適化を使って、利潤を最大にする価格を求める

# ベータードのみ
optimprice3 <- stats::optimise(
  calculate_revenue,
  interval = c(0.3, 5),
  maximum = TRUE, 
  data = data,
  datalist = datalist,
  result2 = result2,
  parameter = parameter,
  option = "ownpi"
)

print(optimprice3)


# ベータードのみの利潤+日評自動車の他の車種の収入
optimprice4 <- stats::optimise(
  calculate_revenue,
  interval = c(0.3, 5),
  maximum = TRUE, 
  data = data,
  datalist = datalist, 
  result2 = result2,
  parameter = parameter, 
  option = "totalpi"
)

print(optimprice4)

# ベータードのみの利潤を最大にするような価格における、「ベータードのみの利潤+日評自動車の他の車種の収入」
calculate_revenue(
  price = optimprice3$maximum,
  data = data, 
  datalist = datalist, 
  result2 = result2,
  parameter = parameter, 
  option = "totalpi"
) -> total_rev

# 最適価格のまとめ

sink(here("02_BLP_Ch03_04_05/output/Ch04_opt_price.txt"))
print("ベータードの利潤のみ最適化")
print("ベータード最適価格と目的関数(ベータード利潤)")
print(optimprice3)
print("ベータードのみの利潤を最大にするような価格における、「ベータードのみの利潤+日評自動車の他の車種の収入」")
print(total_rev)

print("ベータード利潤+他の収入最大化")
print("ベータード最適価格と、目的関数(ベータード利潤+その他収入)")
print(optimprice4)
sink()

