# 変数消去
rm(list = ls())

# パッケージの読み込み
library(tidyverse)
library(fixest)
library(summarytools)
library(sjmisc)
library(here)
library(showtext)

source(here("02_BLP_CH03_04_05/function_Ch03_04.R"))

# 日本語のPDF出力
showtext_auto()

# データの読み込み ----

# 自動車データ
data <- readr::read_csv(file = here("02_BLP_CH03_04_05/data/CleanData_20180222_nippyo.csv"))

# 家計数（潜在的な市場規模）データ
dataHH <- readr::read_csv(file = here("02_BLP_CH03_04_05/data/HHsize.csv"))

# 消費者物価指数
dataCPI <- readr::read_csv(file = here("02_BLP_CH03_04_05/data/zni2015s.csv"), locale = locale(encoding = "shift-jis"))

dataCPI <- dataCPI[6:56, ]

dataCPI <- dataCPI %>% 
  dplyr::rename(
    year = '類・品目',
    cpi = '総合'
  ) %>% 
  dplyr::select(year, cpi) %>% 
  dplyr::mutate(
    year = as.numeric(year), 
    cpi = as.numeric(cpi)
  )

# データクリーニング ----

# 自動車データを整理しつつ、家計データをマージする

# 必要な変数のみをキープ
data <- data %>%
  dplyr::select(
    Maker, Type, Name, Year, Sales, Model,
    Nippyo, price, kata,
    weight, capacity, FuelType, FuelEfficiency, HorsePower,
    overall_length, overall_width, overall_height
  ) %>%
  dplyr::rename(year = Year)

# 家計サイズをマージする
data <- data %>%
  dplyr::left_join(dataHH, by = "year")

# CPIをマージする
data <- data %>% 
  dplyr::left_join(dataCPI, by = "year") 

# 燃費が欠損しているデータがある。今回は簡便な処理として観測から落とす。
data <- data %>%
  dplyr::filter(is.na(FuelEfficiency) == 0)

# 価格の実質化を行う。ここでは、2016年を基準年とする。
# また、価格の単位を100万円にする。元のデータは1万円
cpi2016 <- dataCPI %>% 
  dplyr::filter(year == 2016) %>% 
  dplyr::select(cpi) %>% 
  as.double()

data <- data %>% 
  dplyr::mutate(price = price / (cpi / cpi2016)) %>% 
  dplyr::mutate(price = price / 100) %>%
  dplyr::select(-cpi)

# サイズ(高さ＊幅＊長さ)、燃費の重量に対する比率を定義する
data <- data %>%
  dplyr::mutate(size = (overall_length / 1000) * (overall_width / 1000) * (overall_height / 1000)) %>%
  dplyr::mutate(hppw = HorsePower / weight) %>%
  dplyr::select(-HorsePower, -weight) %>% 
  dplyr::select(-starts_with("overall")) 

# 自動車の車種IDを作成する
data <- data %>%
  dplyr::group_by(Name) %>%
  dplyr::mutate(NameID = dplyr::cur_group_id()) %>%
  dplyr::ungroup() %>%
  dplyr::relocate(NameID, year)

# マーケットシェアとOutside option shareを定義する
# ここでは、マーケット(年)ごとに総販売台数(変数`inside_total`)を定義し、それを用いてアウトサイドオプションのシェアを計算している。
data <- data %>%
  dplyr::group_by(year) %>%
  dplyr::mutate(inside_total = sum(Sales)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(outside_total = HH - inside_total) %>%
  dplyr::mutate(share = Sales / HH) %>%
  dplyr::mutate(share0 = outside_total / HH) %>%
  dplyr::select(-inside_total, -outside_total)


# 操作変数の構築
# ここでは、BLP型操作変数と、Gandhi-HoudeのDifferentiation IVを定義する


#### 解説 
#コードが若干複雑なので、構築のアイデアを先に説明しよう。
# まず、Gandhi-HoudeのDifferentiation IVは以下のように定義される。
#
# $$
# \begin{split}Z_{jtk}^{\text{Quad,Other}}(X)=\sum_{k\in J_{ft}\setminus\{j\}}d_{jkt\ell}^{2},\\
# Z_{jtk}^{\text{Quad,Rival}}(X)=\sum_{k\notin J_{ft}}d_{jkt\ell}^{2}.
# \end{split}
# $$
#
# まず$Z_{jtk}^{\text{Quad,Other}}$の作成方法を考える．
# 具体的なイメージを持つべく、以下のように，市場と企業でグループ化したうえで，同じ企業が生産する財$1,2,3,4$と，財の特質$X$を持ったデータ(data)を考える．
#
# | 市場 | 企業 | 財  | 特性$X$ | 特性の差の二乗$d^{2}$                                            |
# |-------|-------|-------|-------|---------------------------------------------|
# | 1    | 1    | 1   | A       | $\left(A-B\right)^{2}+\left(A-C\right)^{2}+\left(A-D\right)^{2}$ |
# | 1    | 1    | 2   | B       | $\left(B-A\right)^{2}+\left(B-C\right)^{2}+\left(B-D\right)^{2}$ |
# | 1    | 1    | 3   | C       | $\left(C-A\right)^{2}+\left(C-B\right)^{2}+\left(C-D\right)^{2}$ |
# | 1    | 1    | 4   | D       | $\left(D-A\right)^{2}+\left(D-B\right)^{2}+\left(D-C\right)^{2}$ |
#
# ここで，財$1$についての$d^{2}$を考えると，以下のように書ける．
# $$
# \left(A-B\right)^{2}+\left(A-C\right)^{2}+\left(A-D\right)^{2}=3A^{2}+\left(B^{2}+C^{2}+D^{2}\right)-2A\left(B+C+D\right)
# $$
# これを一般化した形で書くと、
# $$
# d^{2}=\left(\mathrm{nrow(data)}-1\right)*X^{2}+\left(\mathrm{sum}(X^{2})-X^{2}\right)-2*X*\left(\mathrm{sum}(X)-X\right)
# $$
# となる。
# 以下では、このアイデアをコードに落とし込んだものを見ていこう。


# まず、マーケット・企業レベルにおける、各製品属性の和と自乗和を計算する。
# ここでacross関数は、最初に文字列ベクトルで指定した変数について、後ろにリスト内で定義した操作を適用している。
data <- data %>%
  dplyr::group_by(year, Maker) %>%
  dplyr::mutate(
    dplyr::across(c("hppw", "FuelEfficiency", "size"),
                  list(sum_own = ~ sum(.x, na.rm = TRUE))),
    dplyr::across(c("hppw", "FuelEfficiency", "size"),
                  list(sqr_sum_own = ~ sum(.x^2, na.rm = TRUE))),
    group_n = n()
  ) %>%
  dplyr::ungroup()

# 次に、マーケットレベルでの、各製品属性の和を計算する。
data <- data %>% 
  dplyr::group_by(year) %>%
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
    iv_BLP_own_hppw = hppw_sum_own - hppw,
    iv_BLP_own_FuelEfficiency = FuelEfficiency_sum_own - FuelEfficiency,
    iv_BLP_own_size = size_sum_own - size,
    iv_BLP_other_hppw = hppw_sum_mkt - hppw_sum_own,
    iv_BLP_other_FuelEfficiency = FuelEfficiency_sum_mkt - FuelEfficiency_sum_own,
    iv_BLP_other_size = size_sum_mkt - size_sum_own
  ) 

# 続いて、Differentiation IVを構築する。
data <- data %>% 
  dplyr::mutate(
    iv_GH_own_hppw = (group_n - 1) * hppw^2 + 
      (hppw_sqr_sum_own - hppw^2) - 2 * hppw * (hppw_sum_own - hppw),
    iv_GH_own_FuelEfficiency = (group_n - 1) * FuelEfficiency^2 + 
      (FuelEfficiency_sqr_sum_own - FuelEfficiency^2) - 
      2 * FuelEfficiency * (FuelEfficiency_sum_own - FuelEfficiency),
    iv_GH_own_size = (group_n - 1) * size^2 + 
      (size_sqr_sum_own - size^2) - 2 * size * (size_sum_own - size),
    iv_GH_other_hppw = (mkt_n - group_n) * hppw^2 + 
      (hppw_sqr_sum_mkt - hppw_sqr_sum_own) - 2 * hppw * (hppw_sum_mkt - hppw_sum_own),
    iv_GH_other_FuelEfficiency = (mkt_n - group_n) * FuelEfficiency^2 + 
      (FuelEfficiency_sqr_sum_mkt - FuelEfficiency_sqr_sum_own) - 
      2 * FuelEfficiency * (FuelEfficiency_sum_mkt - FuelEfficiency_sum_own),
    iv_GH_other_size = (mkt_n - group_n) * size^2 + 
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

# 保存
write.csv(data, here("02_BLP_CH03_04_05/intermediate/data_cleaned.csv"), row.names = FALSE)

# 記述統計と基礎的な分析 ----

# 日評自動車のデータセットを用意する。
data_NIPPYO <- data %>% 
  dplyr::filter(Nippyo == 1) %>%
  dplyr::select(
    Sales, price, hppw, FuelEfficiency, size
  ) %>% 
  dplyr::mutate(
    log_sales = log(Sales),
    log_price = log(price)
  ) 

# イントロダクションにおける回帰式をOLSで推定
ols_intro <- fixest::feols(
  log_sales ~ log_price + hppw + FuelEfficiency + size, data = data_NIPPYO
)

# 結果のアウトプット
tbl_ols_intro <- fixest::etable(
  "OLS_NIPPYO" = ols_intro, 
  se = "hetero",
  signif.code = NA, 
  fitstat = c("r2", "n" ), 
  dict = c(
    log_price = "log(価格)",
    hppw = "馬力／重量",
    FuelEfficiency = "燃費(キロメートル／ 1 リットル)",
    size = "サイズ",
    `(Intercept)` = "定数項"
  ),
  digits = 2,
  digits.stats = 2,
  depvar = FALSE) %>% {
    txt <- capture.output(knitr::kable(., format = "simple"))
    writeLines(txt, con = here("02_BLP_CH03_04_05/output/tbl_ols_intro.txt"), useBytes = TRUE)
    }

# 参考：価格と販売台数の散布図

# 日評自動車のデータのみ取り出す。
data_graph <- data_NIPPYO %>% dplyr::select(price, Sales)

g_1 <- data_graph %>% 
  ggplot2::ggplot(aes(price, Sales)) +
  ggplot2::scale_x_continuous(trans = "log10") +
  ggplot2::scale_y_continuous(trans = "log10") +
  ggplot2::geom_point() +
  ggplot2::geom_smooth(method = "lm", formula = y ~ x)

plot(g_1)


## 記述統計 ----

# データセット全体の記述統計を確認
tbl_data_summary <- data %>% 
  dplyr::select(
    Sales, price, FuelEfficiency, size, hppw
  ) %>% 
  summarytools::descr(
    transpose = TRUE,
    stats = c("mean", "sd", "min", "q1", "med", "q3", "max"),
    order = "preserve"
  ) %>% {
    txt <- capture.output(knitr::kable(., format = "simple"))
    writeLines(txt, con = here("02_BLP_CH03_04_05/output/tab3_1_data_summary.txt"), useBytes = TRUE)
  }

# ロジットモデルの推定とその応用 ----

# 推定に際しては、OLS、BLP操作変数による結果、そしてDifferentiation IVによる結果の3通りについて比較する。

# まず被説明変数を定義する。
data <- data %>%
  dplyr::mutate(logit_share = log(share) - log(share0))

# OLSの結果
ols <- fixest::feols(logit_share ~ price + hppw + FuelEfficiency + size, data = data)

# BLP操作変数を用いた結果
iv_BLP <- fixest::feols(
  logit_share ~ hppw + FuelEfficiency + size | 0 |
    price ~ iv_BLP_own_hppw + iv_BLP_own_FuelEfficiency + iv_BLP_own_size +
    iv_BLP_other_hppw + iv_BLP_other_FuelEfficiency + iv_BLP_other_size,
  data = data
)

# Differentiation IVを用いた結果
iv_GH <- fixest::feols(
  logit_share ~ hppw + FuelEfficiency + size | 0 |
    price ~ iv_GH_own_hppw + iv_GH_own_FuelEfficiency + iv_GH_own_size +
    iv_GH_other_hppw + iv_GH_other_FuelEfficiency + iv_GH_other_size,
  data = data
)

# 推定結果をレポートする。
tbl_logit_iv <- fixest::etable(
  list("OLS" = ols, "iv_BLP" = iv_BLP, "iv_GH" = iv_GH),  
  se = "hetero",
  fitstat = c("r2", "n", "ivf" ), 
  signif.code = NA, 
  dict = c(
    price = "自動車価格",
    hppw = "馬力／重量",
    FuelEfficiency = "燃費(キロメートル／ 1 リットル)",
    size = "サイズ",
    `(Intercept)` = "定数項"
  ),
  digits = 2,
  digits.stats = 2,
  depvar = FALSE) %>% {
    txt <- capture.output(knitr::kable(., format = "simple"))
    writeLines(txt, con = here("02_BLP_CH03_04_05/output/tab3_2_tbl_logit_iv.txt"), useBytes = TRUE)
  }


## 1st stage ----

# BLP操作変数を用いた結果
iv1st_BLP <- fixest::feols(
  price ~ hppw + FuelEfficiency + size +
    iv_BLP_own_hppw + iv_BLP_own_FuelEfficiency + iv_BLP_own_size +
    iv_BLP_other_hppw + iv_BLP_other_FuelEfficiency + iv_BLP_other_size,
  data = data
)

# Differentiation IVを用いた結果
iv1st_GH <- fixest::feols(
  price ~ hppw + FuelEfficiency + size +
    iv_GH_own_hppw + iv_GH_own_FuelEfficiency + iv_GH_own_size +
    iv_GH_other_hppw + iv_GH_other_FuelEfficiency + iv_GH_other_size,
  data = data
)

# 推定結果をレポートする。
tbl_logit_iv_1st <- fixest::etable(
  list("iv1st_BLP" = iv1st_BLP, "iv1st_GH" = iv1st_GH),  
  se = "hetero",
  fitstat = c("r2", "n"), 
  signif.code = NA, 
  dict = c(
    price = "自動車価格",
    hppw = "馬力／重量",
    FuelEfficiency = "燃費(キロメートル／ 1 リットル)",
    size = "サイズ",
    `(Intercept)` = "定数項"
  ),
  digits = 2,
  digits.stats = 2,
  depvar = FALSE
  ) %>% {
    txt <- capture.output(knitr::kable(., format = "simple"))
    writeLines(txt, con = here("02_BLP_CH03_04_05/output/tbl_logit_iv_1st.txt"), useBytes = TRUE)
}

## 自己価格弾力性の計算 ----
data_elas <- data %>% 
  dplyr::mutate(
    own_elas_ols = ols$coefficients["price"] * price * (1 - share), 
    own_elas_ivblp = iv_BLP$coefficients["fit_price"] * price * (1 - share), 
    own_elas_ivgh  = iv_GH$coefficients["fit_price"] * price * (1 - share)
  )

tbl_own_elas <- data_elas %>% 
  dplyr::select(starts_with("own_elas")) %>% 
  summarytools::descr(
    transpose = TRUE, 
    stats = c("mean", "sd", "med", "min", "max"),
    order = "preserve"
  ) %>% {
    txt <- capture.output(knitr::kable(., format = "simple"))
    writeLines(txt, con = here("02_BLP_CH03_04_05/output/tab3_3_own_elas.txt"), useBytes = TRUE)
  }

# 推定結果の応用 ----

## 需要曲線と収入曲線を描く ----

# 需要関数の推定結果にもとづいて、需要曲線を書こう。
# 1.  まず、推定結果から$\xi_{jt}$をゲットする。
# 2.  関数を作成する。
# 3.  どこか市場を固定する。
# 4.  製品を一つ固定する。
# 5.  その製品の価格を変えたときに、どの製品のSalesがどう変わるかをPredictする。その際、他の製品価格はデータのものに固定しておく。


dt_application <- data %>% 
  dplyr::select(NameID, year, Sales, price, FuelEfficiency, size, hppw, HH, share) %>% 
  dplyr::mutate(xi_fit = resid(iv_GH))

# NameID==87という特定の車種に着目する

NameID_target <- data %>%
  dplyr::filter(
    Nippyo == 1,
    Maker == "Toyota",
    Name == "アルファード"
  ) %>%
  dplyr::distinct(NameID) %>%
  dplyr::pull()

dt_application %>% 
  dplyr::filter(NameID == NameID_target & year == 2016) %>% 
  head()

estparam <- iv_GH$coefficients

# 価格の範囲として、0.3(30万円)から5(500万円)を考える。なお、もとの価格は3.198 (319.8万円)
pricevec <- seq(from = 0.3, to = 5,　by = 0.05)
quantvec <- numeric(length(pricevec))
revenuevec <- numeric(length(pricevec))

for (i in 1:length(pricevec)) {
  quantvec[i] <- calculate_sales(
    price_cand = pricevec[i], year = 2016,
    modelID_target = NameID_target, dt = dt_application, est_param = estparam
  )
  revenuevec[i] <- pricevec[i] * quantvec[i]
}

# 需要曲線をプロット
data.frame(x = quantvec, y = pricevec) %>%
  ggplot2::ggplot(aes(x = x / 10000, y = y * 100)) +
  ggplot2::geom_line() +
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
  ggplot2::scale_x_continuous(breaks = seq(0, 20, 2)) + 
  ggplot2::scale_y_continuous(breaks = seq(0, 500, 100)) + 
  ggplot2::labs(
    y = "価格(万円)",
    x = "販売台数(万台)"
  ) 

ggsave(file = here("02_BLP_CH03_04_05/output/fig3_1_demand.pdf"), width = 9, height = 6, device = cairo_pdf)

# 収入と価格の関係をプロット
data.frame(x = pricevec, y = revenuevec) %>%
  ggplot2::ggplot(aes(x = x * 100, y = y * 100 / 10000)) +
  ggplot2::geom_line() +
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
  ggplot2::scale_x_continuous(breaks = seq(0, 500, 100)) + 
  ggplot2::scale_y_continuous(breaks = seq(500, 1500, 250)) + 
  ggplot2::labs(
    y = "収入(億円)",,
    x = "価格(万円)"
  ) 

ggsave(file = here("02_BLP_CH03_04_05/output/fig3_1_revenue.pdf"), width = 9, height = 6, device = cairo_pdf)


# 最適化アルゴリズムを使って、収入を最大にする価格を求める。
result <- optimise(
  maximize_revenue, 
  interval = c(0.3, 3), 
  maximum = TRUE, 
  year = 2016, 
  modelID_target = NameID_target,
  dt = dt_application, 
  est_param = estparam
)

print(result)
