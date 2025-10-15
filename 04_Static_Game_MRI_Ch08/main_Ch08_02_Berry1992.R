## ----setup, include=FALSE------------------------------------------------
# チャンクの設定
knitr::opts_chunk$set(echo = TRUE,
                      warning = FALSE)


## ------------------------------------------------------------------------
# ワークスペースを初期化
rm(list = ls())

# パッケージを読み込む
require(tidyverse)
require(tictoc)
require(skimr)
require(modelsummary)
require(kableExtra)
require(doParallel)
require(doRNG)


## ------------------------------------------------------------------------
# 病院データを読み込む
data_raw <- readr::read_csv(file = "data_hospital_chapter6.csv", 
                            locale = locale(encoding = "UTF-8"))
  


## ------------------------------------------------------------------------
# ダミー変数を作成
data_cleaned <- data_raw %>% 
  # 大学（公立、国立、私立大学）ダミー
  dplyr::mutate(DaigakuDum = if_else(
    Management %in% c('公立大学法人','国（国立大学法人）','私立学校法人'),1, 0)) %>% 
  # 病床数0ダミー
  dplyr::mutate(ZeroBedDum = if_else(NumBeds == 0, 1, 0))


## ------------------------------------------------------------------------
# 単位の調整と、対数変換した変数の作成
data_cleaned <- data_cleaned %>% 
  dplyr::mutate(
    NumBeds = NumBeds/100,
    LogNumBeds = log(NumBeds+0.01),
    Population = Population/10^6,
    Menseki = Menseki/100,
    TaxableIncome = TaxableIncome/1000,
    LogPop = log(Population),
    LogIncome = log(TaxableIncome))


## ------------------------------------------------------------------------
data_cleaned %>% 
  dplyr::group_by(MRIOwnDum) %>% 
  dplyr::select(MRIOwnDum, Kyukyu, Kinou, Sien, Hyoka, 
         DepNeurology ,DepNeurosurgery, NumBeds,
         ZeroBedDum, DaigakuDum) %>% 
  skimr::skim() %>% 
  skimr::yank("numeric") %>% 
  dplyr::select(skim_variable, MRIOwnDum, mean, sd) %>% 
  as.data.frame() %>% 
  tidyr::pivot_wider(names_from = MRIOwnDum,
                     values_from = c(mean, sd)) %>% 
  dplyr::select(skim_variable, mean_1, sd_1, mean_0, sd_0) %>% 
  dplyr::mutate(skim_variable = c('救急', '機能', '支援', 
                                  '評価', '神経内科', '神経外科', 
                                  '病床数', '病床0', '大学病院')) %>% 
  kableExtra::kable(format = "html", 
                    booktabs = TRUE,
                    caption = "記述統計", 
                    col.names = c('変数', '平均', '標準偏差',
                                  '平均', '標準偏差'),
                    digits = 3) %>% 
  kableExtra::add_header_above(
    c(" ", "MRIを購買した病院" = 2, 
      "MRIを購買していない病院" = 2)) %>%
  kableExtra::kable_styling() 


## ------------------------------------------------------------------------
# 以下の表を作成するための関数
per_MRI_cal <- function(df){
  df %>% 
    dplyr::summarise(Total = n(),
                     MRIHos = sum(MRIOwnDum),
                     PerMRI = round(MRIHos / Total * 100, 2))
}

# 病院のカテゴリー毎に、総病院数、MRI所有病院数、その割合を確認
list(data_cleaned,
     data_cleaned %>% dplyr::filter(Kyukyu == 1),
     data_cleaned %>% dplyr::filter(Sien == 1),
     data_cleaned %>% dplyr::filter(Hyoka == 1),
     data_cleaned %>% dplyr::filter(DepNeurology == 1),
     data_cleaned %>% dplyr::filter(DepNeurosurgery == 1),
     # 元データの上位25%が1.2
     data_cleaned %>% dplyr::filter(LogNumBeds >= log(1.2)),
     data_cleaned %>% dplyr::filter(DaigakuDum == 1),
     data_cleaned %>% dplyr::filter(ZeroBedDum == 1)
     ) %>%
  purrr::map_dfr(.f =  per_MRI_cal) %>% 
  dplyr::mutate(Category = 
                  c('All', 'Kyukyu', 'Sien', 
                    'Hyoka', 'DepNeurology', 'DepNeurosurgery', 
                    'Top25%_NumBeds', 'DaigakuDum', 'ZeroBedDum')) %>% 
  dplyr::relocate(Category, .before = Total) %>% 
  kableExtra::kable(format = "html", 
                    booktabs = TRUE) %>% 
  kableExtra::kable_styling() 


## ------------------------------------------------------------------------
# Step 0: 欠損のあるデータを削除
data_processed_pre <- 
  data_cleaned %>% 
  # 必要な変数のみ抽出
  dplyr::select(CityCode, Kyukyu, Kinou, Sien, Hyoka,
         DepNeurology, DepNeurosurgery, LogNumBeds, 
         ZeroBedDum, DaigakuDum,
         Menseki, LogPop, LogIncome, 
         MRIOwnDum) %>% 
  # 欠損を含む観察を削除
  na.omit()


## ------------------------------------------------------------------------
# Step 1: シミュレーションに必要な変数を作成、参入の順番に並び替え
berry_process_data <- function(df){
  data_processed <- 
    df %>% 
    dplyr::arrange(CityCode) %>% 
    # 各病院に乱数を振る
    dplyr::mutate(TieEntryOrder = runif(nrow(df))) %>% 
    # 自治体毎の潜在的参入病院数、観察された参入病院数の変数を作成
    #（本誌のIm, Nmに対応）
    dplyr::group_by(CityCode) %>% 
    dplyr::mutate(NumPotenHos = as.numeric(n()),
                  NumEntryObs = sum(MRIOwnDum)) %>% 
    # データを参入の順番に並び替える
    # 病床数が多い順に並び替える。タイの場合は乱数を参照。
    # 注意: 参入順序の仮定は推定には不要だが、モデルによる予測や反実仮想分析では
    # 必要となるため、前もって初期の並び順を指定している。
    dplyr::arrange(desc(LogNumBeds), desc(TieEntryOrder),
                 .by_group = TRUE) %>% 
    # incumbents move first の仮定の場合
    # dplyr::arrange(desc(MRIOwnDum), desc(LogNumBeds), desc(TieEntryOrder),
    #                .by_group = TRUE) %>% 
    # 自治体毎に、参入した順番に沿った病院idを振る
    dplyr::mutate(EntryOrderId = as.numeric(1:NumPotenHos)) %>% 
    dplyr::ungroup()
  return(data_processed)
}


## ------------------------------------------------------------------------
# 以降の結果が変わらないように、乱数を固定
set.seed(1)
# Step 1で作成した関数の実行
data_processed <- berry_process_data(data_processed_pre)
# 不要なオブジェクトを削除
rm(data_processed_pre)


## ------------------------------------------------------------------------
# Step 2: 潜在的参入病院数を絞ってサブサンプル
NumPotenHos_max <- 4
data_processed <- data_processed %>%
  # NumEntryObs_max 以下の観察に絞る
  dplyr::filter(EntryOrderId <= NumPotenHos_max) %>%
  # NumEntryObs, NumPotenHos の値を更新
  dplyr::group_by(CityCode) %>% 
  dplyr::mutate(NumEntryObs = sum(MRIOwnDum),
                NumPotenHos = as.numeric(n())) %>% 
  dplyr::ungroup()


## ------------------------------------------------------------------------
# 病院のカテゴリー毎に、総病院数、MRI所有病院数、その割合を確認
list(data_processed,
     data_processed %>% dplyr::filter(Kyukyu == 1),
     data_processed %>% dplyr::filter(Sien == 1),
     data_processed %>% dplyr::filter(Hyoka == 1),
     data_processed %>% dplyr::filter(DepNeurology == 1),
     data_processed %>% dplyr::filter(DepNeurosurgery == 1),
     # 元データの上位25%が1.2
     data_processed %>% dplyr::filter(LogNumBeds >= log(1.2)),
     data_processed %>% dplyr::filter(DaigakuDum == 1),
     data_processed %>% dplyr::filter(ZeroBedDum == 1)
     ) %>%
  purrr::map_dfr(.f =  per_MRI_cal) %>% 
  dplyr::mutate(Category = 
                  c('All', 'Kyukyu', 'Sien', 
                    'Hyoka', 'DepNeurology', 'DepNeurosurgery', 
                    'Top25%_NumBeds', 'DaigakuDum', 'ZeroBedDum')) %>% 
  dplyr::relocate(Category, .before = Total) %>% 
  kableExtra::kable(format = "html", 
                    booktabs = TRUE) %>% 
  kableExtra::kable_styling() 


## ------------------------------------------------------------------------
# Step 3: シミュレーションを行いやすくするためにデータを縦に長く加工
# データの行数は \Sigma_{m = 1}^{M} (Im * Im) * nsとなる
berry_expand_data <- function(df, ns){
  # 自治体の数
  M <- df %>% 
    dplyr::select(CityCode) %>% 
    dplyr::n_distinct() 
  # 病院の数（本誌のImのmについての和に対応）
  NumHos <- df %>%
    dplyr::distinct(CityCode, NumPotenHos) %>%
    dplyr::summarise(sum(NumPotenHos)) %>%
    as.numeric()
  
  # 誤差項をシミュレーション
  u_m0 <- rnorm(M * ns, mean = 0, sd = 1)
  u_mIm <- rnorm(NumHos * ns, mean = 0, sd = 1)
  u_df <- df %>% 
    # 自治体、病院idのみのdataframeを作成
    dplyr::select(CityCode, EntryOrderId) %>% 
    dplyr::group_by(CityCode) %>% 
    nest() %>% 
    # シミュレーション番号を振る（列数はns倍）
    tidyr::crossing(s = 1:ns) %>% 
    # 自治体毎の誤差項の変数を作成
    dplyr::mutate(u_m0) %>% 
    tidyr::unnest(cols = c(data)) %>% 
    # 病院毎の誤差項の変数を作成 
    dplyr::mutate(u_mIm)
  # 誤差項のdataframeを元のdataframeに統合
  data_expand <- df %>% 
    dplyr::left_join(u_df, by = c('CityCode', 'EntryOrderId'))
  data_expand
  # シュミレーションで用いる、市場に存在する企業数の変数を作成（本誌のNに対応）
  # 自治体毎、病院毎に1からNumPotenHosまでの数を振るため、データが縦に長くなることに注意。
  num_entry_comb_df <-
    df %>% 
    group_by(CityCode, EntryOrderId) %>%
    summarise(n_cand = 1:NumPotenHos, 
              .groups = 'drop')
  data_expand <-
    data_expand %>% 
    dplyr::left_join(num_entry_comb_df, 
                     by = c('CityCode', 'EntryOrderId'))
  
  return(data_expand)
}


## ------------------------------------------------------------------------
# Step 3 で作成した関数の実行
ns <- 100
data_expand <- 
  berry_expand_data(
    data_processed, 
    ns)


## ------------------------------------------------------------------------
# データとパラメターを所与としたBerry (1992) の目的関数の値を計算する関数を作成
berry_obj <- function(params, df){
  # パラメター
  alpha <- params[1:8]
  beta <- params[9:12]
  delta <- params[13]
  rho <- params[14]
  # 可変利潤、固定費用の作成（本誌のvとphiに対応）
  df <- df %>% 
    dplyr::mutate(
      var_profit = 
        beta[1] + beta[2] * Menseki + 
        beta[3] * LogPop + beta[4] * LogIncome +
        - delta * log(n_cand) + 
        rho * u_m0,
      fixed_cost = 
        alpha[1] * Kyukyu + 
        alpha[2] * Sien + alpha[3] * Hyoka +
        alpha[4] * DepNeurology + alpha[5] * DepNeurosurgery + 
        alpha[6] * LogNumBeds + alpha[7] * ZeroBedDum + 
        alpha[8] * DaigakuDum +
        - sqrt(1 - rho ^2) * u_mIm
    ) 
  # Berry (1992) による手法の目的関数を計算する
  obj <- df %>% 
    # Berry (1992) equation (12)に基づいて期待参入企業数を計算
    # 自治体、s、病院id、について、
    # 参入するかどうかの論理値型の列を作成
    dplyr::mutate(EntryDecision = var_profit > fixed_cost) %>% 
    dplyr::group_by(CityCode, NumEntryObs, s, n_cand) %>% 
    # 参入する企業数を計算
    dplyr::summarise(n_sim = sum(EntryDecision), 
                     .groups = 'drop') %>% 
    # 均衡条件を満たすNの候補を計算
    dplyr::mutate(eq_cond = ifelse(n_sim >= n_cand, 1, 0),
                  n_star = n_cand * eq_cond) %>% 
    # 条件を満たす中で最も大きい値をN*とする
    dplyr::group_by(CityCode, NumEntryObs, s) %>%
    dplyr::summarise(n_star = max(n_star), 
                     .groups = 'drop') %>% 
    # 自治体について、期待参入企業数を計算
    dplyr::group_by(CityCode, NumEntryObs) %>%
    dplyr::summarise(n_exp = mean(n_star), 
                     .groups = 'drop') %>% 
    # 実際の参入病院数と比較し、誤差の二乗を計算
    dplyr::summarise(mean((NumEntryObs - n_exp)^2)) %>% 
    as.numeric()
  return(obj)
}


## ------------------------------------------------------------------------
# 係数の名前
est_names <- 
  c('Kyukyu', 'Sien', 'Hyoka',
    'DepNeurology', 'DepNeurosurgery', 'LogNumBeds',
    'ZeroBedDum', 'DaigakuDum', 'Constant',
    'Menseki', 'LogPop', 'LogIncome',
    'delta', 'rho')


## ------------------------------------------------------------------------
# 初期値
param_init <- 
  c(-0.587617,-1.298752,   -0.556557,
    -0.502602,   -1.058680,   -0.919110,  
    -2.832720,   -2.726932,    1.937292,
    0.009940,    0.206644,   -1.617755,
    0.550239,    0.132202)
# Obj: 0.3666109
# berry_obj(params = param_init,
#           df = data_expand)


## ------------------------------------------------------------------------
tic.clearlog() 
tic()
berry_result <- 
  optim(params <- param_init,
        berry_obj, 
        method = "Nelder-Mead",
        df = data_expand)
# rhoの符号を正しくする
berry_result$par <- berry_result$par %>% 
  setNames(est_names)
berry_result$par[['rho']] <- abs(berry_result$par[['rho']])
berry_result$par <- berry_result$par %>% as.numeric()
toc()


## ----eval=FALSE----------------------------------------------------------
## # 並列計算のための設定を行う。
## cores <- 20 #並列計算に使うCPUコアの個数。
## registerDoParallel(cores)
## registerDoRNG(12345)
## 


## ----eval=FALSE----------------------------------------------------------
## # ブートストラップの回数
## iter <- 100
## # 自治体の数
## M <- data_processed %>%
##   dplyr::select(CityCode) %>%
##   dplyr::n_distinct()
## tic()
## # 繰り返し推定を行う
## bs_berry_result <-
##   foreach::foreach (
##     i = 1:iter,
##     .export = ls(envir=parent.frame()),
##       .packages = 'tidyverse'
##     ) %dopar% {
##       # リサンプルを行う
##       bs_data <- data_processed %>%
##         dplyr::distinct(CityCode) %>%
##         dplyr::slice_sample(n = M, replace = TRUE) %>%
##         dplyr::arrange(CityCode) %>%
##         dplyr::mutate(id = 1:M) %>%
##         dplyr::left_join(data_processed, by = 'CityCode') %>%
##         # 推定ではCityCodeにより自治体を判別しているため名前を修正
##         dplyr::rename(CityCodeOrigin = CityCode, CityCode = id)
##       # リサンプルしたデータを推定のために加工
##       bs_data_processed <-
##         berry_process_data(bs_data)
##       # 先程と同様に推定のためにデータを縦に広げる
##       bs_data_expand <-
##         berry_expand_data(
##           bs_data_processed,
##           ns)
##       # 先程と同様に最適化により推定
##       bs_berry_est <-
##         optim(params <- param_init,
##               berry_obj,
##               method = "Nelder-Mead",
##               df = bs_data_expand
##               )$par
##       # rhoの符号を正しくする
##       bs_berry_est <- bs_berry_est %>%
##         setNames(est_names)
##       bs_berry_est[['rho']] <- abs(bs_berry_est[['rho']])
##       return(bs_berry_est)
##   }
## toc()
## 
## # 結果を保存する。
## saveRDS(bs_berry_result, file = "BootStrap_Berry_Result.rds")
## 


## ------------------------------------------------------------------------
bs_berry_result <- readRDS("BootStrap_Berry_Result.rds")

# 標準誤差を計算
berry_se <- bs_berry_result %>%
  map_dfr(~ unlist(as.data.frame(.))) %>%
  summarise(across(everything(), sd)) %>%
  as.numeric()



## ----message=FALSE-------------------------------------------------------
# 推定結果
tibble(
  name = c(est_names, 'Value of obj'),
  BR_est = c(berry_result$par, berry_result$value),
  # BR_se = berry_result_1129$BR_se
  BR_se = c(berry_se, NA)
  ) %>% 
  dplyr::mutate(name = c('救急', '支援', '評価', '神経内科', '神経外科', 
                         'log(病床数)', '病床0', '大学病院',
                         '定数項', '面積', 'log(人口)', 'log(課税所得)',
                         'delta', 'rho', '目的関数の値')) %>% 
  kableExtra::kable(format = "html", 
                    booktabs = TRUE,
                    caption = "Berry (1992) による推定結果", 
                    col.names = c('変数', '推定値', '標準誤差'),
                    digits = 3) %>%
  kable_styling()


## ------------------------------------------------------------------------
# モデルによる参入予測確率を含めたdataframeを作成
predict_model <- function(df){
  data_predicted_pre <- 
    df %>% 
    # Berry (1992) equation (12)に基づいて期待参入企業数を計算
    # 自治体、s、病院id、について、
    # 参入するかどうかの論理値型の列を作成
    dplyr::mutate(EntryDecision = var_profit > fixed_cost) %>% 
    dplyr::group_by(CityCode, NumEntryObs, s, n_cand) %>% 
    # 参入する企業数を計算
    dplyr::summarise(n_sim = sum(EntryDecision), 
                     .groups = 'drop') %>% 
    # 均衡条件を満たすNの候補を計算
    dplyr::mutate(eq_cond = ifelse(n_sim >= n_cand, 1, 0),
                  n_star = n_cand * eq_cond) %>% 
    # 条件を満たす中で最も大きい値をN*とする
    dplyr::group_by(CityCode, NumEntryObs, s) %>%
    dplyr::summarise(n_star = max(n_star), 
                     .groups = 'drop') %>% 
    dplyr::select(-NumEntryObs) %>% 
    # 元のdata.frameにn_starを含むdata.frameをjoin
    dplyr::right_join(
      df %>% 
        dplyr::select(-n_cand) %>% 
        dplyr::distinct(across(everything())),
      by = c('CityCode', 's')
    ) %>% 
    # 参入順序の仮定に基づいて、それぞれの企業が参入するかどうかを決める
    dplyr::mutate(EntryCond = ifelse(n_star >= EntryOrderId, 1, 0)) %>% 
    # 参入の予測確率を計算
    dplyr::group_by(CityCode, EntryOrderId) %>% 
    dplyr::summarise(EntryProb = mean(EntryCond),
                     .groups = 'drop') %>% 
    dplyr::mutate(EntryPred = ifelse(EntryProb >= 0.5 ,1 ,0))
  
  # 各企業が参入するかどうかを示す列を含むdata.frameを、元のdata.frameにmergeする
  data_predicted <-
    df %>% 
    dplyr::select(-c(s, u_m0, u_mIm, n_cand, var_profit, fixed_cost)) %>% 
    dplyr::distinct(across(everything())) %>% 
    dplyr::left_join(data_predicted_pre, 
                     by = c('CityCode', 'EntryOrderId'))
  return(data_predicted)
}


## ------------------------------------------------------------------------
# 推定したパラメターを取得
berry_est <- berry_result$par
# パラメター
alpha <- berry_est[1:8]
beta <- berry_est[9:12]
delta <- berry_est[13]
rho <- berry_est[14]

data_expand <-
  data_expand %>% 
  # 可変利潤、固定費用の作成（本誌のvとphiに対応）
  dplyr::mutate(
    var_profit = 
      beta[1] + beta[2] * Menseki + 
      beta[3] * LogPop + beta[4] * LogIncome +
      - delta * log(n_cand) + 
      rho * u_m0,
    fixed_cost = 
      alpha[1] * Kyukyu  + 
      alpha[2] * Sien + alpha[3] * Hyoka +
      alpha[4] * DepNeurology + alpha[5] * DepNeurosurgery + 
      alpha[6] * LogNumBeds + alpha[7] * ZeroBedDum + 
      alpha[8] * DaigakuDum +
      - sqrt(1 - rho ^2) * u_mIm,
  ) 


## ------------------------------------------------------------------------
data_predicted <- 
  data_expand %>% 
  predict_model()


## ------------------------------------------------------------------------
# 参入病院数毎に、自治体がいくつ観察されるかを集計
data_predicted %>% 
  dplyr::group_by(CityCode) %>% 
  dplyr::summarise(Actual = sum(MRIOwnDum),
                   Predict = sum(EntryPred)) %>% 
  tidyr::pivot_longer(cols = c('Actual', 'Predict')) %>% 
  ggplot2::ggplot(aes(x = value, fill = name)) + 
  ggplot2::geom_histogram(binwidth = 0.5,
                          position = "dodge") + 
  coord_flip() + 
  scale_x_reverse()


## ------------------------------------------------------------------------
# 以下の表を作成するための関数
actual_pred_cal <- function(df){
  df %>% 
    dplyr::summarise(Actual = sum(MRIOwnDum),
                     Predict = sum(EntryPred))
}


## ------------------------------------------------------------------------
# 病院のタイプ別の観察された参入病院数と予測された参入病院数
list(data_predicted,
     data_predicted %>% dplyr::filter(Kyukyu == 1),
     data_predicted %>% dplyr::filter(Sien == 1),
     data_predicted %>% dplyr::filter(Hyoka == 1),
     data_predicted %>% dplyr::filter(DepNeurology == 1),
     data_predicted %>% dplyr::filter(DepNeurosurgery == 1),
     # 元データの上位25%が1.2
     data_predicted %>% dplyr::filter(LogNumBeds >= log(1.2)),
     data_predicted %>% dplyr::filter(DaigakuDum == 1),
     data_predicted %>% dplyr::filter(ZeroBedDum == 1)
     ) %>%
  purrr::map_dfr(.f =  actual_pred_cal) %>% 
  dplyr::mutate(Category = 
                  c('All', 'Kyukyu', 'Sien', 
                    'Hyoka', 'DepNeurology', 'DepNeurosurgery', 
                    'Top25%_NumBeds', 'DaigakuDum', 'ZeroBedDum')) %>% 
  tidyr::pivot_longer(cols = c('Actual', 'Predict')) %>%
  ggplot2::ggplot(aes(x = Category, y = value, fill = name)) + 
  geom_bar(position="dodge", stat = "identity") + 
  coord_flip() 
  # theme_minimal()


## ------------------------------------------------------------------------
# 特定の病院について平均二乗誤差を表示
tibble(
  category = c('Total', 'DepNeurology', 'DepNeurosurgery', 'ZeroBedDum', 'DaigakuDum'),
  mse = c(data_predicted %>% 
            dplyr::summarise(mean((MRIOwnDum - EntryProb) ^ 2)) %>% 
            as.numeric() %>% round(3),
          data_predicted %>% 
            dplyr::filter(DepNeurology == 1) %>% 
            dplyr::summarise(mean((MRIOwnDum - EntryProb) ^ 2)) %>% 
            as.numeric() %>% round(3),
          data_predicted %>% 
            dplyr::filter(DepNeurosurgery == 1) %>% 
            dplyr::summarise(mean((MRIOwnDum - EntryProb) ^ 2)) %>% 
            as.numeric() %>% round(3),
          data_predicted %>% 
            dplyr::filter(ZeroBedDum == 1) %>% 
            dplyr::summarise(mean((MRIOwnDum - EntryProb) ^ 2)) %>% 
            as.numeric() %>% round(3),
          data_predicted %>% 
            dplyr::filter(DaigakuDum == 1) %>% 
            dplyr::summarise(mean((MRIOwnDum - EntryProb) ^ 2)) %>% 
            as.numeric() %>% round(3)
          )
)


## ------------------------------------------------------------------------
# 神経内科と神経外科を持つ病院はMRIの所有を義務付けられたとする
# MRI所有義務化病院が既に参入しているとして、他の病院が参入するかどうかをシミュレート
data_cf <-  
  data_expand %>%
  dplyr::group_by(CityCode, s, n_cand) %>% 
  # MRI所有義務化病院が先になるように並び替える
  dplyr::arrange(desc(DepNeurology), desc(DepNeurosurgery), 
                 desc(LogNumBeds), desc(TieEntryOrder),
                 .by_group = TRUE) %>% 
  # 参入する順番を更新する
  dplyr::mutate(EntryOrderId = as.numeric(1:NumPotenHos)) %>% 
  dplyr::ungroup() %>% 
  # 反実仮想分析の処理: MRI所有義務化病院が必ず参入するために大学病院の利潤を大きな値に置き換え
  dplyr::mutate(var_profit = 
                  ifelse(DepNeurology == 1 | DepNeurosurgery == 1
                         , 10^5, var_profit)) %>% 
  predict_model()


## ------------------------------------------------------------------------
# 以下の表を作成するための関数
per_decline_cal <- function(df){
  df %>% 
    dplyr::summarise(
      Total = n(),
      Actual = sum(MRIOwnDum),
      Fit_Pred = sum(Fit_Pred),
      CF_Pred = sum(EntryPred),
      CF_Pred_am_actual = sum(MRIOwnDum * EntryPred)
      )
}


## ------------------------------------------------------------------------
# 病院のカテゴリー毎に、総病院数(Total)、
# 元のモデルにより予測されたMRI所有病院数(Fit_Pred)、
# 反実仮想分析により予測されたMRI所有病院数(CF_Pred)、
# MRI所有病院数の変化率(PerChange)、を確認
generate_cf_result <- function(data_predicted, data_cf){
  data_cf_sub <- 
    data_cf %>% 
    # 元のモデルでの参入の予測を表す変数を加える
    dplyr::left_join(data_predicted %>% 
                       dplyr::select(CityCode, TieEntryOrder, EntryPred) %>% 
                       dplyr::rename(Fit_Pred = EntryPred),
                     by = c("CityCode", "TieEntryOrder")) %>%
    # 政策の対象かどうかを表すダミー
    dplyr::mutate(PolicyTarget = 
                    ifelse(
                      (DepNeurology == 1 | DepNeurosurgery == 1) & Fit_Pred != 1,
                      #DepNeurology == 1 | DepNeurosurgery == 1) & MRIOwnDum != 1,
                      1, 0)
    ) %>% 
    dplyr::filter(PolicyTarget == 0)
  
  list(data_cf_sub,
       data_cf_sub %>% dplyr::filter(Kyukyu == 1),
       data_cf_sub %>% dplyr::filter(Sien == 1),
       data_cf_sub %>% dplyr::filter(Hyoka == 1),
       data_cf_sub %>% dplyr::filter(DepNeurology == 1),
       data_cf_sub %>% dplyr::filter(DepNeurosurgery == 1),
       data_cf_sub %>% dplyr::filter(!(DepNeurology == 1 | DepNeurosurgery == 1)),
       # 元データの上位25%が1.2
       data_cf_sub %>% dplyr::filter(LogNumBeds >= log(1.2)),
       data_cf_sub %>% dplyr::filter(DaigakuDum == 1),
       data_cf_sub %>% dplyr::filter(ZeroBedDum == 1)
  ) %>%
    purrr::map_dfr(.f =  per_decline_cal) %>% 
    dplyr::mutate(Category = 
                    c('All', 'Kyukyu', 'Sien', 'Hyoka',
                      'DepNeurology', 'DepNeurosurgery', 
                      'Not_(DepNeurology_and_DepNeurosurgery)',
                      'Top25%_NumBeds', 'DaigakuDum', 'ZeroBedDum')) %>% 
    dplyr::select(Category, Fit_Pred, CF_Pred) %>% 
    dplyr::mutate(PerChange = round((CF_Pred/Fit_Pred - 1) * 100, 2)) %>%
    kableExtra::kable(format = "html",
                      booktabs = TRUE,
                      digits = 3) %>%
    kableExtra::kable_styling() 
}


## ------------------------------------------------------------------------
generate_cf_result(data_predicted, data_cf)


## ------------------------------------------------------------------------
data_predicted_order2 <- 
  data_expand %>% 
  dplyr::group_by(CityCode, s, n_cand) %>% 
  # 既存企業、病床数の順に並び替え
  dplyr::arrange(desc(MRIOwnDum), 
                 desc(LogNumBeds), desc(TieEntryOrder),
                 .by_group = TRUE) %>% 
  # 参入する順番を更新する
  dplyr::mutate(EntryOrderId = as.numeric(1:NumPotenHos)) %>% 
  dplyr::ungroup() %>% 
  predict_model()


## ------------------------------------------------------------------------
data_cf_order2 <-  
  data_expand %>%
  dplyr::group_by(CityCode, s, n_cand) %>% 
  # MRI所有義務化病院が先になるように並び替える
  dplyr::arrange(desc(DepNeurology), desc(DepNeurosurgery), 
                 desc(MRIOwnDum), 
                 desc(LogNumBeds), desc(TieEntryOrder),
                 .by_group = TRUE) %>% 
  # 参入する順番を更新する
  dplyr::mutate(EntryOrderId = as.numeric(1:NumPotenHos)) %>% 
  dplyr::ungroup() %>% 
  # 反実仮想分析の処理: MRI所有義務化病院が必ず参入するために大学病院の利潤を大きな値に置き換え
  dplyr::mutate(var_profit = 
                  ifelse(DepNeurology == 1 | DepNeurosurgery == 1
                         , 10^5, var_profit)) %>% 
  predict_model()


## ------------------------------------------------------------------------
generate_cf_result(data_predicted_order2, data_cf_order2)


## ------------------------------------------------------------------------
random_entry_df <-
  data_expand %>% 
  dplyr::distinct(CityCode, EntryOrderId) %>%
  mutate(RandomEntryOrder = runif(n(),0,1))


## ------------------------------------------------------------------------
data_predicted_order3 <- 
  data_expand %>% 
  dplyr::left_join(random_entry_df, 
                   by = c('CityCode', 'EntryOrderId')) %>% 
  dplyr::group_by(CityCode, s, n_cand) %>% 
  # ランダムに並び替え
  dplyr::arrange(desc(RandomEntryOrder),
                 .by_group = TRUE) %>% 
  # 参入する順番を更新する
  dplyr::mutate(EntryOrderId = as.numeric(1:NumPotenHos)) %>% 
  dplyr::ungroup() %>% 
  predict_model()


## ------------------------------------------------------------------------
data_cf_order3 <-  
  data_expand %>%
  dplyr::left_join(random_entry_df, 
                   by = c('CityCode', 'EntryOrderId')) %>% 
  dplyr::group_by(CityCode, s, n_cand) %>% 
  # MRI所有義務化病院が先になるように並び替える
  dplyr::arrange(desc(DepNeurology), desc(DepNeurosurgery), 
                 desc(RandomEntryOrder),
                 .by_group = TRUE) %>% 
  # 参入する順番を更新する
  dplyr::mutate(EntryOrderId = as.numeric(1:NumPotenHos)) %>% 
  dplyr::ungroup() %>% 
  # 反実仮想分析の処理: MRI所有義務化病院が必ず参入するために大学病院の利潤を大きな値に置き換え
  dplyr::mutate(var_profit = 
                  ifelse(DepNeurology == 1 | DepNeurosurgery == 1
                         , 10^5, var_profit)) %>% 
  predict_model()


## ------------------------------------------------------------------------
generate_cf_result(data_predicted_order3, data_cf_order3)

