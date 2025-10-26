## 全体の流れ
# 1.  Rに関する下準備
# 2.  データの準備
# 3.  企業レベルの変数の記述統計
# 4.  Berry (1992)による推定  
# 5.  推定結果  
# 6.  反実仮想分析


# 1.  Rに関する下準備 ----

# ワークスペースを初期化
rm(list = ls())

# パッケージを読み込む
library(tidyverse)
library(tictoc)
library(skimr)
library(modelsummary)
library(kableExtra)
library(doParallel)
library(doRNG)
library(here)

# 自作関数を読み込む
source("04_Static_Game_MRI_Ch08/function_Berry1992.R")

# Boostrapをするかしないかのオプション
# Boostrapにかかる計算時間が比較的長いので、Bootstrapを行う場合はTRUE,
# すでに回して得られた結果をロードする場合にはFALSEを指定する。 
isBoostrap = FALSE

# 2.  データの準備----

# 病院データを読み込む
data_raw <- readr::read_csv(file = here("04_Static_Game_MRI_Ch08/data/data_hospital_Berry1992.csv"), 
                            locale = locale(encoding = "UTF-8"))
  
# ダミー変数を作成
data_cleaned <- data_raw %>% 
  # 大学（公立、国立、私立大学）ダミー
  dplyr::mutate(DaigakuDum = if_else(
    Management %in% c('公立大学法人','国（国立大学法人）','私立学校法人'),1, 0)) %>% 
  # 病床数0ダミー
  dplyr::mutate(ZeroBedDum = if_else(NumBeds == 0, 1, 0))

# 単位の調整と、対数変換した変数の作成
# - 以下のように単位を調整する。
# - NumBeds: 100病床
# - Population: 100万人
# - Menseki: 100km2
# - TaxableIncome: 100万円
# - また、NumBedsとPopulation、TaxableIncomeは対数変換した変数を作成する。
data_cleaned <- data_cleaned %>% 
  dplyr::mutate(
    NumBeds = NumBeds/100,
    LogNumBeds = log(NumBeds+0.01),
    Population = Population/10^6,
    Menseki = Menseki/100,
    TaxableIncome = TaxableIncome/1000,
    LogPop = log(Population),
    LogIncome = log(TaxableIncome))

# 3.  企業レベルの変数の記述統計----

# MRIを購買した病院と購買していない病院の記述統計
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
                                  '病床数', '病床0', '大学病院')) -> tbl_data

tbl_data %>% 
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

# CSVで保存
write.csv(tbl_data,
          file = here("04_Static_Game_MRI_Ch08/output/Tab8_6_descriptive_stats.csv"),
          row.names = FALSE,
          fileEncoding = "UTF-8")


# 各カテゴリーの病院のMRI保有割合

# 病院のカテゴリー毎に、総病院数、MRI所有病院数、その割合を確認
# 自作関数 `per_MRI_cal` を利用
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


# 4. Berry (1992)による推定 -----

## 推定のためのデータ加工 ----
# - Berry (1992) の推定を行うために、以下のような手順でデータを加工する。
# - Step 0: 欠損のあるデータを削除
# - Step 1: シミュレーションに必要な変数を作成、参入の順番に並び替え
# - Step 2: 潜在的参入病院数を絞ってサブサンプル
# - Step 3: シミュレーションを行いやすくするためにデータを縦に長く加工
# - Step 1, 3 については、それぞれの手順を関数として保存している。なぜなら、これらの手順は標準誤差を計算するためのブートストラップ法でも必要となるからである。
# - Step 2は本来不要だが、今回は推定にかかる時間を短縮するためにサンプルを絞っている。

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

# Step 1: シミュレーションに必要な変数を作成、参入の順番に並び替え
# 以降の結果が変わらないように、乱数を固定
set.seed(1)
# Step 1で作成した関数の実行
data_processed <- berry_process_data(data_processed_pre)
# 不要なオブジェクトを削除
rm(data_processed_pre)

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

# 各カテゴリーの病院のMRI保有割合（サブサンプル後のデータ）
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


# Step 3: データを加工する。
ns <- 100
data_expand <- berry_expand_data(data_processed, ns)
 
## 推定 ----

# 係数の名前
est_names <- 
  c('Kyukyu', 'Sien', 'Hyoka',
    'DepNeurology', 'DepNeurosurgery', 'LogNumBeds',
    'ZeroBedDum', 'DaigakuDum', 'Constant',
    'Menseki', 'LogPop', 'LogIncome',
    'delta', 'rho')

# 初期値
param_init <- 
  c(-0.587617,-1.298752,   -0.556557,
    -0.502602,   -1.058680,   -0.919110,  
    -2.832720,   -2.726932,    1.937292,
    0.009940,    0.206644,   -1.617755,
    0.550239,    0.132202)

# 以上の初期値の元で推定を行う。
# 自作関数 berry_obj が目的関数。
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

## Bootstrapによる標準誤差の計算 ----

# **留意点：** このブートストラップ計算には比較的時間がかかる点に留意されたい。
# 並列計算に用いるコア数やCPUによっても異なるが、
# 再現パッケージに含まれている`BootStrap_Berry_Result.rds`ファイルを用いて、
# 以下のコードを読み進めてもらえれば幸いである。

if (isBoostrap == TRUE){
 
  # 最初のチャンクでは、並列計算のための設定を行う。
  # ここで変数`cores`において、並列計算に使うCPUコア数を設定する。
  # なお、各ブートストラップサンプルにおける推定において乱数を利用するため、
  # 前と同様に乱数を固定する必要がある。
  # しかしながら、並列計算の場合には`set.seed()`では乱数を固定することができない。
  # ここではdoRNGパッケージにおける`registerDoRNG()`を利用する。
  
  ## # 並列計算のための設定を行う。
  cores <- 20 #並列計算に使うCPUコアの個数。
  registerDoParallel(cores)
  registerDoRNG(12345)

  # ブートストラップの回数
  iter <- 100
  # 自治体の数
  M <- data_processed %>%
    dplyr::select(CityCode) %>%
    dplyr::n_distinct()
  
  # 繰り返し推定を行う
  # 並列計算には`foreach`を使う。
  tic()
  bs_berry_result <-
    foreach::foreach (
      i = 1:iter,
      .export = ls(envir=parent.frame()),
        .packages = 'tidyverse'
      ) %dopar% {
        # リサンプルを行う
        bs_data <- data_processed %>%
          dplyr::distinct(CityCode) %>%
          dplyr::slice_sample(n = M, replace = TRUE) %>%
          dplyr::arrange(CityCode) %>%
          dplyr::mutate(id = 1:M) %>%
          dplyr::left_join(data_processed, by = 'CityCode') %>%
          # 推定ではCityCodeにより自治体を判別しているため名前を修正
          dplyr::rename(CityCodeOrigin = CityCode, CityCode = id)
        # リサンプルしたデータを推定のために加工
        bs_data_processed <-
          berry_process_data(bs_data)
        # 先程と同様に推定のためにデータを縦に広げる
        bs_data_expand <-
          berry_expand_data(
            bs_data_processed,
            ns)
        # 先程と同様に最適化により推定
        bs_berry_est <-
          optim(params <- param_init,
                berry_obj,
                method = "Nelder-Mead",
                df = bs_data_expand
                )$par
        # rhoの符号を正しくする
        bs_berry_est <- bs_berry_est %>%
          setNames(est_names)
        bs_berry_est[['rho']] <- abs(bs_berry_est[['rho']])
        return(bs_berry_est)
    }
  toc()

  # 結果を保存する。
  saveRDS(bs_berry_result, 
          file = here("04_Static_Game_MRI_Ch08/output/BootStrap_Berry_Result.rds"))
  
  
} else {
  
  bs_berry_result <- readRDS(here("04_Static_Game_MRI_Ch08/output/BootStrap_Berry_Result.rds"))
  
}

# 標準誤差を計算
berry_se <- bs_berry_result %>%
  map_dfr(~ unlist(as.data.frame(.))) %>%
  summarise(across(everything(), sd)) %>%
  as.numeric()

# 5. 推定結果 ----

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
                         'delta', 'rho', '目的関数の値')) -> tbl_est

tbl_est %>% 
  kableExtra::kable(format = "html", 
                    booktabs = TRUE,
                    caption = "Berry (1992) による推定結果", 
                    col.names = c('変数', '推定値', '標準誤差'),
                    digits = 3) %>%
  kable_styling()

tbl_est %>% 
  rename(Berry_est = BR_est,
         Berry_SE = BR_se) %>% 
  write.csv(file = here("04_Static_Game_MRI_Ch08/output/Tab8_7_Berry1992_Est.csv"))

## モデルによる予測の確認

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


# モデルの予測を計算
data_predicted <- 
  data_expand %>% 
  predict_model()

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
  scale_x_reverse() -> fig_pred

plot(fig_pred)

# 病院のタイプ別の観察された参入病院数と予測された参入病院数
# 自作関数 `actual_pred_cal` を利用
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
  coord_flip() + 
  theme_minimal() -> fig_pred_by_type

plot(fig_pred_by_type)

ggsave(fig_pred_by_type, filename = here("04_Static_Game_MRI_Ch08/output/Fig8_3_Fit.png"))


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


# 6. 反実仮想分析 ----

# 推定では企業の参入順序を仮定する必要がなかったが、
# 反実仮想分析においては、どの企業が参入するかを観察したいため、
# 企業の参入順序を仮定する必要がある。
# そこで、以下では３つの参入順序を考えて分析を行っている。
# なお書籍においては、「病床数」のケースのみをレポートしている。

## 参入順序: 病床数 ----

# 病床数が多い順に参入すると仮定する。タイはランダムとしている。

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

# Counterfactualの結果をレポート
tbl1 <- generate_cf_result(data_predicted, data_cf)

# print and save in csv
tbl1 %>%
  kableExtra::kable(format = "html",
                    booktabs = TRUE,
                    digits = 3) %>%
  kableExtra::kable_styling() 

write.csv(tbl1, file = here("04_Static_Game_MRI_Ch08/output/Tab8_8_CF_Simu.csv"))



## 参入順序: 既存企業 ----
# データ上で参入が観察された企業から行動すると仮定する。
# タイはランダムとしている。

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

tbl2 <- generate_cf_result(data_predicted_order2, data_cf_order2)

# print and save in csv
tbl2 %>%
  kableExtra::kable(format = "html",
                    booktabs = TRUE,
                    digits = 3) %>%
  kableExtra::kable_styling() 


## 参入順序: ランダム ----
# ランダムに参入すると仮定する。
set.seed(1234567)

random_entry_df <-
  data_expand %>% 
  dplyr::distinct(CityCode, EntryOrderId) %>%
  mutate(RandomEntryOrder = runif(n(),0,1))

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

# Table
tbl3 <- generate_cf_result(data_predicted_order3, data_cf_order3)

# print and save in csv
tbl3 %>%
  kableExtra::kable(format = "html",
                    booktabs = TRUE,
                    digits = 3) %>%
  kableExtra::kable_styling() 


