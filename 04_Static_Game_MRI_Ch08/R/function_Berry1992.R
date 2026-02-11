
per_MRI_cal <- function(df){
  #表作成のために利用する関数
  df %>% 
    dplyr::summarise(Total = n(),
                     MRIHos = sum(MRIOwnDum),
                     PerMRI = round(MRIHos / Total * 100, 2))
}


actual_pred_cal <- function(df){
  #表作成のために利用する関数
  df %>% 
    dplyr::summarise(Actual = sum(MRIOwnDum),
                     Predict = sum(EntryPred))
}


berry_process_data <- function(df){
  
  # 推定のためにデータを加工する関数
  
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


berry_expand_data <- function(df, ns){
  
  # シミュレーションを行いやすくするためにデータを縦に長く加工
  # データの行数は \Sigma_{m = 1}^{M} (Im * Im) * nsとなる
  
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



berry_obj <- function(params, df){
  
  # データとパラメターを所与としたBerry (1992) の目的関数の
  # 値を計算する関数を作成
  
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
  
  # 細かいコメント: 
  # 書籍では、均衡参入企業数の計算にあたって、意思決定の順番の仮定を
  # 用いて数値計算する方法を紹介している。
  # しかしながら、以下のコードでは、意思決定の順番を用いない方法で計算している。
  # 1: 均衡参入数の候補N0 を所与として、各企業が参入するかを決定する。
  # その結果の参入数をN1とする。
  # 2: N0 > N1ならば、N0は均衡における参入数ではない。
  # 3; N0 <= N1が成立するようなN0の中で、最も大きいN0を「均衡参入数」とする。
  # なお、Berry (1992) のモデルでは、
  # 「意思決定の順にかかわらず、参入企業の数は一意」に決まるため、
  # 仮に参入意思決定の順番を用いて均衡参入数を計算しても同じ値が得られる。
  
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


predict_model <- function(df){
  
  # モデルによる参入予測確率を含めたdataframeを作成
  
  # 細かい補足: 参入企業数の計算方法は、`berry_obj`と同じである。
  # そのうえで、EntryOrderIdという、参入意思決定の順番の仮定に基づいて、
  # どの病院が参入(MRIを導入)するかを決定する。
   
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



generate_cf_result <- function(data_predicted, data_cf){
  
  # 病院のカテゴリー毎に、総病院数(Total)、
  # 元のモデルにより予測されたMRI所有病院数(Fit_Pred)、
  # 反実仮想分析により予測されたMRI所有病院数(CF_Pred)、
  # MRI所有病院数の変化率(PerChange)、を確認
  
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
    dplyr::mutate(PerChange = round((CF_Pred/Fit_Pred - 1) * 100, 2)) -> tbl
  
  return(tbl) 
}

per_decline_cal <- function(df){
  
  # generate_cf_resultで用いられる、表を作成するための関数
  
  df %>% 
    dplyr::summarise(
      Total = n(),
      Actual = sum(MRIOwnDum),
      Fit_Pred = sum(Fit_Pred),
      CF_Pred = sum(EntryPred),
      CF_Pred_am_actual = sum(MRIOwnDum * EntryPred)
    )
}


