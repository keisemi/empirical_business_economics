
library("tidyverse")
library("knitr")
library("magrittr")
library("mlogit") # ロジットモデル推定のためのパッケージ
library("stargazer") # 推定結果の表作成のためのパッケージ
library("optimx") # 最適化のためのパッケージ
library("here") 

source(here("01_Discrete_Choice_Ch02/function.R"))

KT_2024 <- read_csv(here("01_Discrete_Choice_Ch02/data/KinokoTakenokoSurvey_raw.csv"))

# 列名変更
new_col_name <- c("ID","experience", "Q1", "Q2", "Q3", "Q4", "Q5", "age", "gender", "region", "familyhouse")

KT_2024 <- KT_2024 %>%
  mutate(ID = row_number()) %>% 
  select(ID, c(8:17)) %>% # 必要な変数のみ選択
  set_colnames(new_col_name) %>% 
  filter(experience != "4 : 食べたことがない" & gender != "3 : 回答したくない") %>% # 扱いにくいデータを除外
  drop_na() # NA(未回答？)の行を除外

# 記述統計 ----

N <- length(KT_2024$ID)

# 回答数の比率を表にするためのデータ整形
datafig <- KT_2024 %>%
  select(ID, Q1, Q2, Q3, Q4, Q5)  %>% # IDと回答のみ抽出
  gather(key = Q, value  = choice, Q1, Q2, Q3, Q4, Q5)

# 対応する選択肢に置き換える
datafig <- datafig %>% 
  mutate(Q = ifelse(Q == "Q1", "Q1: (200円, 200円)" ,Q),  
         Q = ifelse(Q == "Q2", "Q2: (180円, 200円)" ,Q),  
         Q = ifelse(Q == "Q3", "Q3: (200円, 170円)" ,Q),  
         Q = ifelse(Q == "Q4", "Q4: (220円, 200円)" ,Q),  
         Q = ifelse(Q == "Q5", "Q5: (190円, 210円)" ,Q))

# 回答数の比率の表を作成
tab_choice <- datafig %>% 
  group_by(Q, choice) %>% 
  tally() %>% 
  mutate(n = n/N) %>% 
  pivot_wider( id_cols = "Q", names_from = "choice", values_from = "n")

write.table(tab_choice,
            file = here("01_Discrete_Choice_Ch02/output/tab2_1_tab_choice.txt"),
            quote = FALSE,  
            col.names = NA) 

# 回答者の属性の比率を表にするためのデータ整形
data_atr <- KT_2024 %>%
  mutate(gender =　case_when(gender == "1 : 男性" ~ "男性", 
                            gender == "2 : 女性" ~ "女性" ),
         exp    =  case_when(experience == "1 : 過去半年以内" ~ "過去半年以内", 
                             experience == "2 : 過去半年から１年以内" ~ "過去半年から１年以内", 
                             experience == "3 : １年以上前" ~ "１年以上前"),
         adult = if_else(age>=20, "成人", "未成年"),
         region =  case_when(region == "1 : 北海道地方" ~ "関東", 
                             region == "2 : 東北地方" ~ "関東", 
                             region == "3 : 関東地方" ~ "関東" ,
                             region == "4 : 中部地方" ~ "関東",
                             region == "5 : 近畿地方" ~ "関西",
                             region == "6 : 中国地方" ~ "関西",
                             region == "7 : 四国地方" ~ "関西",
                             region == "8 : 九州地方(沖縄含む)" ~ "関西",
                             region == "9 : 海外" ~ "海外"),
         familyhouse = if_else(familyhouse==1,"実家暮らし","実家暮らしでない"))

# 出身地方の割合の表を作成
region_fig <- data_atr %>%
  group_by(region) %>%
  tally() %>%
  mutate(割合 = n/N) %>% 
  select (-n) %>%
  mutate(変数 = c("出身地方","","")) %>%
  rename(属性 = region) %>%
  select(変数, everything())


# 食事経験の割合の表を作成
exp_fig <- data_atr %>%
   group_by(exp) %>%
   tally() %>%
   mutate(割合 = n/N) %>% 
   pivot_wider(names_from = "exp", values_from = "n") %>% 
   mutate("食事経験" = " ") %>%
   select("食事経験", everything())


# 性別の割合の表を作成
gender_fig <- data_atr %>%
  group_by(gender) %>%
  tally() %>%
  mutate(割合 = n/N) %>% 
  select(-n) %>%
  mutate(変数 = c("性別","")) %>%
  rename(属性 = gender) %>%
  select(変数, everything())


# 実家暮らしかどうかの割合の表を作成
fam_fig <- data_atr %>%
  group_by(familyhouse) %>%
  tally() %>%
  mutate(割合 = n/N) %>% 
  select(-n) %>%
  mutate(変数 = c("実家暮らしかどうか","")) %>%
  rename(属性 = familyhouse) %>%
  select(変数, everything(),)


# 成人かどうかの割合の表を作成
adult_fig <- data_atr %>%
  group_by(adult) %>%
  tally() %>%
  mutate(割合 = n/N) %>% 
  select(-n) %>%
  mutate(変数 = c("成人かどうか","")) %>%
  rename(属性 = adult) %>%
  select(変数, everything())


tab_atr <- rbind(gender_fig,adult_fig, fam_fig, region_fig)

print(tab_atr)


# データ整形 ----

# 分析用のデータ整形
data_for_estimation <- KT_2024 %>% 
  gather(key = "occasion", value = choice, starts_with("Q"))

# 提示される価格の組み合わせ
pricedata <- data.frame(
  occasion = c("Q1","Q2", "Q3", "Q4", "Q5" ),
  price_0 = numeric(5), 
  price_1 = c(200, 180, 200, 220, 190), 
  price_2 = c(200, 200, 170, 200, 210)) %>% 
  as_tibble()

data_for_estimation <- data_for_estimation %>% 
  left_join(pricedata) %>%
  mutate(Kinoko_0 = 0, 
         Kinoko_1 = 1, 
         Kinoko_2 = 0, 
         Takenoko_0 = 0, 
         Takenoko_1 = 0, 
         Takenoko_2 = 1) %>%
  arrange(ID, occasion)


data_for_estimation <- data_for_estimation %>% 
  mutate(choice =  case_when(choice == "1 : きのこの山を買う" ~ 1, 
                              choice == "2 : たけのこの里を買う" ~ 2, 
                              choice == "3 : どちらも買わない" ~ 0 ),
         gender =　case_when(gender == "1 : 男性" ~ "0", 
                             gender == "2 : 女性" ~ "female" ),
         exp    =  case_when(experience == "1 : 過去半年以内" ~ "halfyear", 
                              experience == "2 : 過去半年から１年以内" ~ "half_to_1year", 
                              experience == "3 : １年以上前" ~ "0" ),
         adult = if_else(age>=20, 1, 0),
         region =  case_when(region == "1 : 北海道地方" ~ "_Kanto", 
                              region == "2 : 東北地方" ~ "_Kanto", 
                              region == "3 : 関東地方" ~ "_Kanto" ,
                              region == "4 : 中部地方" ~ "_Kanto",
                              region == "5 : 近畿地方" ~ "Kansai",
                              region == "6 : 中国地方" ~ "Kansai",
                              region == "7 : 四国地方" ~ "Kansai",
                              region == "8 : 九州地方(沖縄含む)" ~ "Kansai",
                              region == "9 : 海外" ~ "Oversea",
          )
  )

data_for_estimation$choiceid <- 1:nrow(data_for_estimation)

# 分析用のデータフレーム作成
datalogit <- dfidx(data = as.data.frame(data_for_estimation),
                   choice = "choice", 
                   varying = c(9:17), 
                   sep = "_",
                   idx = list(c("choiceid", "ID")), 
                   idnames = c("chid", "alt"),
                   opposite = c("price")) 


# 回答者属性を含まない回帰式
no_atr <- formula(choice ~ price + Kinoko + Takenoko | 0)

# 回答者属性を含む回帰式
with_atr <- formula(choice ~ price + Kinoko + Takenoko + 
                      price:gender + Takenoko:gender + Kinoko:gender + 
                      price:familyhouse + Takenoko:familyhouse + Kinoko:familyhouse + 
                      price:adult + Takenoko:adult + Kinoko:adult + 
                      price:region + Takenoko:region + Kinoko:region | 0 | 0)


# パッケージを用いた推定 ----

# 消費者属性を入れない多項ロジット
ml <- mlogit(formula = no_atr,
             data = datalogit,
             reflevel = '0')

# 消費者属性を入れないランダム係数ロジット
rcdc <- mlogit(formula = no_atr, 
               data = datalogit, 
               panel = TRUE, 
               rpar = c(price = "n", Kinoko = "n", Takenoko = "n") , 
               R = 2000, 
               correlation = FALSE)

# 消費者属性を入れた多項ロジット
ml_consumer_atr <- mlogit(formula = with_atr, 
                          data = datalogit,
                          reflevel = '0')

# 消費者属性を入れたランダム係数ロジット
rcdc_consumer_atr <- mlogit(formula = with_atr, 
                            data = datalogit, 
                            panel = TRUE, 
                            rpar = c(price = "n", Kinoko = "n", Takenoko = "n") , 
                            R = 2000, 
                            correlation = FALSE) 

# 推定結果
stargazer::stargazer(ml,
                     rcdc,
                     ml_consumer_atr,
                     rcdc_consumer_atr,
                     column.labels=c("mlogit","rcdc","mlogit","rcdc"),
                     single.row=TRUE,
                     type="text",
                     align = TRUE,
                     out = here("01_Discrete_Choice_Ch02/output/result.txt"))


# 消費者属性を考慮した選好パラメタおよびWTPの分布の計算 ----

# 全てダミー変数にする
data_for_param <- KT_2024 %>%
  mutate(gender =　case_when(gender == "1 : 男性" ~ 0,
                            gender == "2 : 女性" ~ 1 ),
         adult = if_else(age>=20, 1, 0),
         kansai =  case_when(region == "1 : 北海道地方" ~ 0,
                             region == "2 : 東北地方" ~ 0,
                             region == "3 : 関東地方" ~ 0 ,
                             region == "4 : 中部地方" ~ 0,
                             region == "5 : 近畿地方" ~ 1,
                             region == "6 : 中国地方" ~ 1,
                             region == "7 : 四国地方" ~ 1,
                             region == "8 : 九州地方(沖縄含む)" ~ 1,
                             region == "9 : 海外" ~ 0),
         oversea =  case_when(region == "1 : 北海道地方" ~ 0,
                              region == "2 : 東北地方" ~ 0,
                              region == "3 : 関東地方" ~ 0 ,
                              region == "4 : 中部地方" ~ 0,
                              region == "5 : 近畿地方" ~ 0,
                              region == "6 : 中国地方" ~ 0,
                              region == "7 : 四国地方" ~ 0,
                              region == "8 : 九州地方(沖縄含む)" ~ 0,
                              region == "9 : 海外" ~ 1))

dist_Kinoko <- rpar(rcdc_consumer_atr,'Kinoko')
dist_Takenoko <- rpar(rcdc_consumer_atr,'Takenoko')
dist_price <- rpar(rcdc_consumer_atr,'price')

# 消費者属性を考慮するため、各個人の属性のパラメタを計算し、その値と選好パラメタの平均(消費者属性を考慮していない)との和を取る。
param_Takenoko <- dist_Takenoko$mean + 
  data_for_param$gender * rcdc_consumer_atr$coefficients[5] +
  data_for_param$familyhouse * rcdc_consumer_atr$coefficients[8] +
  data_for_param$adult * rcdc_consumer_atr$coefficients[11] +
  data_for_param$kansai * rcdc_consumer_atr$coefficients[15] +
  data_for_param$oversea * rcdc_consumer_atr$coefficients[16]

param_Kinoko <- dist_Kinoko$mean +
  data_for_param$gender * rcdc_consumer_atr$coefficients[6] +
  data_for_param$familyhouse *rcdc_consumer_atr$coefficients[9] +
  data_for_param$adult * rcdc_consumer_atr$coefficients[12] +
  data_for_param$kansai * rcdc_consumer_atr$coefficients[17] +
  data_for_param$oversea * rcdc_consumer_atr$coefficients[18]

param_price <- dist_price$mean +
  data_for_param$gender * rcdc_consumer_atr$coefficients[4] +
  data_for_param$familyhouse * rcdc_consumer_atr$coefficients[7] +
  data_for_param$adult * rcdc_consumer_atr$coefficients[10] +
  data_for_param$kansai * rcdc_consumer_atr$coefficients[13] +
  data_for_param$oversea * rcdc_consumer_atr$coefficients[14]

# 236人の各個人について1000回ドローして、計236000個の個人選好パラメタのドローを得る

# 格納用
alpha_Kinoko_vec <- numeric(236000)
alpha_Takenoko_vec <- numeric(236000)
beta_vec <- numeric(236000)

# ドロー
for(i in 1:length(param_Kinoko)){
  set.seed(100+i)
  
  alpha_Kinoko_draw = dist_Kinoko$sigma*rnorm(1000, mean = 0, sd = 1)+param_Kinoko[i]
  alpha_Kinoko_vec[(1000*(i-1)+1) : (1000*i)] = alpha_Kinoko_draw
  
  alpha_Takenoko_draw = dist_Takenoko$sigma*rnorm(1000, mean = 0, sd = 1)+param_Takenoko[i]
  alpha_Takenoko_vec[(1000*(i-1)+1) : (1000*i)] = alpha_Takenoko_draw
  
  beta_price_draw = dist_price$sigma*rnorm(1000, mean = 0, sd = 1)+param_price[i]
  beta_vec[(1000*(i-1)+1) : (1000*i)] = beta_price_draw
}


data_for_distribution <-
  tibble(alpha_Kinoko_vec, 
         alpha_Takenoko_vec,  
         beta_vec) 

# きのことたけのこの係数の差の分布
KT_dist_hist <- ggplot(data = data_for_distribution, aes(x = alpha_Kinoko_vec - alpha_Takenoko_vec)) +
  geom_histogram(binwidth = 0.75) +
  ylab("frequency") +
  xlab("Kinoko - Takenoko") +
  ggtitle("Preference Kinoko over Takenoko")

plot(KT_dist_hist)

ggsave(filename = here("01_Discrete_Choice_Ch02/output/fig2_1_LHS_hist_preference_KT.pdf"), plot = KT_dist_hist, width = 6, height = 6)

KT_dist_density <- ggplot(data = data_for_distribution, aes(x = alpha_Kinoko_vec - alpha_Takenoko_vec)) +
  geom_density() +
  ylab("density") + 
  xlab("Kinoko - Takenoko") +
  ggtitle("Preference Kinoko over Takenoko")

plot(KT_dist_density)

ggsave(filename = here("01_Discrete_Choice_Ch02/output/fig_dist_preference_KT.pdf"), plot = KT_dist_density, width = 6, height = 6)

# 価格係数の分布
p_coef_dist_hist <- ggplot(data = data_for_distribution, aes(beta_vec)) +
  geom_histogram(binwidth = 0.0075)+
  ylab("frequency") + 
  xlab("beta") +
  ggtitle("Beta (price coefficient)")

plot(p_coef_dist_hist)

ggsave(filename = here("01_Discrete_Choice_Ch02/output/fig2_1_RHS_hist_price_coef.pdf"), plot = p_coef_dist_hist, width = 6, height = 6)

p_coef_dist_density <- ggplot(data = data_for_distribution, aes(beta_vec)) +
  geom_density()+
  ylab("density") + 
  xlab("beta") + 
  ggtitle("Beta (price coefficient)")

plot(p_coef_dist_density)

ggsave(filename = here("01_Discrete_Choice_Ch02/output/fig_dist_price_coef.pdf"), plot = p_coef_dist_density, width = 6, height = 6)


# WTPを求める

# WTPの計算
data_for_distribution <- data_for_distribution %>%
  mutate(brand_Kinoko = alpha_Kinoko_vec / beta_vec, 
         brand_Takenoko = alpha_Takenoko_vec / beta_vec,
         love_Kinoko = ifelse(alpha_Kinoko_vec > alpha_Takenoko_vec, 1, 0))


dt <- data_for_distribution %>% 
  select(brand_Kinoko, brand_Takenoko) %>% 
  rename(WTP_Kinoko = brand_Kinoko, 
         WTP_Takenoko = brand_Takenoko) %>% 
  pivot_longer(cols = everything(), names_to = "Name", values_to = "WTP")

# WTPの分布
fig_brand <- dt %>% 
  ggplot(aes(x = WTP, colour = Name)) +
  theme_bw() + 
  geom_density() +
  xlim(c(0,500))

plot(fig_brand)

ggsave(filename = here("01_Discrete_Choice_Ch02/output/fig2_2_dist_wtp.pdf"), plot = fig_brand, width = 12, height = 6)


# 需要曲線と収入曲線を求める

# たけのこの価格を200円で固定する
price_Takenoko = 200 

# きのこの価格を100円から250円まで5円ずつ動かす
price_vec = seq(from = 100, to = 250, by = 5)
kinoko_vec = numeric(length(price_vec))

# きのこの価格に対応する各オプションの需要を計算する
for (i in 1:length(price_vec)){
  result <- f_demand(alpha_Kinoko_vec,
                     alpha_Takenoko_vec, 
                     beta_vec,
                     price_vec[i],
                     price_Takenoko = 200)
  kinoko_vec[i] <- result[1]
}

# データフレームに結果を保存
data_demand_kinoko = tibble(price = price_vec, 
                            demand = kinoko_vec,
                            revenue = price_vec*kinoko_vec)


# 需要曲線の描画
fig_demand <- ggplot(data = data_demand_kinoko, aes(x = demand, y = price)) +
  geom_line() + 
  xlab("Demand of Kinoko") + 
  ylab("Price of Kinoko in JPY") + 
  ggtitle("Demand Curve of Kinoko when Takenoko's price = JPY 200")

plot(fig_demand)

ggsave(filename = here("01_Discrete_Choice_Ch02/output/demand_Kinoko.pdf"), plot = fig_demand, width = 8, height = 6)

# x軸をシェアにした場合の需要曲線
fig_share <- ggplot(data = data_demand_kinoko, aes(x = demand/236000, y = price)) +
  geom_line() + 
  xlab("Share of Kinoko") + 
  ylab("Price of Kinoko in JPY") + 
  ggtitle("Demand Curve(Share) of Kinoko when Takenoko's price = JPY 200")

plot(fig_share)

ggsave(filename = here("01_Discrete_Choice_Ch02/output/share_kinoko.pdf"), plot = fig_share, width = 8, height = 6)

# 収入曲線の描画
fig_rev <- ggplot(data = data_demand_kinoko, aes(x = price, y = revenue/10000)) +
  geom_line() + 
  xlab("Price of Kinoko in JPY") + 
  ylab("Revenue of Kinoko in 10000 JPY") + 
  ggtitle("Revenue Curve of Kinoko when Takenoko's price = JPY 200")

plot(fig_rev)

ggsave(filename = here("01_Discrete_Choice_Ch02/output/Revenue_Kinoko.pdf"), plot = fig_rev, width = 8, height = 6)


# スクラッチで推定 ----

## 消費者属性を入れない多項ロジット ----

# ここでf_likelihood_logitを定義(function.Rに移行)

# 初期パラメタ
ini_ml <- c(5, 5, -0.01)

# 最適化
result <- optimx(par = ini_ml, 
                 fn = f_likelihood_logit, 
                 method="Nelder-Mead", 
                 data = data_for_estimation, 
                 control = list(fnscale=-1), 
                 hessian = TRUE)

# 推定値の取得
est_vec <- as.numeric(result[1,1:3])
est_vec <- est_vec[c(3, 1, 2)]
# パッケージの推定値の取
est_vec_package <- ml$coefficients
names(est_vec) <- c(attr(ml$coefficients, "names"))

# 比較
print(rbind(est_vec, est_vec_package))

# 対数尤度の取得
loglik_ml <- result[4]$value
# パッケージの対数尤度の取得
loglik_ml_package <- as.numeric(ml$logLik)
# 比較
print(rbind(loglik_ml, loglik_ml_package))

# ヘシアンを計算し、標準誤差を確認する
Hessian <- gHgen(par = as.numeric(result[1,1:3]), 
                 fn = f_likelihood_logit, 
                 data = data_for_estimation)
se_vec <- round(sqrt(diag(solve(-Hessian$Hn))), 3)
se_vec <- se_vec[c(3, 1, 2)]
print(rbind(est_vec, se_vec))


## 消費者属性を入れた多項ロジット ----

# ダミー変数に直す
data_for_estimation_with_atr <- data_for_estimation %>%
  mutate(gender = case_when(gender == "female"~ 1,
                            gender == "0" ~ 0),
         kansai =  case_when(region == "_Kanto" ~ 0,
                             region == "Oversea" ~ 0,
                             region == "Kansai" ~ 1),
         oversea =  case_when(region == "_Kanto" ~ 0,
                              region == "Oversea" ~ 1,
                              region == "Kansai" ~ 0))


# 初期値
ini_ml_with_atr <- c(0, 10, 10, 0, 1, 1, 0, 0, 0, 0, -1, -1, 0, 0, 1, 1, 1, 1)

# 最適化
# Nelder-Mead法では初期値にかなり依存して推定されてしまうためBFGSを用いる
result_with_atr <- optimx(par = ini_ml_with_atr, 
                          fn = f_likelihood_logit_with_atr, 
                          method="BFGS", 
                          data = data_for_estimation_with_atr, 
                          control = list(fnscale=-1), 
                          hessian = TRUE)


# 推定値の取得
est_vec_ml_atr <- as.numeric(result_with_atr[1,1:18])
# パッケージの推定値の取得
est_vec_package_ml_atr <- ml_consumer_atr$coefficients[1:18]
# 比較
print(rbind(est_vec_ml_atr, est_vec_package_ml_atr))

# 対数尤度の取得
loglik_ml_atr <- result_with_atr[19]$value
# パッケージの対数尤度の取得
loglik_ml_atr_package <-as.numeric(ml_consumer_atr$logLik)
# 比較
print(rbind(loglik_ml_atr, loglik_ml_atr_package))


# ヘシアンを計算し、標準誤差を確認する
Hessian <- gHgen(par = as.numeric(result_with_atr[1,1:18]), 
                 fn = f_likelihood_logit_with_atr, 
                 data = data_for_estimation_with_atr)
se_vec_ml_atr <- round(sqrt(diag(solve(-Hessian$Hn))), 3)

names(est_vec_ml_atr) = c(attr(ml_consumer_atr$coefficients, "names"))
print(rbind(est_vec_ml_atr, se_vec_ml_atr))

## 消費者属性を入れないランダム係数ロジット ----

# 事前にドローを用意する

# ドロー回数
times_draw =1000
# 乱数の格納用
tau<- matrix(nrow = times_draw*236, ncol = 3)

for(i in 1:236){
  set.seed(500+i)
  # 回答者iがドローする乱数
  tau[(times_draw*(i-1)+1):(times_draw*i), 1] <- rnorm(times_draw, mean = 0, sd = 1)
  tau[(times_draw*(i-1)+1):(times_draw*i), 2] <- rnorm(times_draw, mean = 0, sd = 1)
  tau[(times_draw*(i-1)+1):(times_draw*i), 3] <- rnorm(times_draw, mean = 0, sd = 1)}


# 初期値
ini_rcdc <- c(0.6,20,20,	0.1,1,1)
# 最大化
result_rcdc <- optimx(par = ini_rcdc, 
                      fn = f_likelihood_rcdclogit, 
                      method="BFGS", 
                      data = data_for_estimation,
                      control = list(fnscale=-1), 
                      hessian = TRUE)

# 推定値の取得
est_vec_rcdc <- as.numeric(result_rcdc[1,1:6])
# パッケージの推定値の取得
est_vec_package_rcdc_package <- rcdc$coefficients[1:6]
# 比較
print(rbind(est_vec_rcdc, est_vec_package_rcdc_package))

# 対数尤度の取得
loglik_rcdc <- result_rcdc[7]$value
# パッケージの対数尤度の取得
loglik_rcdc_package <- as.numeric(rcdc$logLik)
# 比較
print(rbind(loglik_rcdc, loglik_rcdc_package))

# ヘシアンを計算し、標準誤差を確認する
Hessian <- gHgen(par = as.numeric(result_rcdc[1,1:6]), 
                 fn = f_likelihood_rcdclogit, 
                 data = data_for_estimation)
se_vec_rcdc <- round(sqrt(diag(solve(-Hessian$Hn))), 3)

names(est_vec_rcdc) = c(attr(rcdc$coefficients, "names"))
print(rbind(est_vec_rcdc, se_vec_rcdc))


## 消費者属性を入れたランダム係数ロジット ----

# ドロー回数
times_draw =1000
# 乱数の格納用
tau<- matrix(nrow = times_draw*236, ncol = 3)

for(i in 1:236){
  set.seed(500+i)
  # 回答者iがドローする乱数
  tau[(times_draw*(i-1)+1):(times_draw*i), 1] <- rnorm(times_draw, mean = 0, sd = 1)
  tau[(times_draw*(i-1)+1):(times_draw*i), 2] <- rnorm(times_draw, mean = 0, sd = 1)
  tau[(times_draw*(i-1)+1):(times_draw*i), 3] <- rnorm(times_draw, mean = 0, sd = 1)}


# 初期値
ini_rcdc_with_atr <- c(0.6, 20, 20, 0.1, 1, 1, 0, 3, 3, 0, 0, 0, 0, -3, -3, 0, 0, 1, 1, 1, 1 )

# 最大化
result_rcdc_with_atr <- optimx(par = ini_rcdc_with_atr, 
                               fn = f_likelihood_rcdclogit_with_atr, 
                               method="BFGS", 
                               data = data_for_estimation_with_atr, 
                               control = list(fnscale=-1), 
                               hessian = TRUE)

# 推定値の取得
est_vec_rcdc_with_atr <- as.numeric(result_rcdc_with_atr[1,1:21])
# パッケージの推定値の取得
est_vec_package_rcdc_with_atr <- rcdc_consumer_atr$coefficients[c(1:3,19:21,4:18)]
names(est_vec_rcdc_with_atr) <- c(attr(rcdc_consumer_atr$coefficients, "names")[c(1:3,19:21,4:18)])
# 比較
print(rbind(est_vec_rcdc_with_atr, est_vec_package_rcdc_with_atr))

# 対数尤度の取得
loglik_rcdc_atr <- result_rcdc_with_atr[22]$value
# パッケージの対数尤度の取得
loglik_rcdc_atr_package <- as.numeric(rcdc_consumer_atr$logLik)
# 比較
print(rbind(loglik_rcdc_atr, loglik_rcdc_atr_package))

# ヘシアンを計算し、標準誤差を確認する
Hessian <- gHgen(par = as.numeric(result_rcdc_with_atr[1,1:21]), 
                 fn = f_likelihood_rcdclogit_with_atr, 
                 data = data_for_estimation_with_atr)
se_vec_with_atr <- round(sqrt(diag(solve(-Hessian$Hn))), 3)

print(rbind(est_vec_rcdc_with_atr, se_vec_with_atr))

