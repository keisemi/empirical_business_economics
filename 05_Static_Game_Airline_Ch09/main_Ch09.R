
library(tidyverse)
library(here)
library(ROCR)
library(splines)
library(nleqslv)
library(showtext)
library(knitr)

source(here("05_Static_Game_Airline_Ch09/function.R"))

# データの読み込み ----
y = read.csv(here("05_Static_Game_Airline_Ch09/data/y.csv"))
x_ANA = read.csv(here("05_Static_Game_Airline_Ch09/data/x_ANA.csv"))
x_JAL = read.csv(here("05_Static_Game_Airline_Ch09/data/x_JAL.csv"))
id_port_dyad_long = read.csv(here("05_Static_Game_Airline_Ch09/data/id_port_dyad_long.csv"))

# 二つのデータが同じであることを確認
identical(select(x_ANA, -Flight), select(x_JAL, -Flight)) # TrueであればOK

# マージ
df = y %>%
  mutate(id = 1:n()) %>% 
  bind_cols(x_ANA %>% 
              rename(Flight_ANA = Flight) 
  ) %>% 
  bind_cols(x_JAL %>%
              rename(Flight_JAL = Flight) %>% 
              select(Flight_JAL)) %>% 
  relocate(id)


# 記述統計 ----
df %>% 
  select(y_ANA, y_JAL, Distance:Flight_JAL) %>% 
  summarise_all(mean) %>% 
  bind_rows(
    df %>% 
      select(y_ANA, y_JAL, Distance:Flight_JAL) %>% 
      summarise_all(sd)
  ) %>% 
  `rownames<-`( c("Mean", "SD")) %>% 
  t() %>% 
  as.data.frame() %>% 
  `rownames<-`(c("y_ANA", "y_JAL", "距離", "人口", "人口2乗", "新幹線", "フライトANA", "フライトJAL")) %>%
  knitr::kable(format = "pipe", digits = 3) %>%
  writeLines(here("05_Static_Game_Airline_Ch09/output/tab9_1_descriptive_statistics.txt"))


# 推定 ----

## ステップ1 ----

### 1. 多項式 ----

#ANA
r_sq_vec = tibble()

for(i in 1:5){
  model = df %>% 
    lm(y_ANA ~ Train*poly(Distance, Population, Flight_ANA, Flight_JAL, degree=i, raw=F), data=.) %>% 
    summary()
  
  r_sq_vec = r_sq_vec %>% 
    bind_rows(tibble(degree = i, Rsq = model$adj.r.squared))
}


#JAL
r_sq_vec2 = tibble()

for(i in 1:5){
  model = df %>% 
    lm(y_JAL ~ Train*poly(Distance, Population, Flight_ANA, Flight_JAL, degree=i, raw=F), data=.) %>% 
    summary()
  
  r_sq_vec2 = r_sq_vec2 %>% 
    bind_rows(tibble(degree = i, Rsq = model$adj.r.squared))
}

# スプライン次数（degree）ごとの決定係数を可視化
r_sq_plot <- r_sq_vec %>% 
  mutate(comp = "ANA") %>% 
  bind_rows(r_sq_vec2 %>% mutate(comp = "JAL")) 

p <- ggplot(r_sq_plot, aes(degree, Rsq, col = comp)) +
  geom_line() +
  geom_point()

ggsave(filename = here("05_Static_Game_Airline_Ch09/output/fig_Rsquared.pdf"), plot=p)



# 交差検証
set.seed(825)

K = 10

df = df %>% 
  mutate(k_id = ceiling(runif(nrow(df))*K))

df %>% 
  count(k_id)

# 予測パフォーマンス

degrees = 1:4
AUCs1 = matrix(0,nrow = K,ncol = length(degrees))
AUCs2 = AUCs1
MAE1 = AUCs1
MAE2 = AUCs1
AIC1 = AUCs1
AIC2 = AUCs1
BIC1 = AUCs1
BIC2 = AUCs1

for(i in 1:length(degrees)){
  for(k in 1:K){
    model1 = df %>% 
      filter(k_id != k) %>% 
      lm(y_ANA ~ Train*poly(Distance, Population, Flight_ANA, Flight_JAL, degree=i), data=.)
    
    model2 = df %>% 
      filter(k_id != k) %>% 
      glm(y_ANA ~ Train*poly(Distance, Population, Flight_ANA, Flight_JAL, degree=i),
          family = binomial, data = .)
    
    s1 = df %>% 
      filter(k_id == k) %>% 
      mutate(pred = predict(model1, newdata = ., type = "response"),
             pred = pmax(pmin(pred,1),0)
      ) %>% 
      pull(pred)
    
    s2 = df %>% 
      filter(k_id == k) %>% 
      predict(model2, newdata = ., type = "response")
    
    if(i < 4){
      AUCs1[k,i] = calculate_AUC(s1, df %>% filter(k_id == k) %>% pull(y_ANA))
      AUCs2[k,i] = calculate_AUC(s2, df %>% filter(k_id == k) %>% pull(y_ANA))
    }
    
    
    MAE1[k,i] = mean(abs(s1 - (df %>% filter(k_id == k) %>% pull(y_ANA))))
    MAE2[k,i] = mean(abs(s2 - (df %>% filter(k_id == k) %>% pull(y_ANA))))
    
    AIC1[k,i] = AIC(model1)
    AIC2[k,i] = AIC(model2)
    
    BIC1[k,i] = BIC(model1)
    BIC2[k,i] = BIC(model2)
  }
}


### 2. B-spline ----

vec = c("Distance", "Population", "Flight_ANA", "Flight_JAL")

degrees = 1:3
AUCs1_bs = matrix(0, nrow = K, ncol = length(degrees))
AUCs2_bs = AUCs1_bs
MAE1_bs = AUCs1_bs
MAE2_bs = AUCs1_bs
AIC1_bs = AUCs1_bs
AIC2_bs = AUCs1_bs
BIC1_bs = AUCs1_bs
BIC2_bs = AUCs1_bs


for(i in 1:length(degrees)){
  for(k in 1:K){
    if(i == 1){
      var_i = vec
    } 
    
    model1 = df %>% 
      filter(k_id != k) %>% 
      lm(generate_bspline_formula(var_i, "ANA", degree = i), data=.)
    
    model2 = df %>% 
      filter(k_id != k) %>% 
      glm(generate_bspline_formula(var_i, "ANA", degree = i), family = binomial, data = .)
    
    s1 = df %>% 
      filter(k_id == k) %>% 
      mutate(pred = predict(model1, newdata = ., type = "response"),
             pred = pmax(pmin(pred, 1), 0)
      ) %>% 
      pull(pred)%>% 
      as.vector()
    
    s2 = df %>% 
      filter(k_id == k) %>% 
      predict(model2, newdata = ., type = "response") %>% 
      as.vector()
    
    AUCs1_bs[k,i] = calculate_AUC(s1, df %>% filter(k_id == k) %>% pull(y_ANA))
    AUCs2_bs[k,i] = calculate_AUC(s2, df %>% filter(k_id == k) %>% pull(y_ANA))
    
    MAE1_bs[k,i] = mean(abs(s1 - (df %>% filter(k_id == k) %>% pull(y_ANA))))
    MAE2_bs[k,i] = mean(abs(s2 - (df %>% filter(k_id == k) %>% pull(y_ANA))))
    
    AIC1_bs[k,i] = AIC(model1)
    AIC2_bs[k,i] = AIC(model2)
    
    BIC1_bs[k,i] = BIC(model1)
    BIC2_bs[k,i] = BIC(model2)
  }
}

AUCs1 %>% 
  .[, 1:3] %>% 
  data.frame() %>% 
  `colnames<-`(c(paste0("d", 1:3))) %>% 
  summarise_all(mean) %>% 
  bind_rows(
    AUCs2 %>% 
      .[, 1:3] %>% 
      data.frame() %>% 
      `colnames<-`(c(paste0("d", 1:3))) %>% 
      summarise_all(mean)
  ) %>% 
  bind_rows(
    AUCs1_bs %>% 
      data.frame() %>% 
      `colnames<-`(c(paste0("d", 1:3))) %>% 
      summarise_all(mean)
  ) %>% 
  bind_rows(
    AUCs2_bs %>% 
      data.frame() %>% 
      `colnames<-`(c(paste0("d", 1:3))) %>% 
      summarise_all(mean)
  ) %>% 
  `rownames<-`(c("Linear", "Logit", "Linear spline", "Logit spline")) %>%
  knitr::kable(format = "pipe", digits = 3) %>%
  writeLines(here("05_Static_Game_Airline_Ch09/output/tab_AUC.txt"))


# 平均絶対誤差（MAE）の集計
MAE1 %>% 
  data.frame() %>% 
  `colnames<-`(c(paste0("d", 1:4))) %>% 
  summarise_all(mean) %>% 
  bind_rows(
    MAE2 %>% 
      data.frame() %>% 
      `colnames<-`(c(paste0("d", 1:4))) %>% 
      summarise_all(mean)
  ) %>% 
  bind_rows(
    MAE1_bs %>% 
      data.frame() %>% 
      `colnames<-`(c(paste0("d", 1:3))) %>% 
      summarise_all(mean)
  ) %>% 
  bind_rows(
    MAE2_bs %>% 
      data.frame() %>% 
      `colnames<-`(c(paste0("d", 1:3))) %>% 
      summarise_all(mean)
  ) %>% 
  `rownames<-`(c("Linear", "Logit", "Linear spline", "Logit spline")) %>%
  knitr::kable(format = "pipe", digits = 3) %>%
  writeLines(here("05_Static_Game_Airline_Ch09/output/tab_MAE.txt"))


# AICの集計
AIC1 %>%
  data.frame() %>%
  `colnames<-`(c(paste0("d", 1:4))) %>%
  summarise_all(mean) %>%
  bind_rows(
    AIC2 %>%
      data.frame() %>%
      `colnames<-`(c(paste0("d", 1:4))) %>%
      summarise_all(mean)
  ) %>%
  bind_rows(
    AIC1_bs %>%
      data.frame() %>%
      `colnames<-`(c(paste0("d", 1:3))) %>%
      summarise_all(mean)
  ) %>%
  bind_rows(
    AIC2_bs %>%
      data.frame() %>%
      `colnames<-`(c(paste0("d", 1:3))) %>%
      summarise_all(mean)
  ) %>%
  `rownames<-`(c("Linear", "Logit", "Linear spline", "Logit spline")) %>%
  knitr::kable(format = "pipe", digits = 3) %>%
  writeLines(here("05_Static_Game_Airline_Ch09/output/tab_AIC.txt"))


# BICの集計
BIC1 %>%
  data.frame() %>%
  `colnames<-`(c(paste0("d", 1:4))) %>%
  summarise_all(mean) %>%
  bind_rows(
    BIC2 %>%
      data.frame() %>%
      `colnames<-`(c(paste0("d", 1:4))) %>%
      summarise_all(mean)
  ) %>%
  bind_rows(
    BIC1_bs %>%
      data.frame() %>%
      `colnames<-`(c(paste0("d", 1:3))) %>%
      summarise_all(mean)
  ) %>%
  bind_rows(
    BIC2_bs %>%
      data.frame() %>%
      `colnames<-`(c(paste0("d", 1:3))) %>%
      summarise_all(mean)
  ) %>%
  `rownames<-`(c("Linear", "Logit", "Linear spline", "Logit spline")) %>%
  knitr::kable(format = "pipe", digits = 3) %>%
  writeLines(here("05_Static_Game_Airline_Ch09/output/tab_BIC.txt"))


# 1st stage specification 1: linear degree 1
m1st_lin_ANA = df %>% 
  lm(y_ANA ~ Train*poly(Distance, Population, Flight_ANA, Flight_JAL, degree=1), data=.)

m1st_lin_JAL = df %>% 
  lm(y_JAL ~ Train*poly(Distance, Population, Flight_ANA, Flight_JAL, degree=1), data=.)

# 1st stage specification 2: logit with b-spline degree 2
m1st_logit_spline_ANA = df %>% 
  glm(generate_bspline_formula(vec, "ANA", degree = 2), family = binomial, data = .)

m1st_logit_spline_JAL = df %>% 
  glm(generate_bspline_formula(vec, "JAL", degree = 2), family = binomial, data = .)

# 1st stage specification 3: logit degree 1
m1st_logit_ANA = df %>% 
  glm(y_ANA ~ Train*poly(Distance, Population, Flight_ANA, Flight_JAL, degree=1),
      family = binomial,
      data = .)

m1st_logit_JAL = df %>% 
  glm(y_JAL ~ Train*poly(Distance, Population, Flight_ANA, Flight_JAL, degree=1),
      family = binomial,
      data = .)

# 各モデル（線形・ロジット・ロジット＋スプライン）で ANA/JALの参入確率を予測
# 線形モデルの予測値が 0〜1 の範囲外に出た場合は切り詰め
df = df %>% 
  mutate(p_lin_ANA = predict(m1st_lin_ANA, newdata = ., type = "response"),
         p_lin_ANA_max = max(ifelse(p_lin_ANA > 1 | p_lin_ANA < 0, NA, p_lin_ANA), na.rm = T),
         p_lin_ANA_min = min(ifelse(p_lin_ANA > 1 | p_lin_ANA < 0, NA, p_lin_ANA), na.rm = T),
         p_lin_ANA = pmax(pmin(p_lin_ANA, p_lin_ANA_max), p_lin_ANA_min),
         p_lin_JAL = predict(m1st_lin_JAL, newdata = ., type = "response"),
         p_lin_JAL_max = max(ifelse(p_lin_JAL > 1 | p_lin_JAL < 0, NA, p_lin_JAL), na.rm = T),
         p_lin_JAL_min = min(ifelse(p_lin_JAL > 1 | p_lin_JAL < 0, NA, p_lin_JAL), na.rm = T),
         p_lin_JAL = pmax(pmin(p_lin_JAL, p_lin_JAL_max), p_lin_JAL_min),
         p_logit_sp_ANA = predict(m1st_logit_spline_ANA, newdata = ., type = "response") ,
         p_logit_sp_JAL = predict(m1st_logit_spline_JAL, newdata = ., type = "response"),
         p_logit_ANA = predict(m1st_logit_ANA, newdata = ., type = "response") ,
         p_logit_JAL = predict(m1st_logit_JAL, newdata = ., type = "response"),
  ) %>% 
  select(-ends_with("max"), -ends_with("min")) 

# long形式に整理
df_long = df %>% 
  mutate(p_lin_ANA_opp = p_lin_ANA,
         p_lin_JAL_opp = p_lin_JAL,
         p_logit_ANA_opp = p_logit_ANA,
         p_logit_JAL_opp = p_logit_JAL,
         p_logit_sp_ANA_opp = p_logit_sp_ANA,
         p_logit_sp_JAL_opp = p_logit_sp_JAL) %>% 
  pivot_longer(cols = c(p_lin_ANA, p_lin_JAL,
                        p_logit_ANA, p_logit_JAL,
                        p_logit_sp_ANA, p_logit_sp_JAL), values_to = "value")%>% 
  mutate(comp = str_sub(name, start = -3),
         spec = str_sub(name, start = 3, end = -5),
  ) %>% 
  select(-name) %>% 
  pivot_wider(values_from = value, names_from = spec) %>% 
  mutate(p_lin_opp = ifelse(comp == "ANA", p_lin_JAL_opp, p_lin_ANA_opp),
         p_logit_opp = ifelse(comp == "ANA", p_logit_JAL_opp, p_logit_ANA_opp),
         p_logit_sp_opp = ifelse(comp == "ANA", p_logit_sp_JAL_opp, p_logit_sp_ANA_opp)) %>% 
  select(-p_lin_JAL_opp, -p_lin_ANA_opp, -p_logit_JAL_opp, -p_logit_ANA_opp,
         -p_logit_sp_JAL_opp, -p_logit_sp_ANA_opp)



## 結果 (linear & logit)
texreg::screenreg(list(m1st_lin_ANA, m1st_logit_ANA, m1st_lin_JAL, m1st_logit_JAL),
                  custom.header = list("y_ANA" = 1:2, "y_JAL" = 3:4),
                  custom.model.names = c("linear", "logit", "linear", "logit"),
                  custom.coef.names = c("Intercept", "Train", "Distance", "Population", 
                                        "Flight_ANA",  "Flight_JAL", "Train:Distance", "Train:Population",
                                        "Train:Flight_ANA", "Train:Flight_JAL"),
                  file = here("05_Static_Game_Airline_Ch09/output/tab_result_linear_logit.txt"))


# ステップ1の可視化
# 便数と確率の関係
flight_probability_plot <- df %>% 
  summarise_all(mean) %>% 
  select(-Flight_ANA) %>%
  expand_grid(Flight_ANA = seq(min(df$Flight_ANA), max(df$Flight_ANA), 0.1)) %>%
  mutate(pred_logit_ANA = predict(m1st_logit_ANA, newdata = ., type = "response"),
         pred_lin_ANA = predict(m1st_lin_ANA, newdata = ., type = "response"),
         pred_lin_ANA = pmin(1, pmax(0, pred_lin_ANA)),
         pred_logit_sp_ANA = predict(m1st_logit_spline_ANA, newdata = ., type = "response"),
         pred_logit_JAL = predict(m1st_logit_JAL, newdata = ., type = "response"),
         pred_lin_JAL = predict(m1st_lin_JAL, newdata = ., type = "response"),
         pred_lin_JAL = pmin(1, pmax(0, pred_lin_JAL)),
         pred_logit_sp_JAL = predict(m1st_logit_spline_JAL, newdata = ., type = "response")
  ) %>%
  pivot_longer(cols = c(pred_logit_ANA, pred_lin_ANA, pred_logit_sp_ANA, 
                        pred_logit_JAL, pred_lin_JAL, pred_logit_sp_JAL), 
               values_to = "value") %>%
  mutate(name_comp = str_sub(name,start = -3) %>% paste0("y_",.),
         flight_comp = "FLight_ANA"
  )　%>% 
  rename(Flight = Flight_ANA) %>%
  bind_rows(
    df %>% 
      summarise_all(mean) %>% 
      select(-Flight_JAL) %>%
      expand_grid(Flight_JAL = seq(min(df$Flight_JAL), max(df$Flight_JAL), 0.1)) %>%
      mutate(pred_logit_ANA = predict(m1st_logit_ANA, newdata = ., type = "response"),
             pred_lin_ANA = predict(m1st_lin_ANA, newdata = ., type = "response"),
             pred_lin_ANA = pmin(1, pmax(0, pred_lin_ANA)),
             pred_logit_sp_ANA = predict(m1st_logit_spline_ANA, newdata = ., type = "response"),
             pred_logit_JAL = predict(m1st_logit_JAL, newdata = ., type = "response"),
             pred_lin_JAL = predict(m1st_lin_JAL, newdata = ., type = "response"),
             pred_lin_JAL = pmin(1, pmax(0, pred_lin_JAL)),
             pred_logit_sp_JAL = predict(m1st_logit_spline_JAL, newdata = ., type = "response")
      ) %>%
      pivot_longer(cols = c(pred_logit_ANA, pred_lin_ANA, pred_logit_sp_ANA, 
                            pred_logit_JAL, pred_lin_JAL, pred_logit_sp_JAL),
                   values_to = "value") %>%
      mutate(name_comp = str_sub(name,start = -3) %>% paste0("y_",.),
             flight_comp = "FLight_JAL"
      ) %>% 
      rename(Flight = Flight_JAL)
  ) %>% 
  mutate(name = str_sub(name, end = -5))

p <- ggplot(flight_probability_plot) +
  geom_line(aes(Flight, value, col = name)) +
  facet_wrap(flight_comp ~ name_comp) +
  labs(y = "y") +
  theme_minimal()

ggsave(filename = here("05_Static_Game_Airline_Ch09/output/fig_flight_probability.pdf"), plot = p)


# 予測値と実測の関係
first_stage_plot <- df %>% 
  pivot_longer(cols = c(p_lin_ANA, p_logit_ANA, p_logit_sp_ANA, y_ANA,
                        p_lin_JAL, p_logit_JAL, p_logit_sp_JAL, y_JAL,
  ),values_to = "value") %>%
  pivot_longer(cols = c(Flight_ANA, Flight_JAL), values_to = "Flight", names_to = "name_flight") %>%
  mutate(val_comp = str_sub(name, start = -3) %>% paste0("y_", .),
         flight_comp = str_sub(name_flight, start = -3) %>% paste0("Flight_", .),
         y_val = ifelse(str_detect(name, "y_"), value, NA),
         name = str_sub(name, end = -5),
  ) 

p <- ggplot(first_stage_plot) +
  geom_point(aes(Flight, value, col = name)) +
  geom_smooth(aes(Flight, y_val), se = F, col = "red") +
  facet_wrap(flight_comp ~ val_comp)+
  labs(y = "y") +
  theme_minimal()

ggsave(filename = here("05_Static_Game_Airline_Ch09/output/fig_first_stage.pdf"), plot = p)


## ステップ2 ----

### 1. 重み付き最小二乗法（Weighted OLS） ----

df_long = df_long %>% 
  mutate(
    Pi_lin = log(lin) - log(1 - lin),
    Pi_logit = log(logit) - log(1 - logit),
    Pi_logit_sp = log(logit_sp) - log(1 - logit_sp)
  )

# 構造パラメータの推定

df_long = df_long %>%
  mutate(Flight = ifelse(comp == "ANA", Flight_ANA, Flight_JAL),
         Flight_opp = ifelse(comp == "ANA", Flight_JAL, Flight_ANA)
  ) 

#linear
est1_1 = df_long %>%
  lm(Pi_lin ~ Train*(Distance + Population + Pop_Square + Train + Flight) + p_lin_opp, data = .) 

#logit
est1_2 = df_long %>%
  lm(Pi_logit ~ Train*(Distance + Population + Pop_Square + Train + Flight) + p_logit_opp, data = .) 

#logit with spline
est1_3 = df_long %>%
  lm(Pi_logit_sp ~ Train*(Distance + Population + Pop_Square + Train + Flight) + p_logit_sp_opp, data = .)

texreg::screenreg(list(est1_1, est1_2, est1_3),
                  custom.model.names = c("linear","logit","logit with spline"),
                  custom.coef.map = list("(Intercept)" = NA,
                                         "Distance" = NA,
                                         "Train" = NA,
                                         "Population" = NA,
                                         "Flight" = NA,
                                         "Pop_Square" = NA,
                                         "Train:Distance" = NA,
                                         "Train:Population" = NA,
                                         "Train:Flight" = NA,
                                         "Train:Pop_Square" = NA,
                                         "p_lin_opp" = "Delta",
                                         "p_logit_opp" = "Delta",
                                         "p_logit_sp_opp" = "Delta"),
                  file = here("05_Static_Game_Airline_Ch09/output/tab_result_linear_logit_spline.txt"))


### 2. GMM  ----

df_long = df_long %>% 
  mutate(y = ifelse(comp == "ANA", y_ANA, y_JAL),
         Distance_Train = Distance*Train,
         Population_Train = Population*Train,
         Pop_Square_Train = Pop_Square*Train,
         Flight_Train = Flight*Train
  ) %>%
  select(-c(y_ANA, y_JAL))

Z = df_long %>% 
  select(Constant, Distance , Population , Pop_Square , Train , Flight , 
         Distance_Train, Population_Train, Pop_Square_Train, Flight_Train, Flight_opp) %>%
  as.matrix()

X_lin = df_long %>% 
  select(Constant, Distance , Population , Pop_Square , Train , Flight , 
         Distance_Train, Population_Train, Pop_Square_Train, Flight_Train, p_lin_opp) %>%
  as.matrix()

X_logit = df_long %>% 
  select(Constant, Distance , Population , Pop_Square , Train , Flight , 
         Distance_Train, Population_Train, Pop_Square_Train, Flight_Train, p_logit_opp) %>% 
  as.matrix()

X_logit_sp = df_long %>% 
  select(Constant, Distance , Population , Pop_Square , Train , Flight , 
         Distance_Train, Population_Train, Pop_Square_Train, Flight_Train, p_logit_sp_opp) %>% 
  as.matrix()

y_ = df_long %>% pull(y)


param0 = rep(0,ncol(X_lin))
names(param0) = c("(Intercept)", "Distance", "Population", "Pop_Square", "Train", "Flight",
                  "Train:Distance", "Train:Population", "Train:Pop_Square", "Train:Flight", "p_opp")

#not IV
est2_1_gmm = nleqslv(
  fn = function(x) calculate_moment(par = x, X = X_lin, Z = X_lin, y = y_),
  x = param0,
)

est2_2_gmm = nleqslv(
  fn = function(x) calculate_moment(par = x, X = X_logit, y = y_, Z = X_logit),
  x = param0,
)

est2_3_gmm = nleqslv(
  fn = function(x) calculate_moment(par = x, X = X_logit_sp, y = y_, Z = X_logit_sp),
  x = param0,
)


#IV estimator
est2_1_gmm2 = nleqslv(
  fn = function(x) calculate_moment(par = x, X = X_lin, y = y_, Z = Z),
  x = param0,
)

est2_2_gmm2 = nleqslv(
  fn = function(x) calculate_moment(par = x, X = X_logit, y = y_, Z = Z),
  x = param0,
)

est2_3_gmm2 = nleqslv(
  fn = function(x) calculate_moment(par = x, X = X_logit_sp, y = y_, Z = Z),
  x = param0,
)


###  3. 最尤推定 ----

#linear
est3_1 = df_long %>%
  glm(y ~ Train*(Distance + Population + Pop_Square + Flight) + p_lin_opp,
      family = binomial(link = "logit"), data = .)


#logit
est3_2 = df_long %>%
  glm(y ~  Train*(Distance + Population + Pop_Square + Flight) + p_logit_opp,
      family = binomial(link = "logit"), data = .)

#logit with spline
est3_3 = df_long %>%
  glm(y ~ Train*(Distance + Population + Pop_Square + Flight) + p_logit_sp_opp,
      family = binomial(link = "logit"), data = .)

# 最尤推定の結果の比較
texreg::screenreg(list(est3_1, est3_2, est3_3),
                  custom.model.names = c("linear","logit","logit with spline"),
                  custom.coef.map = list("(Intercept)" = NA,
                                         "Distance" = NA,
                                         "Train" = NA,
                                         "Population" = NA,
                                         "Flight" = NA,
                                         "Pop_Square" = NA,
                                         "Train:Distance" = NA,
                                         "Train:Population" = NA,
                                         "Train:Flight" = NA,
                                         "Train:Pop_Square" = NA,
                                         "p_lin_opp" = "Delta",
                                         "p_logit_opp" = "Delta",
                                         "p_logit_sp_opp" = "Delta"),
                  file = here("05_Static_Game_Airline_Ch09/output/tab_result_MLE.txt"))

# 重み付き最小二乗法と最尤推定の結果の比較
# 英語
texreg::screenreg(list(est1_1, est1_2, est3_1, est3_2),
               custom.header = list("OLS (y = Pi)" = 1:2, "MLE (y = action)" = 3:4),
               custom.model.names = c("linear", "logit", "linear", "logit"),
               custom.coef.map = list("(Intercept)" = NA,
                                      "Distance" = NA,
                                      "Train" = NA,
                                      "Population" = NA,
                                      "Flight" = NA,
                                      "Pop_Square" = NA,
                                      "Train:Distance" = NA,
                                      "Train:Population" = NA,
                                      "Train:Flight" = NA,
                                      "Train:Pop_Square" = NA,
                                      "p_lin_opp" = "Delta",
                                      "p_logit_opp" = "Delta",
                                      "p_logit_sp_opp" = "Delta"),
               file = here("05_Static_Game_Airline_Ch09/output/tab_result_OLS_MLE_Eng.txt"))



# 日本語
texreg::screenreg(list(est1_1, est1_2, est3_1, est3_2),
               custom.header = list("OLS (y = Pi)" = 1:2, "MLE (y = action)" = 3:4),
               custom.model.names = c("linear","logit", "linear", "logit"),
               custom.coef.map = list("(Intercept)" = "定数項" ,
                                      "Distance" = "距離",
                                      "Train" = "新幹線",
                                      "Population" = "人口",
                                      "Pop_Square" = "人口2乗",
                                      "Flight" = "フライト",
                                      "Train:Distance" = "新幹線:距離",
                                      "Train:Population" = "新幹線:人口",
                                      "Train:Pop_Square" = "新幹線:人口2乗",
                                      "Train:Flight" = "新幹線:フライト",
                                      "p_lin_opp" = "競争",
                                      "p_logit_opp" = "競争",
                                      "p_logit_sp_opp" = "競争"),
               file = here("05_Static_Game_Airline_Ch09/output/tab_result_OLS_MLE_Jpn.txt"))



# moment
tibble(
  var = c(colnames(X_lin), "value"),
  est = c(est2_1_gmm$x, (t(est2_1_gmm$fvec) %*% est2_1_gmm$fvec)[1]),
  se = c(calculate_standard_error(est2_1_gmm$x, X_lin, X_lin, y_),NA),
  f_model = "linear",
  s_model = "Moment"
) %>% 
  bind_rows(
    tibble(
      var = c(colnames(X_lin), "value"),
      est = c(est2_2_gmm$x,(t(est2_2_gmm$fvec) %*% est2_2_gmm$fvec)[1]),
      se = c(calculate_standard_error(est2_2_gmm$x,X_logit,X_logit, y_),NA),
      f_model = "logit",
      s_model = "Moment"
    )
  ) %>%
  bind_rows(
    tibble(
      var = c(colnames(X_lin), "value"),
      est = c(est2_3_gmm$x, (t(est2_3_gmm$fvec) %*% est2_3_gmm$fvec)[1]),
      se = c(calculate_standard_error(est2_3_gmm$x, X_logit_sp, X_logit_sp, y_), NA),
      f_model = "logit_spline",
      s_model = "Moment"
    )
  )  %>% 
  #IV moment
  bind_rows(
    tibble(
      var = c(colnames(X_lin),"value"),
      est = c(est2_1_gmm2$x, (t(est2_1_gmm2$fvec) %*% est2_1_gmm2$fvec)[1]),
      se = c(calculate_standard_error(est2_1_gmm2$x, X_logit, Z, y_), NA),
      f_model = "linear",
      s_model = "IV_Moment"
    )
  ) %>%  
  bind_rows(
    tibble(
      var = c(colnames(X_lin), "value"),
      est = c(est2_2_gmm2$x, (t(est2_2_gmm2$fvec) %*% est2_2_gmm2$fvec)[1]),
      se = c(calculate_standard_error(est2_2_gmm2$x, X_logit,Z, y_), NA),
      f_model = "logit",
      s_model = "IV_Moment"
    )
  ) %>% 
  bind_rows(
    tibble(
      var = c(colnames(X_lin), "value"),
      est = c(est2_3_gmm2$x, (t(est2_3_gmm2$fvec) %*% est2_3_gmm2$fvec)[1]),
      se = c(calculate_standard_error(est2_3_gmm2$x, X_logit_sp, Z, y_), NA),
      f_model = "logit_spline",
      s_model = "IV_Moment"
    )
  ) %>% 
  #OLS
  bind_rows(
    tibble(
      var = names(summary(est1_1)$coefficients[,1]),
      est = summary(est1_1)$coefficients[,1],
      se = summary(est1_1)$coefficients[,2],
      f_model = "linear",
      s_model = "OLS"
    )
  ) %>%
  bind_rows(
    tibble(
      var = names(summary(est1_2)$coefficients[,1]),
      est = summary(est1_2)$coefficients[,1],
      se = summary(est1_2)$coefficients[,2],
      f_model = "logit",
      s_model = "OLS"
    )
  ) %>%
  bind_rows(
    tibble(
      var = names(summary(est1_3)$coefficients[,1]),
      est = summary(est1_3)$coefficients[,1],
      se = summary(est1_3)$coefficients[,2],
      f_model = "logit_spline",
      s_model = "OLS"
    )
  ) %>%
  #MLE
  bind_rows(
    tibble(
      var = names(summary(est3_1)$coefficients[,1]),
      est = summary(est3_1)$coefficients[,1],
      se = summary(est3_1)$coefficients[,2],
      f_model = "linear",
      s_model = "MLE"
    )
  ) %>%
  bind_rows(
    tibble(
      var = names(summary(est3_2)$coefficients[,1]),
      est = summary(est3_2)$coefficients[,1],
      se = summary(est3_2)$coefficients[,2],
      f_model = "logit",
      s_model = "MLE"
    )
  ) %>%
  bind_rows(
    tibble(
      var = names(summary(est3_3)$coefficients[,1]),
      est = summary(est3_3)$coefficients[,1],
      se = summary(est3_3)$coefficients[,2],
      f_model = "logit_spline",
      s_model = "MLE"
    )
  ) %>%
  filter(var != "value") %>% 
  mutate(var = 
           # 英語
           case_when(
             var == "Constant" ~ "(Intercept)",
             var == "Distance_Train" ~ "Train:Distance",
             var == "Population_Train" ~ "Train:Population",
             var == "Pop_Square_Train" ~ "Train:Pop_Square",
             var == "Flight_Train" ~ "Train:Flight",
             var %in% c("p_lin_opp","p_logit_opp","p_logit_sp_opp") ~ "Delta",
             TRUE ~ var
           ),
         # 日本語
         var = case_when(
           var == "(Intercept)" ~ "定数",
           var == "Distance" ~ "距離",
           var == "Train" ~ "新幹線",
           var == "Population" ~ "人口",
           var == "Pop_Square" ~ "人口2乗",
           var == "Flight" ~ "フライト",
           var == "Delta" ~ "競争",
           var == "Train:Distance" ~ "新幹線:距離",
           var == "Train:Population" ~ "新幹線:人口",
           var == "Train:Pop_Square" ~ "新幹線:人口2乗",
           var == "Train:Flight" ~ "新幹線:フライト",
           TRUE ~ var
         )
  ) %>% 
  pivot_longer(cols = c(est,se),names_to = "type",values_to = "value") %>% 
  mutate(
    # 小数点以下2桁に丸める
    # 末尾の数字が0の場合も表示する
    # type が se のときは () を付ける
    value = ifelse(type == "est", sprintf("%.2f",value), paste0(" (", sprintf("%.2f", value), ")")),
  ) %>% 
  pivot_wider(names_from = c(s_model, f_model), values_from = c(value)) %>% 
  mutate(var = ifelse(type == "se", "", var)) %>% 
  select(var,
         starts_with("OLS"),
         starts_with("MLE"),
         starts_with("Moment"),
         starts_with("IV")
  ) %>% 
  bind_rows(
    tibble(
      var = c("Num. obs.","R^2","Log Likelihood"),
      OLS_linear = c(nobs(est1_1), 
                     round(summary(est1_1)$r.squared,2), 
                     caluculate_log_likelihood(est1_1$coefficients)),
      OLS_logit = c(nobs(est1_2), 
                    round(summary(est1_2)$r.squared,2), 
                    caluculate_log_likelihood(est1_2$coefficients)),
      OLS_logit_spline = c(nobs(est1_3),
                           round(summary(est1_3)$r.squared,2),
                           caluculate_log_likelihood(est1_3$coefficients)),
      MLE_linear = c(nobs(est3_1),
                     NA,
                     caluculate_log_likelihood(est3_1$coefficients)),
      MLE_logit = c(nobs(est3_2),
                    NA,
                    caluculate_log_likelihood(est3_2$coefficients)),
      MLE_logit_spline = c(nobs(est3_3),
                           NA,
                           caluculate_log_likelihood(est3_3$coefficients)),
      Moment_linear = c(nrow(df_long),
                        NA,
                        caluculate_log_likelihood(par = est2_1_gmm$x, f_spec = "p_lin_opp")),
      Moment_logit = c(nrow(df_long),
                       NA,
                       caluculate_log_likelihood(est2_2_gmm$x, f_spec = "p_logit_opp")),
      Moment_logit_spline = c(nrow(df_long),
                              NA,
                              caluculate_log_likelihood(est2_3_gmm$x, f_spec = "p_logit_sp_opp")),
      IV_Moment_linear = c(nrow(df_long),
                           NA,
                           caluculate_log_likelihood(est2_1_gmm2$x, f_spec = "p_lin_opp")),
      IV_Moment_logit = c(nrow(df_long),
                          NA,
                          caluculate_log_likelihood(est2_2_gmm2$x, f_spec = "p_logit_opp")),
      IV_Moment_logit_spline = c(nrow(df_long),
                                 NA,
                                 caluculate_log_likelihood(est2_3_gmm2$x, f_spec = "p_logit_sp_opp"))
    ) %>% 
      mutate(
        across(c(-var), ~  ifelse(var == "Num. obs.", as.character(.x), sprintf("%.2f", .x)))
      )
  ) %>% 
  knitr::kable(format = "pipe", digits = 3) %>%
  writeLines(here("05_Static_Game_Airline_Ch09/output/tab9_2_estimation.txt"))


# 均衡の計算 ----

df_long = df_long %>% 
  arrange(comp)


df_long = df_long %>% 
  mutate(
    new_eq_ana = find_equilibrium_ANA_JAL(unique(df$id), "ANA", est2_2_gmm2$x, df_long = df_long),
    new_eq_jal = find_equilibrium_ANA_JAL(unique(df$id), "JAL", est2_2_gmm2$x, df_long = df_long),
    dif = abs(new_eq_ana - new_eq_jal)
  )

# 複数均衡の確認
p <- df_long %>% 
  ggplot() +
  geom_point(aes(x = new_eq_ana, y = new_eq_jal,color = dif))+
  geom_abline(intercept = 0, slope = 1)+
  scale_color_gradient2(low = "blue", high = "red", midpoint = .01)+
  facet_wrap(~comp) +
  theme_minimal()

ggsave(filename = here("05_Static_Game_Airline_Ch09/output/fig_equilibrium_multiplicity.pdf"), plot=p)

# ステップ1のモデルフィットの確認

df_long_fit_plot <- df_long %>% 
  select(id,comp,new_eq_jal, p_lin_opp) %>% 
  pivot_wider(names_from = comp, values_from = c(new_eq_jal, p_lin_opp)) %>%
  rename(ANA_eq = new_eq_jal_ANA, JAL_eq = new_eq_jal_JAL, 
         JAL_p_1st = p_lin_opp_ANA, ANA_p_1st = p_lin_opp_JAL) %>% 
  pivot_longer(cols = c(ANA_eq, JAL_eq), names_to = "comp", values_to = "eq") %>% 
  pivot_longer(cols = c(JAL_p_1st, ANA_p_1st), names_to = "comp_1st", values_to = "p_1st") %>% 
  mutate(comp = str_replace(comp, "_eq", ""),
         comp_1st = str_replace(comp_1st, "_p_1st", "")) %>%
  filter(comp == comp_1st) %>% 
  left_join(
    df_long %>% 
      mutate(p_2nd = calculate_predict_probability(coef(est1_2))) %>% 
      select(id,comp, p_2nd), by = c("id","comp")
  ) %>%
  pivot_longer(cols = c(p_2nd,p_1st), names_to = "p_name", values_to = "p_value")

p <- df_long_fit_plot %>% 
  ggplot() +
  geom_point(aes(x = eq, y = p_value, color = abs(p_value - eq)))+
  geom_abline(intercept = 0, slope = 1)+
  scale_color_gradient(low = "blue", high = "red")+
  facet_wrap(p_name~comp) +
  theme_minimal() 

ggsave(filename = here("05_Static_Game_Airline_Ch09/output/fig_model_fit.pdf"), plot=p)


# 反実仮想: 北陸新幹線の通る空港 ----

df_long = df_long %>% 
  left_join(id_port_dyad_long, by = "id")


df_long_hokuriku = df_long %>% 
  mutate(
    hokuriku = case_when(
      port1 %in% c("小松","富山","松本","羽田","伊丹") & port2 %in% c("小松","富山","松本","羽田","伊丹") ~ 1,
      TRUE ~ 0
    )
  ) %>% 
  mutate(Train_orig = Train,
         Train = pmin(1, Train + hokuriku),
         Distance_Train = Distance*Train,
         Population_Train = Population*Train,
         Pop_Square_Train = Pop_Square*Train,
         Flight_Train = Flight*Train)


# シミュレーション
# パラメータはOLS-ロジットの推定結果（est1_2）を利用
df_long_hokuriku = df_long_hokuriku %>% 
  mutate(
    orig_eq_ana = find_equilibrium_ANA_JAL(
      unique(df$id),
      "ANA",
      est1_2$coefficients,
      df_long = df_long),
    orig_eq_jal = find_equilibrium_ANA_JAL(
      unique(df$id),
      "JAL",
      est1_2$coefficients,
      df_long = df_long),
    hoku_eq_ana = find_equilibrium_ANA_JAL(
      unique(df$id),
      "ANA",
      est1_2$coefficients,
      df_long = df_long_hokuriku),
    hoku_eq_jal = find_equilibrium_ANA_JAL(
      unique(df$id),
      "JAL",
      est1_2$coefficients,
      df_long = df_long_hokuriku)
  )


# port_orderを作成
p2_order <- id_port_dyad_long %>%
  filter(port2 != "静岡") %>%
  distinct(port2) %>%
  pull(port2)

south_tail <- id_port_dyad_long %>%
  filter(port1 != "静岡") %>%
  anti_join(tibble(port2 = p2_order), by = c("port1" = "port2")) %>%
  distinct(port1) %>%
  pull(port1)

port_order <- c(p2_order, south_tail) %>% unique()

# 均衡での参入確率を可視化
df_hokuriku_plot <- df_long_hokuriku %>% 
  filter(hokuriku == 1) %>% 
  mutate(port1_b = str_sub(dyad, regexpr("_", dyad) + 1, nchar(dyad)),
         port2_b = str_sub(dyad, 1, regexpr("_", dyad) - 1)) %>% 
  pivot_longer(c(port1, port1_b), values_to = "port1") %>% 
  mutate(port2 = ifelse(name == "port1", port2, port2_b)) %>% 
  pivot_longer(c(orig_eq_ana, hoku_eq_ana), values_to = "eq", names_to = "counterfactual") %>%
  bind_rows(
    expand_grid(
      tibble(port1 = port_order,
             port2 = port_order),
      eq = NA,
      counterfactual = c("orig_eq_ana","hoku_eq_ana"),
      comp = c("ANA","JAL")
    )%>% 
      mutate(hokuriku = case_when(
        port1 %in% c("小松","富山","松本","羽田","伊丹") & port2 %in% c("小松","富山","松本","羽田","伊丹") ~ 1,
        TRUE ~ 0
      ))
  ) %>%
  filter(hokuriku == 1) %>% 
  mutate(
    port1 = factor(port1, levels = port_order),
    port1_id = as.numeric(factor(port1, levels = port_order)),
    port2 = factor(port2, levels = port_order),
    port2_id = as.numeric(factor(port2, levels = port_order)),
    eq = ifelse(port1_id >= port2_id,eq,NA),
    counterfactual = ifelse(counterfactual == "orig_eq_ana", "original", "hokuriku_shinkansen"),
    counterfactual = factor(counterfactual, levels = c("original", "hokuriku_shinkansen")),
  ) 

p <- df_hokuriku_plot %>%
  ggplot(aes(port1, port2, fill = eq)) + 
  geom_tile(color = "white") +
  facet_wrap(counterfactual~comp,labeller = labeller(counterfactual = 
                                                       c("original" = "北陸新幹線なし",
                                                         "hokuriku_shinkansen" = "北陸新幹線あり")
  )) +
  labs(x = "", y = "") +
  scale_fill_gradient2(name = "均衡参入確率") + 
  theme_gray(base_family = "HiraKakuPro-W3")

ggsave(here("05_Static_Game_Airline_Ch09/output/fig_hokuriku_heatmap.pdf"), plot = p, device = cairo_pdf)


df_long_hokuriku %>% 
  filter(hokuriku == 1) %>% 
  select(dyad,comp, orig_eq_ana, hoku_eq_ana) %>% 
  rename("ルート" = dyad,
         "企業" = comp,
         "北陸新幹線あり" = hoku_eq_ana,
         "北陸新幹線なし" = orig_eq_ana) %>% 
  pivot_wider(names_from = "企業", values_from = c("北陸新幹線あり", "北陸新幹線なし")) %>% 
  relocate(ルート, 北陸新幹線なし_ANA, 北陸新幹線あり_ANA, 北陸新幹線なし_JAL, 北陸新幹線あり_JAL) %>%
  knitr::kable(format = "pipe", digits = 3) %>%
  writeLines(here("05_Static_Game_Airline_Ch09/output/tab9_3_counter_factual.txt"))
  


