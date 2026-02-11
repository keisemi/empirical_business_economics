"""Generate main_ch09.ipynb for Chapter 9: Static Game - Airline Entry."""
import nbformat as nbf

nb = nbf.v4.new_notebook()
nb.metadata.update({
    "kernelspec": {
        "display_name": "Python 3",
        "language": "python",
        "name": "python3"
    },
    "language_info": {
        "name": "python",
        "version": "3.9.0"
    }
})

cells = []

# ============================================================
# Cell: Title
# ============================================================
cells.append(nbf.v4.new_markdown_cell(
    "# 第9章 静学ゲーム：航空会社の参入ゲーム\n"
    "\n"
    "2社（ANA, JAL）の航空路線への参入ゲームを分析する。\n"
    "2段階推定法により構造パラメータを推定し、北陸新幹線の反実仮想分析を行う。"
))

# ============================================================
# Cell: Setup
# ============================================================
cells.append(nbf.v4.new_markdown_cell("## 1. Pythonに関する下準備"))

cells.append(nbf.v4.new_code_cell("""\
import numpy as np
import pandas as pd
from scipy import optimize
from scipy.optimize import fsolve
from sklearn.preprocessing import PolynomialFeatures
from sklearn.metrics import roc_auc_score, mean_absolute_error
import statsmodels.api as sm
from patsy import dmatrix
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib
import warnings
import time

warnings.filterwarnings('ignore')

# 日本語フォント設定
for font_name in ['Hiragino Maru Gothic Pro', 'Hiragino Sans',
                   'IPAexGothic', 'Noto Sans CJK JP', 'Yu Gothic']:
    try:
        matplotlib.font_manager.findfont(font_name, fallback_to_default=False)
        plt.rcParams['font.family'] = font_name
        break
    except ValueError:
        continue

# パス設定
base_dir = Path('..')
data_dir = base_dir / 'data'
output_dir = base_dir / 'output'
output_dir.mkdir(exist_ok=True)

print("Setup complete.")
"""))

# ============================================================
# Cell: Data Loading
# ============================================================
cells.append(nbf.v4.new_markdown_cell("## 2. データの読み込みと準備"))

cells.append(nbf.v4.new_code_cell("""\
# データの読み込み
y = pd.read_csv(data_dir / 'y.csv')
x_ANA = pd.read_csv(data_dir / 'x_ANA.csv')
x_JAL = pd.read_csv(data_dir / 'x_JAL.csv')
id_port_dyad_long = pd.read_csv(data_dir / 'id_port_dyad_long.csv')

# 二つのデータが同じであることを確認（Flight列以外）
cols_check = [c for c in x_ANA.columns if c != 'Flight']
assert x_ANA[cols_check].equals(x_JAL[cols_check]), "x_ANA and x_JAL differ!"
print("x_ANA and x_JAL (Flight以外) は同一: OK")

# マージ
df = y.copy()
df['id'] = range(1, len(df) + 1)
df['Constant'] = x_ANA['Constant'].values
df['Distance'] = x_ANA['Distance'].values
df['Population'] = x_ANA['Population'].values
df['Pop_Square'] = x_ANA['Pop_Square'].values
df['Train'] = x_ANA['Train'].values
df['Flight_ANA'] = x_ANA['Flight'].values
df['Flight_JAL'] = x_JAL['Flight'].values

# id を先頭に
df = df[['id'] + [c for c in df.columns if c != 'id']]

print(f"データサイズ: {df.shape}")
df.head()
"""))

# ============================================================
# Cell: Descriptive Statistics (Tab 9.1)
# ============================================================
cells.append(nbf.v4.new_markdown_cell("## 3. 記述統計 (表9.1)"))

cells.append(nbf.v4.new_code_cell("""\
# 記述統計
desc_vars = ['y_ANA', 'y_JAL', 'Distance', 'Population', 'Pop_Square',
             'Train', 'Flight_ANA', 'Flight_JAL']
desc_stats = df[desc_vars].agg(['mean', 'std']).T
desc_stats.columns = ['Mean', 'SD']
desc_stats.index = ['y_ANA', 'y_JAL', '距離', '人口', '人口2乗',
                     '新幹線', 'フライトANA', 'フライトJAL']

print("表9.1: 記述統計")
print(desc_stats.round(3).to_string())

# 保存
desc_stats.round(3).to_csv(output_dir / 'tab9_1_descriptive_statistics.txt', sep='|')
print("\\n保存先: output/tab9_1_descriptive_statistics.txt")
"""))

# ============================================================
# Cell: Helper Functions for Step 1
# ============================================================
cells.append(nbf.v4.new_markdown_cell(
    "## 4. 推定\n\n"
    "### ステップ1: ノンパラメトリックな確率推定\n\n"
    "#### 補助関数の定義"
))

cells.append(nbf.v4.new_code_cell("""\
def create_poly_features(df_sub, degree):
    \"\"\"多項式特徴量をTrain変数との交差項込みで作成する。
    R の y ~ Train*poly(Distance, Population, Flight_ANA, Flight_JAL, degree=i) に対応。
    \"\"\"
    base_vars = ['Distance', 'Population', 'Flight_ANA', 'Flight_JAL']
    X_base = df_sub[base_vars].values
    train = df_sub['Train'].values.reshape(-1, 1)

    poly = PolynomialFeatures(degree=degree, include_bias=False)
    X_poly = poly.fit_transform(X_base)

    # Train と Train*poly の交差項を追加
    X_train_inter = np.hstack([train, train * X_poly])
    X_all = np.hstack([X_poly, X_train_inter])

    return X_all


def create_bspline_features(df_sub, degree, df_full=None):
    \"\"\"B-spline特徴量をTrain変数との交差項込みで作成する。
    R の y ~ Train*(bs(Distance, ...) + bs(Population, ...) + ...) に対応。
    df_full: knot計算に使うデータ（テスト時はトレーニングデータを渡す）
    \"\"\"
    if df_full is None:
        df_full = df_sub
    base_vars = ['Distance', 'Population', 'Flight_ANA', 'Flight_JAL']
    train = df_sub['Train'].values.reshape(-1, 1)

    # 各変数にB-splineを適用（中央値にknotを1つ置く）
    bs_list = []
    for var in base_vars:
        knot = np.median(df_full[var].values)
        formula = f"bs(x, knots=[{knot}], degree={degree}, include_intercept=False) - 1"
        bs_mat = np.asarray(dmatrix(formula, {"x": df_sub[var].values}))
        bs_list.append(bs_mat)

    X_bs = np.hstack(bs_list)
    # Train との交差項を追加
    X_all = np.hstack([X_bs, train, train * X_bs])

    return X_all


print("補助関数の定義完了")
"""))

# ============================================================
# Cell: Step 1 - Polynomial R-squared
# ============================================================
cells.append(nbf.v4.new_markdown_cell("#### 1. 多項式モデルの決定係数"))

cells.append(nbf.v4.new_code_cell("""\
# ANA: 多項式次数ごとの決定係数
r_sq_ANA = []
for deg in range(1, 6):
    X = sm.add_constant(create_poly_features(df, deg))
    model = sm.OLS(df['y_ANA'].values, X).fit()
    r_sq_ANA.append({'degree': deg, 'Rsq': model.rsquared_adj})

# JAL: 多項式次数ごとの決定係数
r_sq_JAL = []
for deg in range(1, 6):
    X = sm.add_constant(create_poly_features(df, deg))
    model = sm.OLS(df['y_JAL'].values, X).fit()
    r_sq_JAL.append({'degree': deg, 'Rsq': model.rsquared_adj})

r_sq_df = pd.DataFrame(r_sq_ANA).assign(comp='ANA')
r_sq_df = pd.concat([r_sq_df, pd.DataFrame(r_sq_JAL).assign(comp='JAL')])

# 可視化
fig, ax = plt.subplots(figsize=(8, 5))
for comp, grp in r_sq_df.groupby('comp'):
    ax.plot(grp['degree'], grp['Rsq'], 'o-', label=comp)
ax.set_xlabel('degree')
ax.set_ylabel('Adjusted R-squared')
ax.legend()
plt.tight_layout()
plt.savefig(output_dir / 'fig_Rsquared.png', dpi=150)
plt.show()
"""))

# ============================================================
# Cell: Cross-Validation
# ============================================================
cells.append(nbf.v4.new_markdown_cell("#### 交差検証 (10-fold CV)"))

cells.append(nbf.v4.new_code_cell("""\
# 交差検証
np.random.seed(825)
K = 10
df['k_id'] = np.ceil(np.random.uniform(size=len(df)) * K).astype(int)

print("各フォールドのサイズ:")
print(df['k_id'].value_counts().sort_index())

degrees = [1, 2, 3, 4]

# 格納用: 多項式
AUCs1 = np.zeros((K, len(degrees)))   # linear polynomial
AUCs2 = np.zeros((K, len(degrees)))   # logit polynomial
MAE1 = np.zeros((K, len(degrees)))
MAE2 = np.zeros((K, len(degrees)))
AIC1 = np.zeros((K, len(degrees)))
AIC2 = np.zeros((K, len(degrees)))
BIC1 = np.zeros((K, len(degrees)))
BIC2 = np.zeros((K, len(degrees)))

for i, deg in enumerate(degrees):
    for k in range(1, K + 1):
        train_mask = df['k_id'] != k
        test_mask = df['k_id'] == k

        X_train = sm.add_constant(create_poly_features(df[train_mask], deg))
        X_test = sm.add_constant(create_poly_features(df[test_mask], deg))
        y_train = df.loc[train_mask, 'y_ANA'].values
        y_test = df.loc[test_mask, 'y_ANA'].values

        # Linear model
        model1 = sm.OLS(y_train, X_train).fit()
        s1 = np.clip(model1.predict(X_test), 0, 1)

        # Logit model
        try:
            model2 = sm.Logit(y_train, X_train).fit(disp=0, method='bfgs', maxiter=200)
            s2 = model2.predict(X_test)
        except Exception:
            s2 = np.full(len(y_test), 0.5)

        # AUC (degree <= 3 only, matching R: i < 4)
        if deg < 4:
            try:
                AUCs1[k-1, i] = roc_auc_score(y_test, s1)
                AUCs2[k-1, i] = roc_auc_score(y_test, s2)
            except ValueError:
                pass

        MAE1[k-1, i] = mean_absolute_error(y_test, s1)
        MAE2[k-1, i] = mean_absolute_error(y_test, s2)

        # AIC, BIC
        AIC1[k-1, i] = model1.aic
        BIC1[k-1, i] = model1.bic
        try:
            AIC2[k-1, i] = model2.aic
            BIC2[k-1, i] = model2.bic
        except Exception:
            pass

print("多項式CV完了")
"""))

# ============================================================
# Cell: B-spline CV
# ============================================================
cells.append(nbf.v4.new_markdown_cell("#### 2. B-splineモデル"))

cells.append(nbf.v4.new_code_cell("""\
# B-spline CV
degrees_bs = [1, 2, 3]
AUCs1_bs = np.zeros((K, len(degrees_bs)))
AUCs2_bs = np.zeros((K, len(degrees_bs)))
MAE1_bs = np.zeros((K, len(degrees_bs)))
MAE2_bs = np.zeros((K, len(degrees_bs)))
AIC1_bs = np.zeros((K, len(degrees_bs)))
AIC2_bs = np.zeros((K, len(degrees_bs)))
BIC1_bs = np.zeros((K, len(degrees_bs)))
BIC2_bs = np.zeros((K, len(degrees_bs)))

for i, deg in enumerate(degrees_bs):
    for k in range(1, K + 1):
        train_mask = df['k_id'] != k
        test_mask = df['k_id'] == k

        X_train = sm.add_constant(
            create_bspline_features(df[train_mask], deg, df_full=df[train_mask]))
        X_test = sm.add_constant(
            create_bspline_features(df[test_mask], deg, df_full=df[train_mask]))
        y_train = df.loc[train_mask, 'y_ANA'].values
        y_test = df.loc[test_mask, 'y_ANA'].values

        # Linear
        model1 = sm.OLS(y_train, X_train).fit()
        s1 = np.clip(model1.predict(X_test), 0, 1)

        # Logit
        try:
            model2 = sm.Logit(y_train, X_train).fit(disp=0, method='bfgs', maxiter=200)
            s2 = model2.predict(X_test)
        except Exception:
            s2 = np.full(len(y_test), 0.5)

        try:
            AUCs1_bs[k-1, i] = roc_auc_score(y_test, s1)
            AUCs2_bs[k-1, i] = roc_auc_score(y_test, s2)
        except ValueError:
            pass

        MAE1_bs[k-1, i] = mean_absolute_error(y_test, s1)
        MAE2_bs[k-1, i] = mean_absolute_error(y_test, s2)

        AIC1_bs[k-1, i] = model1.aic
        BIC1_bs[k-1, i] = model1.bic
        try:
            AIC2_bs[k-1, i] = model2.aic
            BIC2_bs[k-1, i] = model2.bic
        except Exception:
            pass

print("B-spline CV完了")
"""))

# ============================================================
# Cell: CV Summary Tables (AUC, MAE, AIC, BIC)
# ============================================================
cells.append(nbf.v4.new_markdown_cell("#### CVの結果まとめ"))

cells.append(nbf.v4.new_code_cell("""\
# AUC (degree 1-3 のみ)
auc_table = pd.DataFrame({
    'd1': [AUCs1[:, 0].mean(), AUCs2[:, 0].mean(),
           AUCs1_bs[:, 0].mean(), AUCs2_bs[:, 0].mean()],
    'd2': [AUCs1[:, 1].mean(), AUCs2[:, 1].mean(),
           AUCs1_bs[:, 1].mean(), AUCs2_bs[:, 1].mean()],
    'd3': [AUCs1[:, 2].mean(), AUCs2[:, 2].mean(),
           AUCs1_bs[:, 2].mean(), AUCs2_bs[:, 2].mean()],
}, index=['Linear', 'Logit', 'Linear spline', 'Logit spline'])

print("AUC:")
print(auc_table.round(3).to_string())
auc_table.round(3).to_csv(output_dir / 'tab_AUC.txt', sep='|')

# MAE (degree 1-4 for poly, 1-3 for spline)
mae_table = pd.DataFrame({
    'd1': [MAE1[:, 0].mean(), MAE2[:, 0].mean(),
           MAE1_bs[:, 0].mean(), MAE2_bs[:, 0].mean()],
    'd2': [MAE1[:, 1].mean(), MAE2[:, 1].mean(),
           MAE1_bs[:, 1].mean(), MAE2_bs[:, 1].mean()],
    'd3': [MAE1[:, 2].mean(), MAE2[:, 2].mean(),
           MAE1_bs[:, 2].mean(), MAE2_bs[:, 2].mean()],
    'd4': [MAE1[:, 3].mean(), MAE2[:, 3].mean(), np.nan, np.nan],
}, index=['Linear', 'Logit', 'Linear spline', 'Logit spline'])

print("\\nMAE:")
print(mae_table.round(3).to_string())
mae_table.round(3).to_csv(output_dir / 'tab_MAE.txt', sep='|')

# AIC
aic_table = pd.DataFrame({
    'd1': [AIC1[:, 0].mean(), AIC2[:, 0].mean(),
           AIC1_bs[:, 0].mean(), AIC2_bs[:, 0].mean()],
    'd2': [AIC1[:, 1].mean(), AIC2[:, 1].mean(),
           AIC1_bs[:, 1].mean(), AIC2_bs[:, 1].mean()],
    'd3': [AIC1[:, 2].mean(), AIC2[:, 2].mean(),
           AIC1_bs[:, 2].mean(), AIC2_bs[:, 2].mean()],
    'd4': [AIC1[:, 3].mean(), AIC2[:, 3].mean(), np.nan, np.nan],
}, index=['Linear', 'Logit', 'Linear spline', 'Logit spline'])

print("\\nAIC:")
print(aic_table.round(3).to_string())
aic_table.round(3).to_csv(output_dir / 'tab_AIC.txt', sep='|')

# BIC
bic_table = pd.DataFrame({
    'd1': [BIC1[:, 0].mean(), BIC2[:, 0].mean(),
           BIC1_bs[:, 0].mean(), BIC2_bs[:, 0].mean()],
    'd2': [BIC1[:, 1].mean(), BIC2[:, 1].mean(),
           BIC1_bs[:, 1].mean(), BIC2_bs[:, 1].mean()],
    'd3': [BIC1[:, 2].mean(), BIC2[:, 2].mean(),
           BIC1_bs[:, 2].mean(), BIC2_bs[:, 2].mean()],
    'd4': [BIC1[:, 3].mean(), BIC2[:, 3].mean(), np.nan, np.nan],
}, index=['Linear', 'Logit', 'Linear spline', 'Logit spline'])

print("\\nBIC:")
print(bic_table.round(3).to_string())
bic_table.round(3).to_csv(output_dir / 'tab_BIC.txt', sep='|')
"""))

# ============================================================
# Cell: 1st Stage Models (full data)
# ============================================================
cells.append(nbf.v4.new_markdown_cell("#### ステップ1のモデル推定（全データ）"))

cells.append(nbf.v4.new_code_cell("""\
# 1st stage specification 1: linear degree 1
X_lin = sm.add_constant(create_poly_features(df, 1))

m1st_lin_ANA = sm.OLS(df['y_ANA'].values, X_lin).fit()
m1st_lin_JAL = sm.OLS(df['y_JAL'].values, X_lin).fit()

# 1st stage specification 2: logit with b-spline degree 2
X_bs2 = sm.add_constant(create_bspline_features(df, 2))

m1st_logit_spline_ANA = sm.Logit(df['y_ANA'].values, X_bs2).fit(disp=0, method='bfgs', maxiter=200)
m1st_logit_spline_JAL = sm.Logit(df['y_JAL'].values, X_bs2).fit(disp=0, method='bfgs', maxiter=200)

# 1st stage specification 3: logit degree 1
m1st_logit_ANA = sm.Logit(df['y_ANA'].values, X_lin).fit(disp=0)
m1st_logit_JAL = sm.Logit(df['y_JAL'].values, X_lin).fit(disp=0)

# ステップ1の結果表示
print("ステップ1の推定完了")
print(f"  Linear ANA Adj. R-squared: {m1st_lin_ANA.rsquared_adj:.4f}")
print(f"  Linear JAL Adj. R-squared: {m1st_lin_JAL.rsquared_adj:.4f}")
print(f"  Logit ANA Log-Likelihood: {m1st_logit_ANA.llf:.4f}")
print(f"  Logit JAL Log-Likelihood: {m1st_logit_JAL.llf:.4f}")
print(f"  Logit+Spline ANA Log-Likelihood: {m1st_logit_spline_ANA.llf:.4f}")
print(f"  Logit+Spline JAL Log-Likelihood: {m1st_logit_spline_JAL.llf:.4f}")
"""))

# ============================================================
# Cell: Predicted Probabilities
# ============================================================
cells.append(nbf.v4.new_markdown_cell("#### 各モデルの予測確率"))

cells.append(nbf.v4.new_code_cell("""\
# 各モデルで予測値を計算

# linear: 0-1範囲外のものをクリップ (Rのロジックに従い、0-1内の最大最小でクリップ)
p_lin_ANA = m1st_lin_ANA.predict(X_lin)
p_lin_ANA_valid = p_lin_ANA[(p_lin_ANA >= 0) & (p_lin_ANA <= 1)]
p_lin_ANA = np.clip(p_lin_ANA, p_lin_ANA_valid.min(), p_lin_ANA_valid.max())

p_lin_JAL = m1st_lin_JAL.predict(X_lin)
p_lin_JAL_valid = p_lin_JAL[(p_lin_JAL >= 0) & (p_lin_JAL <= 1)]
p_lin_JAL = np.clip(p_lin_JAL, p_lin_JAL_valid.min(), p_lin_JAL_valid.max())

# logit with b-spline
p_logit_sp_ANA = m1st_logit_spline_ANA.predict(X_bs2)
p_logit_sp_JAL = m1st_logit_spline_JAL.predict(X_bs2)

# logit
p_logit_ANA = m1st_logit_ANA.predict(X_lin)
p_logit_JAL = m1st_logit_JAL.predict(X_lin)

df['p_lin_ANA'] = p_lin_ANA
df['p_lin_JAL'] = p_lin_JAL
df['p_logit_sp_ANA'] = p_logit_sp_ANA
df['p_logit_sp_JAL'] = p_logit_sp_JAL
df['p_logit_ANA'] = p_logit_ANA
df['p_logit_JAL'] = p_logit_JAL

print("予測確率の計算完了")
print(f"  p_lin_ANA: [{p_lin_ANA.min():.4f}, {p_lin_ANA.max():.4f}]")
print(f"  p_lin_JAL: [{p_lin_JAL.min():.4f}, {p_lin_JAL.max():.4f}]")
print(f"  p_logit_ANA: [{p_logit_ANA.min():.4f}, {p_logit_ANA.max():.4f}]")
print(f"  p_logit_JAL: [{p_logit_JAL.min():.4f}, {p_logit_JAL.max():.4f}]")
print(f"  p_logit_sp_ANA: [{p_logit_sp_ANA.min():.4f}, {p_logit_sp_ANA.max():.4f}]")
print(f"  p_logit_sp_JAL: [{p_logit_sp_JAL.min():.4f}, {p_logit_sp_JAL.max():.4f}]")
"""))

# ============================================================
# Cell: Long format
# ============================================================
cells.append(nbf.v4.new_markdown_cell("#### long形式に整理"))

cells.append(nbf.v4.new_code_cell("""\
# ANA行とJAL行に分解してlong形式にする
df_ANA = df.copy()
df_ANA['comp'] = 'ANA'
df_ANA['lin'] = df_ANA['p_lin_ANA']
df_ANA['logit'] = df_ANA['p_logit_ANA']
df_ANA['logit_sp'] = df_ANA['p_logit_sp_ANA']
df_ANA['p_lin_opp'] = df_ANA['p_lin_JAL']
df_ANA['p_logit_opp'] = df_ANA['p_logit_JAL']
df_ANA['p_logit_sp_opp'] = df_ANA['p_logit_sp_JAL']

df_JAL = df.copy()
df_JAL['comp'] = 'JAL'
df_JAL['lin'] = df_JAL['p_lin_JAL']
df_JAL['logit'] = df_JAL['p_logit_JAL']
df_JAL['logit_sp'] = df_JAL['p_logit_sp_JAL']
df_JAL['p_lin_opp'] = df_JAL['p_lin_ANA']
df_JAL['p_logit_opp'] = df_JAL['p_logit_ANA']
df_JAL['p_logit_sp_opp'] = df_JAL['p_logit_sp_ANA']

df_long = pd.concat([df_ANA, df_JAL], ignore_index=True)

print(f"df_long サイズ: {df_long.shape}")
df_long.head()
"""))

# ============================================================
# Cell: Step 1 - Results Table (linear & logit)
# ============================================================
cells.append(nbf.v4.new_markdown_cell("#### ステップ1の結果 (linear & logit)"))

cells.append(nbf.v4.new_code_cell("""\
# ステップ1の推定結果（linear & logit）の表示
# R: texreg::screenreg(list(m1st_lin_ANA, m1st_logit_ANA, m1st_lin_JAL, m1st_logit_JAL))
coef_names_1st = ['Intercept', 'Distance', 'Population', 'Flight_ANA',
                  'Flight_JAL', 'Train', 'Train:Distance', 'Train:Population',
                  'Train:Flight_ANA', 'Train:Flight_JAL']

results_1st = pd.DataFrame(index=coef_names_1st)
for label, model in [('ANA_linear', m1st_lin_ANA), ('ANA_logit', m1st_logit_ANA),
                     ('JAL_linear', m1st_lin_JAL), ('JAL_logit', m1st_logit_JAL)]:
    n_coef = min(len(model.params), len(coef_names_1st))
    results_1st[f'{label}_est'] = pd.Series(model.params[:n_coef], index=coef_names_1st[:n_coef])
    results_1st[f'{label}_se'] = pd.Series(model.bse[:n_coef], index=coef_names_1st[:n_coef])

print("ステップ1の結果 (linear & logit):")
for label in ['ANA_linear', 'ANA_logit', 'JAL_linear', 'JAL_logit']:
    print(f"\\n{label}:")
    for var in coef_names_1st:
        e = results_1st.loc[var, f'{label}_est']
        s = results_1st.loc[var, f'{label}_se']
        if pd.notna(e):
            print(f"  {var:20s}: {e:8.4f} ({s:.4f})")

# 保存
results_1st.round(4).to_csv(output_dir / 'tab_result_linear_logit.txt', sep='\\t')
"""))

# ============================================================
# Cell: Step 1 Visualization - Flight vs Probability
# ============================================================
cells.append(nbf.v4.new_markdown_cell("#### ステップ1の可視化: 便数と確率の関係"))

cells.append(nbf.v4.new_code_cell("""\
# 便数と確率の関係を可視化
# R: expand_grid で平均値ベースのデータに Flight_ANA/JAL を変えて予測

def predict_linear(model, X_design):
    \"\"\"線形モデルの予測値をクリップ\"\"\"
    pred = model.predict(X_design)
    return np.clip(pred, 0, 1)


# Flight_ANA を変化させた場合
mean_vals = df.mean(numeric_only=True)
flight_ANA_range = np.arange(df['Flight_ANA'].min(), df['Flight_ANA'].max() + 0.1, 0.1)

rows_fANA = []
for fa in flight_ANA_range:
    row = mean_vals.copy()
    row['Flight_ANA'] = fa
    rows_fANA.append(row)
df_pred_fANA = pd.DataFrame(rows_fANA)

X_pred_lin_fANA = sm.add_constant(create_poly_features(df_pred_fANA, 1))
X_pred_bs2_fANA = sm.add_constant(create_bspline_features(df_pred_fANA, 2, df_full=df))

pred_fANA = pd.DataFrame({
    'Flight': df_pred_fANA['Flight_ANA'].values,
    'pred_lin_ANA': predict_linear(m1st_lin_ANA, X_pred_lin_fANA),
    'pred_logit_ANA': m1st_logit_ANA.predict(X_pred_lin_fANA),
    'pred_logit_sp_ANA': m1st_logit_spline_ANA.predict(X_pred_bs2_fANA),
    'pred_lin_JAL': predict_linear(m1st_lin_JAL, X_pred_lin_fANA),
    'pred_logit_JAL': m1st_logit_JAL.predict(X_pred_lin_fANA),
    'pred_logit_sp_JAL': m1st_logit_spline_JAL.predict(X_pred_bs2_fANA),
})

# Flight_JAL を変化させた場合
flight_JAL_range = np.arange(df['Flight_JAL'].min(), df['Flight_JAL'].max() + 0.1, 0.1)

rows_fJAL = []
for fj in flight_JAL_range:
    row = mean_vals.copy()
    row['Flight_JAL'] = fj
    rows_fJAL.append(row)
df_pred_fJAL = pd.DataFrame(rows_fJAL)

X_pred_lin_fJAL = sm.add_constant(create_poly_features(df_pred_fJAL, 1))
X_pred_bs2_fJAL = sm.add_constant(create_bspline_features(df_pred_fJAL, 2, df_full=df))

pred_fJAL = pd.DataFrame({
    'Flight': df_pred_fJAL['Flight_JAL'].values,
    'pred_lin_ANA': predict_linear(m1st_lin_ANA, X_pred_lin_fJAL),
    'pred_logit_ANA': m1st_logit_ANA.predict(X_pred_lin_fJAL),
    'pred_logit_sp_ANA': m1st_logit_spline_ANA.predict(X_pred_bs2_fJAL),
    'pred_lin_JAL': predict_linear(m1st_lin_JAL, X_pred_lin_fJAL),
    'pred_logit_JAL': m1st_logit_JAL.predict(X_pred_lin_fJAL),
    'pred_logit_sp_JAL': m1st_logit_spline_JAL.predict(X_pred_bs2_fJAL),
})

# プロット (2x2 facets: flight_comp x name_comp)
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

for col_idx, comp_label in enumerate(['ANA', 'JAL']):
    for row_idx, (flight_label, pred_df) in enumerate(
            [('Flight_ANA', pred_fANA), ('Flight_JAL', pred_fJAL)]):
        ax = axes[row_idx, col_idx]
        for model_name, style in [('pred_lin', '-'), ('pred_logit', '--'), ('pred_logit_sp', ':')]:
            col_name = f'{model_name}_{comp_label}'
            if col_name in pred_df.columns:
                ax.plot(pred_df['Flight'], pred_df[col_name], style,
                        label=model_name.replace('pred_', ''))
        ax.set_xlabel(flight_label)
        ax.set_ylabel('y')
        ax.set_title(f'{flight_label} -> y_{comp_label}')
        ax.legend(fontsize=8)

plt.tight_layout()
plt.savefig(output_dir / 'fig_flight_probability.png', dpi=150)
plt.show()
"""))

# ============================================================
# Cell: Step 1 Visualization - Predicted vs Actual
# ============================================================
cells.append(nbf.v4.new_markdown_cell("#### 予測値と実測の関係"))

cells.append(nbf.v4.new_code_cell("""\
# 予測値と実測の関係を可視化
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

for col_idx, comp in enumerate(['ANA', 'JAL']):
    for row_idx, flight_col in enumerate(['Flight_ANA', 'Flight_JAL']):
        ax = axes[row_idx, col_idx]
        flight_vals = df[flight_col].values

        # 実測値
        y_col = f'y_{comp}'
        ax.scatter(flight_vals, df[y_col].values, s=10, alpha=0.3,
                   color='grey', label=y_col, zorder=1)

        # 予測値
        for pred_col, color, label in [
                (f'p_lin_{comp}', 'blue', 'linear'),
                (f'p_logit_{comp}', 'green', 'logit'),
                (f'p_logit_sp_{comp}', 'orange', 'logit_sp')]:
            ax.scatter(flight_vals, df[pred_col].values, s=5, alpha=0.5,
                       color=color, label=label, zorder=2)

        ax.set_xlabel(flight_col)
        ax.set_ylabel('y')
        ax.set_title(f'{flight_col} -> y_{comp}')
        ax.legend(fontsize=7, markerscale=2)

plt.tight_layout()
plt.savefig(output_dir / 'fig_first_stage.png', dpi=150)
plt.show()
"""))

# ============================================================
# Cell: Step 2 - OLS (Pi regression)
# ============================================================
cells.append(nbf.v4.new_markdown_cell(
    "## 5. ステップ2: 構造パラメータの推定\n\n"
    "### 1. 重み付き最小二乗法（OLS）"
))

cells.append(nbf.v4.new_code_cell("""\
# Pi (log-odds変換) を作成
eps_clip = 1e-10
df_long['Pi_lin'] = (np.log(np.clip(df_long['lin'], eps_clip, 1 - eps_clip))
                     - np.log(1 - np.clip(df_long['lin'], eps_clip, 1 - eps_clip)))
df_long['Pi_logit'] = (np.log(np.clip(df_long['logit'], eps_clip, 1 - eps_clip))
                        - np.log(1 - np.clip(df_long['logit'], eps_clip, 1 - eps_clip)))
df_long['Pi_logit_sp'] = (np.log(np.clip(df_long['logit_sp'], eps_clip, 1 - eps_clip))
                           - np.log(1 - np.clip(df_long['logit_sp'], eps_clip, 1 - eps_clip)))

# Flight を企業に合わせる
df_long['Flight'] = np.where(df_long['comp'] == 'ANA', df_long['Flight_ANA'], df_long['Flight_JAL'])
df_long['Flight_opp'] = np.where(df_long['comp'] == 'ANA', df_long['Flight_JAL'], df_long['Flight_ANA'])

# 交差項
df_long['Distance_Train'] = df_long['Distance'] * df_long['Train']
df_long['Population_Train'] = df_long['Population'] * df_long['Train']
df_long['Pop_Square_Train'] = df_long['Pop_Square'] * df_long['Train']
df_long['Flight_Train'] = df_long['Flight'] * df_long['Train']

# y変数
df_long['y'] = np.where(df_long['comp'] == 'ANA', df_long['y_ANA'], df_long['y_JAL'])

# OLS説明変数
X_vars_base = ['Constant', 'Distance', 'Population', 'Pop_Square', 'Train', 'Flight',
               'Distance_Train', 'Population_Train', 'Pop_Square_Train', 'Flight_Train']

# OLS推定 (各ステップ1のspecificationごと)
X_ols_lin = df_long[X_vars_base + ['p_lin_opp']].values
est1_1 = sm.OLS(df_long['Pi_lin'].values, X_ols_lin).fit()

X_ols_logit = df_long[X_vars_base + ['p_logit_opp']].values
est1_2 = sm.OLS(df_long['Pi_logit'].values, X_ols_logit).fit()

X_ols_logit_sp = df_long[X_vars_base + ['p_logit_sp_opp']].values
est1_3 = sm.OLS(df_long['Pi_logit_sp'].values, X_ols_logit_sp).fit()

param_names = ['(Intercept)', 'Distance', 'Population', 'Pop_Square', 'Train', 'Flight',
               'Train:Distance', 'Train:Population', 'Train:Pop_Square', 'Train:Flight', 'Delta']

print("OLS推定結果:")
ols_results = pd.DataFrame({
    'linear': est1_1.params,
    'logit': est1_2.params,
    'logit_sp': est1_3.params,
}, index=param_names)
print(ols_results.round(4).to_string())

# 保存
ols_results.round(4).to_csv(output_dir / 'tab_result_linear_logit_spline.txt', sep='\\t')
"""))

# ============================================================
# Cell: Step 2 - GMM
# ============================================================
cells.append(nbf.v4.new_markdown_cell("### 2. GMM推定"))

cells.append(nbf.v4.new_code_cell("""\
def calculate_moment(par, X, Z, y):
    \"\"\"モーメント条件: g(beta) = Z'(y - Lambda(X*beta)) を計算する。\"\"\"
    p = 1 / (1 + np.exp(-X @ par))
    g = Z.T @ (y - p)
    return g


def calculate_standard_error_gmm(par, X, Z, y):
    \"\"\"GMM推定量の標準誤差（サンドイッチ推定量）。\"\"\"
    p = 1 / (1 + np.exp(-X @ par))
    N = X.shape[0]

    # Omega = (1/N) * sum_i g_i g_i'
    residuals = (p - y).reshape(-1, 1)
    g_i = Z * residuals  # N x K_z
    Omega = g_i.T @ g_i / N

    # G = (1/N) * Z' * diag(p*(1-p)) * X
    weights = (p * (1 - p)).reshape(-1, 1)
    G = Z.T @ (X * weights) / N

    # V = (1/N) * G^{-1} Omega (G^{-1})'
    G_inv = np.linalg.solve(G, np.eye(G.shape[0]))
    V = G_inv @ Omega @ G_inv.T / N

    return np.sqrt(np.diag(V))


# 操作変数行列 (Flight_opp を使用)
Z = df_long[X_vars_base + ['Flight_opp']].values
X_lin_mat = df_long[X_vars_base + ['p_lin_opp']].values
X_logit_mat = df_long[X_vars_base + ['p_logit_opp']].values
X_logit_sp_mat = df_long[X_vars_base + ['p_logit_sp_opp']].values
y_vec = df_long['y'].values

param0 = np.zeros(X_lin_mat.shape[1])

# GMM (not IV: Z = X)
est2_1_gmm = fsolve(lambda x: calculate_moment(x, X_lin_mat, X_lin_mat, y_vec),
                     param0, full_output=True)
est2_2_gmm = fsolve(lambda x: calculate_moment(x, X_logit_mat, X_logit_mat, y_vec),
                     param0, full_output=True)
est2_3_gmm = fsolve(lambda x: calculate_moment(x, X_logit_sp_mat, X_logit_sp_mat, y_vec),
                     param0, full_output=True)

# GMM (IV: Z = instrument matrix)
est2_1_gmm2 = fsolve(lambda x: calculate_moment(x, X_lin_mat, Z, y_vec),
                      param0, full_output=True)
est2_2_gmm2 = fsolve(lambda x: calculate_moment(x, X_logit_mat, Z, y_vec),
                      param0, full_output=True)
est2_3_gmm2 = fsolve(lambda x: calculate_moment(x, X_logit_sp_mat, Z, y_vec),
                      param0, full_output=True)

print("GMM推定結果 (non-IV, logit 1st stage):")
for name, val in zip(param_names, est2_2_gmm[0]):
    print(f"  {name:20s}: {val:.4f}")

print("\\nGMM推定結果 (IV, logit 1st stage):")
for name, val in zip(param_names, est2_2_gmm2[0]):
    print(f"  {name:20s}: {val:.4f}")
"""))

# ============================================================
# Cell: Step 2 - MLE
# ============================================================
cells.append(nbf.v4.new_markdown_cell("### 3. 最尤推定 (MLE)")  )

cells.append(nbf.v4.new_code_cell("""\
# MLE: glm(y ~ ..., family = binomial)
# R の est3_1, est3_2, est3_3 に対応
X_mle_lin = df_long[X_vars_base + ['p_lin_opp']].values
X_mle_logit = df_long[X_vars_base + ['p_logit_opp']].values
X_mle_logit_sp = df_long[X_vars_base + ['p_logit_sp_opp']].values

est3_1 = sm.Logit(y_vec, X_mle_lin).fit(disp=0)
est3_2 = sm.Logit(y_vec, X_mle_logit).fit(disp=0)
est3_3 = sm.Logit(y_vec, X_mle_logit_sp).fit(disp=0)

print("MLE推定結果:")
mle_results = pd.DataFrame({
    'linear': est3_1.params,
    'logit': est3_2.params,
    'logit_sp': est3_3.params,
}, index=param_names)
print(mle_results.round(4).to_string())

# 保存
mle_results.round(4).to_csv(output_dir / 'tab_result_MLE.txt', sep='\\t')
"""))

# ============================================================
# Cell: OLS vs MLE Comparison
# ============================================================
cells.append(nbf.v4.new_markdown_cell("### OLS vs MLE の比較"))

cells.append(nbf.v4.new_code_cell("""\
# OLS (y=Pi) と MLE (y=action) の比較表
compare_df = pd.DataFrame(index=param_names)
compare_df['OLS_linear_est'] = est1_1.params
compare_df['OLS_linear_se'] = est1_1.bse
compare_df['OLS_logit_est'] = est1_2.params
compare_df['OLS_logit_se'] = est1_2.bse
compare_df['MLE_linear_est'] = est3_1.params
compare_df['MLE_linear_se'] = est3_1.bse
compare_df['MLE_logit_est'] = est3_2.params
compare_df['MLE_logit_se'] = est3_2.bse

print("OLS (y=Pi) vs MLE (y=action) の比較:")
for method in ['OLS_linear', 'OLS_logit', 'MLE_linear', 'MLE_logit']:
    print(f"\\n{method}:")
    for var in param_names:
        e = compare_df.loc[var, f'{method}_est']
        s = compare_df.loc[var, f'{method}_se']
        print(f"  {var:20s}: {e:8.4f} ({s:.4f})")

# 保存 (英語)
compare_df.round(4).to_csv(output_dir / 'tab_result_OLS_MLE_Eng.txt', sep='\\t')
"""))

# ============================================================
# Cell: Summary Table (Tab 9.2)
# ============================================================
cells.append(nbf.v4.new_markdown_cell("### 推定結果のまとめ (表9.2)"))

cells.append(nbf.v4.new_code_cell("""\
def calc_predict_probability(par, X):
    \"\"\"予測確率を計算する。\"\"\"
    p = 1 / (1 + np.exp(-X @ par))
    return p


def calc_log_likelihood(par, X, y):
    \"\"\"対数尤度を計算する。\"\"\"
    p = calc_predict_probability(par, X)
    p = np.clip(p, 1e-15, 1 - 1e-15)
    ll = y * np.log(p) + (1 - y) * np.log(1 - p)
    return np.sum(ll)


# 各メソッド・specificationの X 行列を構築
# OLS/MLE/GMM の対数尤度計算では、OLS推定量の場合は
# predict に使う X と対数尤度計算の X で p_opp 列が異なることに注意
# R の calculate_predict_probability は df_long のグローバル変数を参照

# GMM/IV-GMM の対数尤度はそれぞれの p_opp を使って計算
# (R: caluculate_log_likelihood(est2_1_gmm$x, f_spec = "p_lin_opp") etc.)

# 各方法ごとの対数尤度を計算するための X 行列マッピング
# f_spec -> X 行列
X_for_ll = {
    'p_lin_opp': X_lin_mat,
    'p_logit_opp': X_logit_mat,
    'p_logit_sp_opp': X_logit_sp_mat,
}

# 日本語変数名
jp_names = ['定数', '距離', '人口', '人口2乗', '新幹線', 'フライト',
            '新幹線:距離', '新幹線:人口', '新幹線:人口2乗', '新幹線:フライト', '競争']

results_list = []

# OLS
for label, est in [('OLS_linear', est1_1),
                    ('OLS_logit', est1_2),
                    ('OLS_logit_spline', est1_3)]:
    for name, val, se_val in zip(jp_names, est.params, est.bse):
        results_list.append({'var': name, 'est': val, 'se': se_val, 'method': label})

# MLE
for label, est in [('MLE_linear', est3_1),
                    ('MLE_logit', est3_2),
                    ('MLE_logit_spline', est3_3)]:
    for name, val, se_val in zip(jp_names, est.params, est.bse):
        results_list.append({'var': name, 'est': val, 'se': se_val, 'method': label})

# Moment (non-IV)
for label, est_res, X in [('Moment_linear', est2_1_gmm, X_lin_mat),
                           ('Moment_logit', est2_2_gmm, X_logit_mat),
                           ('Moment_logit_spline', est2_3_gmm, X_logit_sp_mat)]:
    se_vals = calculate_standard_error_gmm(est_res[0], X, X, y_vec)
    for name, val, se_val in zip(jp_names, est_res[0], se_vals):
        results_list.append({'var': name, 'est': val, 'se': se_val, 'method': label})

# IV Moment
for label, est_res, X in [('IV_Moment_linear', est2_1_gmm2, X_lin_mat),
                           ('IV_Moment_logit', est2_2_gmm2, X_logit_mat),
                           ('IV_Moment_logit_spline', est2_3_gmm2, X_logit_sp_mat)]:
    se_vals = calculate_standard_error_gmm(est_res[0], X, Z, y_vec)
    for name, val, se_val in zip(jp_names, est_res[0], se_vals):
        results_list.append({'var': name, 'est': val, 'se': se_val, 'method': label})

results_df = pd.DataFrame(results_list)

# 表の整形 (推定値と標準誤差を交互に表示)
rows_out = []
for name in jp_names:
    sub = results_df[results_df['var'] == name]
    row_est = {'var': name}
    row_se = {'var': ''}
    for _, r in sub.iterrows():
        row_est[r['method']] = f"{r['est']:.2f}"
        row_se[r['method']] = f" ({r['se']:.2f})"
    rows_out.append(row_est)
    rows_out.append(row_se)

# 追加行 (Num.obs, R², Log-Likelihood)
n_obs_ols = len(df_long)
summary_rows = []

# Num. obs.
nobs_row = {'var': 'Num. obs.'}
for label, est in [('OLS_linear', est1_1), ('OLS_logit', est1_2), ('OLS_logit_spline', est1_3)]:
    nobs_row[label] = str(int(est.nobs))
for label, est in [('MLE_linear', est3_1), ('MLE_logit', est3_2), ('MLE_logit_spline', est3_3)]:
    nobs_row[label] = str(int(est.nobs))
for label in ['Moment_linear', 'Moment_logit', 'Moment_logit_spline',
              'IV_Moment_linear', 'IV_Moment_logit', 'IV_Moment_logit_spline']:
    nobs_row[label] = str(n_obs_ols)
summary_rows.append(nobs_row)

# R²
rsq_row = {'var': 'R^2'}
for label, est in [('OLS_linear', est1_1), ('OLS_logit', est1_2), ('OLS_logit_spline', est1_3)]:
    rsq_row[label] = f"{est.rsquared:.2f}"
summary_rows.append(rsq_row)

# Log-Likelihood
ll_row = {'var': 'Log Likelihood'}
for label, params, X in [('OLS_linear', est1_1.params, X_lin_mat),
                           ('OLS_logit', est1_2.params, X_logit_mat),
                           ('OLS_logit_spline', est1_3.params, X_logit_sp_mat),
                           ('MLE_linear', est3_1.params, X_lin_mat),
                           ('MLE_logit', est3_2.params, X_logit_mat),
                           ('MLE_logit_spline', est3_3.params, X_logit_sp_mat)]:
    ll_val = calc_log_likelihood(params, X, y_vec)
    ll_row[label] = f"{ll_val:.2f}"

# GMM LL: R code uses calculate_predict_probability with f_spec
for label, est_res, f_spec in [
        ('Moment_linear', est2_1_gmm, 'p_lin_opp'),
        ('Moment_logit', est2_2_gmm, 'p_logit_opp'),
        ('Moment_logit_spline', est2_3_gmm, 'p_logit_sp_opp'),
        ('IV_Moment_linear', est2_1_gmm2, 'p_lin_opp'),
        ('IV_Moment_logit', est2_2_gmm2, 'p_logit_opp'),
        ('IV_Moment_logit_spline', est2_3_gmm2, 'p_logit_sp_opp')]:
    X_ll = X_for_ll[f_spec]
    ll_val = calc_log_likelihood(est_res[0], X_ll, y_vec)
    ll_row[label] = f"{ll_val:.2f}"
summary_rows.append(ll_row)

rows_out.extend(summary_rows)

tab9_2 = pd.DataFrame(rows_out)
# 列の順序を整理
col_order = ['var',
             'OLS_linear', 'OLS_logit', 'OLS_logit_spline',
             'MLE_linear', 'MLE_logit', 'MLE_logit_spline',
             'Moment_linear', 'Moment_logit', 'Moment_logit_spline',
             'IV_Moment_linear', 'IV_Moment_logit', 'IV_Moment_logit_spline']
col_order = [c for c in col_order if c in tab9_2.columns]
tab9_2 = tab9_2[col_order].fillna('')

print("表9.2: 推定結果のまとめ")
print(tab9_2.to_string(index=False))

# 保存
tab9_2.to_csv(output_dir / 'tab9_2_estimation.txt', sep='|', index=False)
print("\\n保存先: output/tab9_2_estimation.txt")
"""))

# ============================================================
# Cell: Equilibrium Calculation - Function Definitions
# ============================================================
cells.append(nbf.v4.new_markdown_cell("## 6. 均衡の計算"))

cells.append(nbf.v4.new_code_cell("""\
def calculate_best_response(p_opp, par, X):
    \"\"\"最適反応を計算する。X にp_oppを結合してlogit確率を返す。\"\"\"
    X_full = np.column_stack([X, p_opp])
    p = 1 / (1 + np.exp(-X_full @ par))
    return p


def calculate_best_response_opponent(p_opp, par, X):
    \"\"\"相手プレイヤーの最適反応を計算する。
    ANA部分とJAL部分を入れ替えて返す。\"\"\"
    n_market = X.shape[0] // 2
    p = calculate_best_response(p_opp, par, X)
    p_opp_new = np.concatenate([p[n_market:], p[:n_market]])
    return p_opp_new


def find_equilibrium(p_init, par, X, tol=1e-6, max_iter=1000):
    \"\"\"最適反応の反復計算で均衡を求める。\"\"\"
    n_market = X.shape[0] // 2
    p_opp = np.concatenate([p_init[n_market:], p_init[:n_market]])

    for i in range(max_iter):
        p_opp_new = calculate_best_response_opponent(p_opp, par, X)
        if np.sum((p_opp_new - p_opp)**2) < tol:
            break
        p_opp = p_opp_new

    # 均衡確率ベクトル
    p_eq = np.concatenate([p_opp[n_market:], p_opp[:n_market]])
    return p_eq


def find_equilibrium_ANA_JAL(ids, optimal, par, df_long_input):
    \"\"\"ANA/JALの2社参入ゲームの均衡確率を計算する。
    par: 11要素 [Intercept, Distance, Population, Pop_Square, Train, Flight,
                 Train:Distance, Train:Population, Train:Pop_Square, Train:Flight, Delta]
    \"\"\"
    df_sub = df_long_input[df_long_input['id'].isin(ids)].copy()

    X = df_sub[['Constant', 'Distance', 'Population', 'Pop_Square', 'Train', 'Flight',
                'Distance_Train', 'Population_Train', 'Flight_Train', 'Pop_Square_Train']].values

    n_market = X.shape[0] // 2

    if optimal == 'ANA':
        p_init = np.concatenate([np.ones(n_market), np.zeros(n_market)])
    else:
        p_init = np.concatenate([np.zeros(n_market), np.ones(n_market)])

    p_eq = find_equilibrium(p_init, par, X)
    return p_eq


print("均衡計算の関数定義完了")
"""))

# ============================================================
# Cell: Compute Equilibrium
# ============================================================
cells.append(nbf.v4.new_code_cell("""\
# comp列でソート (R: arrange(comp))
df_long = df_long.sort_values('comp').reset_index(drop=True)

# IV GMM logit の推定結果を利用 (R: est2_2_gmm2$x)
par_gmm = est2_2_gmm2[0]
ids_unique = df['id'].unique()

new_eq_ana = find_equilibrium_ANA_JAL(ids_unique, 'ANA', par_gmm, df_long)
new_eq_jal = find_equilibrium_ANA_JAL(ids_unique, 'JAL', par_gmm, df_long)

df_long['new_eq_ana'] = new_eq_ana
df_long['new_eq_jal'] = new_eq_jal
df_long['dif'] = np.abs(new_eq_ana - new_eq_jal)

print("均衡計算完了")
print(f"  均衡の差 (|ANA - JAL|) の平均: {df_long['dif'].mean():.6f}")
print(f"  均衡の差の最大値: {df_long['dif'].max():.6f}")
"""))

# ============================================================
# Cell: Equilibrium Multiplicity Plot
# ============================================================
cells.append(nbf.v4.new_markdown_cell("#### 複数均衡の確認"))

cells.append(nbf.v4.new_code_cell("""\
# 複数均衡の確認 (R: facet_wrap(~comp))
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
for idx, comp in enumerate(['ANA', 'JAL']):
    mask = df_long['comp'] == comp
    ax = axes[idx]
    sc = ax.scatter(df_long.loc[mask, 'new_eq_ana'],
                    df_long.loc[mask, 'new_eq_jal'],
                    c=df_long.loc[mask, 'dif'], cmap='RdBu_r',
                    vmin=0, vmax=max(0.02, df_long['dif'].max()),
                    s=10, alpha=0.6)
    ax.plot([0, 1], [0, 1], 'k-', linewidth=0.5)
    ax.set_xlabel('ANA optimal equilibrium')
    ax.set_ylabel('JAL optimal equilibrium')
    ax.set_title(comp)
    plt.colorbar(sc, ax=ax, label='|difference|')
plt.tight_layout()
plt.savefig(output_dir / 'fig_equilibrium_multiplicity.png', dpi=150)
plt.show()
"""))

# ============================================================
# Cell: Model Fit Check
# ============================================================
cells.append(nbf.v4.new_markdown_cell("#### ステップ1のモデルフィットの確認"))

cells.append(nbf.v4.new_code_cell("""\
# ステップ1のモデルフィットの確認
# R: p_2nd = calculate_predict_probability(coef(est1_2))
# est1_2 はOLS-logitの推定結果で、X は df_long の説明変数

# OLS-logit (est1_2) のパラメータで予測確率を計算
X_for_p2nd = df_long[X_vars_base + ['p_logit_opp']].values
p_2nd = calc_predict_probability(est1_2.params, X_for_p2nd)
df_long['p_2nd'] = p_2nd

# R の df_long_fit_plot に対応する整形
# wide形式: id ごとに ANA/JAL の new_eq_jal と p_lin_opp を取得
df_wide = df_long.pivot_table(index='id', columns='comp',
                               values=['new_eq_jal', 'p_lin_opp', 'p_2nd'],
                               aggfunc='first')

# ANA_eq = new_eq_jal from ANA row, JAL_eq = new_eq_jal from JAL row
# JAL_p_1st = p_lin_opp from ANA row (opponent), ANA_p_1st = p_lin_opp from JAL row

fig, axes = plt.subplots(2, 2, figsize=(12, 10))

for col_idx, comp in enumerate(['ANA', 'JAL']):
    eq_vals = df_wide[('new_eq_jal', comp)].values

    # p_1st: 相手のp_lin_opp -> 自分の1st stage予測
    # ANA の場合: p_1st = p_lin_opp from JAL row (= p_lin_ANA)
    p_1st = df_wide[('p_lin_opp', 'JAL' if comp == 'ANA' else 'ANA')].values
    p_2nd_vals = df_wide[('p_2nd', comp)].values

    for row_idx, (p_vals, p_label) in enumerate(
            [(p_1st, 'p_1st'), (p_2nd_vals, 'p_2nd')]):
        ax = axes[row_idx, col_idx]
        color_vals = np.abs(p_vals - eq_vals)
        sc = ax.scatter(eq_vals, p_vals, c=color_vals, cmap='RdBu_r',
                        s=10, alpha=0.6)
        ax.plot([0, 1], [0, 1], 'k-', linewidth=0.5)
        ax.set_xlabel('Equilibrium probability')
        ax.set_ylabel(p_label)
        ax.set_title(f'{p_label} / {comp}')
        plt.colorbar(sc, ax=ax, label='|deviation|')

plt.tight_layout()
plt.savefig(output_dir / 'fig_model_fit.png', dpi=150)
plt.show()
"""))

# ============================================================
# Cell: Counterfactual - Hokuriku Shinkansen
# ============================================================
cells.append(nbf.v4.new_markdown_cell(
    "## 7. 反実仮想分析: 北陸新幹線\n\n"
    "北陸新幹線が開通した場合の航空路線への参入確率の変化を分析する。"
))

cells.append(nbf.v4.new_code_cell("""\
# 路線情報をマージ
df_long = df_long.merge(id_port_dyad_long, on='id', how='left')

# 北陸新幹線に関連する空港
hokuriku_ports = ['小松', '富山', '松本', '羽田', '伊丹']

df_long_hokuriku = df_long.copy()
df_long_hokuriku['hokuriku'] = (
    df_long_hokuriku['port1'].isin(hokuriku_ports) &
    df_long_hokuriku['port2'].isin(hokuriku_ports)
).astype(int)

# 反実仮想: 北陸新幹線が通る路線のTrain=1にする
df_long_hokuriku['Train_orig'] = df_long_hokuriku['Train']
df_long_hokuriku['Train'] = np.minimum(1, df_long_hokuriku['Train'] + df_long_hokuriku['hokuriku'])

# 交差項を更新
df_long_hokuriku['Distance_Train'] = df_long_hokuriku['Distance'] * df_long_hokuriku['Train']
df_long_hokuriku['Population_Train'] = df_long_hokuriku['Population'] * df_long_hokuriku['Train']
df_long_hokuriku['Pop_Square_Train'] = df_long_hokuriku['Pop_Square'] * df_long_hokuriku['Train']
df_long_hokuriku['Flight_Train'] = df_long_hokuriku['Flight'] * df_long_hokuriku['Train']

print(f"北陸新幹線関連路線数 (ANA+JALの行数): {df_long_hokuriku['hokuriku'].sum()}")
"""))

# ============================================================
# Cell: Counterfactual Simulation
# ============================================================
cells.append(nbf.v4.new_code_cell("""\
# OLS-logitの推定結果（est1_2）を利用してシミュレーション
# R: est1_2$coefficients
par_ols_logit = est1_2.params

# 元の均衡 (df_long を使って計算)
orig_eq_ana = find_equilibrium_ANA_JAL(ids_unique, 'ANA', par_ols_logit, df_long)
orig_eq_jal = find_equilibrium_ANA_JAL(ids_unique, 'JAL', par_ols_logit, df_long)

# 北陸新幹線後の均衡 (df_long_hokuriku を使って計算)
hoku_eq_ana = find_equilibrium_ANA_JAL(ids_unique, 'ANA', par_ols_logit, df_long_hokuriku)
hoku_eq_jal = find_equilibrium_ANA_JAL(ids_unique, 'JAL', par_ols_logit, df_long_hokuriku)

df_long_hokuriku['orig_eq_ana'] = orig_eq_ana
df_long_hokuriku['orig_eq_jal'] = orig_eq_jal
df_long_hokuriku['hoku_eq_ana'] = hoku_eq_ana
df_long_hokuriku['hoku_eq_jal'] = hoku_eq_jal

print("反実仮想分析の均衡計算完了")
print(f"  元の均衡 (ANA): [{orig_eq_ana.min():.4f}, {orig_eq_ana.max():.4f}]")
print(f"  北陸新幹線後 (ANA): [{hoku_eq_ana.min():.4f}, {hoku_eq_ana.max():.4f}]")
"""))

# ============================================================
# Cell: Heatmap
# ============================================================
cells.append(nbf.v4.new_markdown_cell("#### ヒートマップの作成"))

cells.append(nbf.v4.new_code_cell("""\
# port_order を作成 (R: port_order)
p2_order = id_port_dyad_long[id_port_dyad_long['port2'] != '静岡']['port2'].unique().tolist()

south_tail_df = id_port_dyad_long[id_port_dyad_long['port1'] != '静岡']
south_tail = [p for p in south_tail_df['port1'].unique() if p not in p2_order]

port_order = list(dict.fromkeys(p2_order + south_tail))

# 北陸新幹線関連路線のヒートマップ (4パネル: counterfactual x comp)
cf_hoku = df_long_hokuriku[df_long_hokuriku['hokuriku'] == 1].copy()

if len(cf_hoku) > 0:
    # 港のリスト (北陸関連のみ)
    hoku_ports_in_data = sorted(
        set(cf_hoku['port1'].tolist() + cf_hoku['port2'].tolist()),
        key=lambda x: port_order.index(x) if x in port_order else len(port_order))

    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    for row_idx, (eq_col, cf_label) in enumerate([
            ('orig_eq_ana', '北陸新幹線なし'),
            ('hoku_eq_ana', '北陸新幹線あり')]):
        for col_idx, comp in enumerate(['ANA', 'JAL']):
            ax = axes[row_idx, col_idx]
            cf_comp = cf_hoku[cf_hoku['comp'] == comp].copy()

            # ピボットテーブルを作成
            pivot_data = pd.DataFrame(
                np.nan, index=hoku_ports_in_data, columns=hoku_ports_in_data)

            for _, row in cf_comp.iterrows():
                p1, p2 = row['port1'], row['port2']
                val = row[eq_col]
                if p1 in hoku_ports_in_data and p2 in hoku_ports_in_data:
                    # 両方向にセット（dyad_long の構造に合わせて）
                    pivot_data.loc[p1, p2] = val
                    pivot_data.loc[p2, p1] = val
                    # dyad名からも逆方向を取得
                    dyad = row.get('dyad', '')
                    if '_' in str(dyad):
                        parts = dyad.split('_')
                        pivot_data.loc[parts[1], parts[0]] = val
                        pivot_data.loc[parts[0], parts[1]] = val

            # 下三角のみ表示 (R: eq = ifelse(port1_id >= port2_id, eq, NA))
            n_ports = len(hoku_ports_in_data)
            for i in range(n_ports):
                for j in range(i + 1, n_ports):
                    pivot_data.iloc[j, i] = np.nan

            im = ax.imshow(pivot_data.values.astype(float),
                           cmap='RdBu_r', vmin=0, vmax=1, aspect='auto')
            ax.set_xticks(range(n_ports))
            ax.set_xticklabels(hoku_ports_in_data, rotation=45, ha='right', fontsize=8)
            ax.set_yticks(range(n_ports))
            ax.set_yticklabels(hoku_ports_in_data, fontsize=8)
            ax.set_title(f'{cf_label} / {comp}')
            plt.colorbar(im, ax=ax, label='均衡参入確率', shrink=0.8)

    plt.tight_layout()
    plt.savefig(output_dir / 'fig_hokuriku_heatmap.png', dpi=150)
    plt.show()
else:
    print("北陸新幹線関連路線のデータが見つかりませんでした")
"""))

# ============================================================
# Cell: Counterfactual Results Table (Tab 9.3)
# ============================================================
cells.append(nbf.v4.new_markdown_cell("#### 反実仮想分析の結果 (表9.3)"))

cells.append(nbf.v4.new_code_cell("""\
# 北陸新幹線関連路線の結果 (Tab 9.3)
cf_result = df_long_hokuriku[df_long_hokuriku['hokuriku'] == 1].copy()

cf_table = cf_result[['dyad', 'comp', 'orig_eq_ana', 'hoku_eq_ana']].copy()
cf_table.columns = ['ルート', '企業', '北陸新幹線なし', '北陸新幹線あり']

# R: pivot_wider(names_from = "企業", values_from = c("北陸新幹線あり", "北陸新幹線なし"))
cf_pivot = cf_table.pivot_table(
    index='ルート',
    columns='企業',
    values=['北陸新幹線なし', '北陸新幹線あり'],
    aggfunc='first'
)

# 列の順序を整理 (R: relocate)
col_order_cf = [('北陸新幹線なし', 'ANA'), ('北陸新幹線あり', 'ANA'),
                ('北陸新幹線なし', 'JAL'), ('北陸新幹線あり', 'JAL')]
col_order_cf = [c for c in col_order_cf if c in cf_pivot.columns]
if col_order_cf:
    cf_pivot = cf_pivot[col_order_cf]

print("表9.3: 反実仮想分析の結果（北陸新幹線関連路線）")
print(cf_pivot.round(3).to_string())

# 保存
cf_pivot.round(3).to_csv(output_dir / 'tab9_3_counter_factual.txt', sep='|')
print("\\n保存先: output/tab9_3_counter_factual.txt")
"""))

# ============================================================
# Cell: Summary
# ============================================================
cells.append(nbf.v4.new_markdown_cell("## 8. まとめ"))

cells.append(nbf.v4.new_code_cell("""\
print("=" * 60)
print("第9章の分析完了")
print("=" * 60)
print("\\n出力ファイル:")
for f in sorted(output_dir.glob('*')):
    print(f"  {f.name}")
"""))

nb.cells = cells
nbf.write(nb, 'main_ch09.ipynb')
print("Generated: main_ch09.ipynb")
