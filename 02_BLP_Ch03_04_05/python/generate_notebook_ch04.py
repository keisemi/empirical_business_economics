"""Generate main_ch04.ipynb for Chapter 4: BLP Estimation."""
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
# Title
# ============================================================
cells.append(nbf.v4.new_markdown_cell(
    "# 第4章 ランダム係数ロジット (BLP) モデルの推定\n\n"
    "ロジットモデル・入れ子型ロジットモデル・ランダム係数ロジット (BLP) モデルの\n"
    "弾力性行列を比較し、BLPモデルによるプライシング応用を行う。"
))

# ============================================================
# Setup
# ============================================================
cells.append(nbf.v4.new_markdown_cell("## 1. Pythonに関する下準備"))

cells.append(nbf.v4.new_code_cell("""\
import numpy as np
import pandas as pd
from scipy import optimize
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
intermediate_dir = base_dir / 'intermediate'
output_dir = base_dir / 'output'
output_dir.mkdir(exist_ok=True)
intermediate_dir.mkdir(exist_ok=True)

print("Setup complete.")
"""))

# ============================================================
# 2SLS helper
# ============================================================
cells.append(nbf.v4.new_markdown_cell("### 2SLS推定関数の定義"))

cells.append(nbf.v4.new_code_cell("""\
def tsls_estimate(y, X_exog, X_endog, Z_instruments):
    \"\"\"
    Manual 2SLS estimation with HC1 robust standard errors.

    Parameters
    ----------
    y : (N,) array - dependent variable
    X_exog : (N, k1) array - exogenous regressors (including constant)
    X_endog : (N, k2) array - endogenous regressors
    Z_instruments : (N, m) array - excluded instruments

    Returns
    -------
    dict with keys: beta, se, resid, V
        beta order: [exog coefs, endog coefs]
    \"\"\"
    N = len(y)
    y = y.reshape(-1, 1) if y.ndim == 1 else y

    # Full instrument set: exogenous + excluded instruments
    Z_full = np.hstack([X_exog, Z_instruments])

    # First stage: project endogenous onto Z_full
    PZ = Z_full @ np.linalg.solve(Z_full.T @ Z_full, Z_full.T)
    X_endog_hat = PZ @ X_endog

    # Second stage regressors
    X_hat = np.hstack([X_exog, X_endog_hat])
    X = np.hstack([X_exog, X_endog])

    # 2SLS estimator: (X_hat' X)^{-1} X_hat' y
    beta = np.linalg.solve(X_hat.T @ X, X_hat.T @ y)

    # Residuals using original X
    resid = y - X @ beta

    # HC1 robust standard errors
    k = X.shape[1]
    bread = np.linalg.inv(X_hat.T @ X)
    meat = np.zeros((k, k))
    for i in range(N):
        xi = X_hat[i:i+1, :].T
        meat += (resid[i, 0] ** 2) * (xi @ xi.T)
    meat *= N / (N - k)

    V = bread @ meat @ bread.T
    se = np.sqrt(np.diag(V))

    return {
        'beta': beta.flatten(),
        'se': se,
        'resid': resid.flatten(),
        'V': V
    }

print("2SLS推定関数を定義しました。")
"""))

# ============================================================
# Data Loading
# ============================================================
cells.append(nbf.v4.new_markdown_cell("## 2. データの読み込み"))

cells.append(nbf.v4.new_code_cell("""\
# Ch03で作成したクリーニング済みのデータを読み込む
data = pd.read_csv(intermediate_dir / 'data_cleaned.csv')

# logit_shareを作成
data['logit_share'] = np.log(data['share']) - np.log(data['share0'])

print(f"データサイズ: {data.shape}")
print(f"年数: {data['year'].nunique()}")
print(f"車種数: {data['NameID'].nunique()}")
"""))

# ============================================================
# Part 1: Logit Elasticity Matrix (Tab 4.1)
# ============================================================
cells.append(nbf.v4.new_markdown_cell(
    "## 3. ロジットモデルにおける価格弾力性行列 (表 4.1)"
))

cells.append(nbf.v4.new_code_cell("""\
# Differentiation IVを用いたロジットモデルの推定 (manual 2SLS)
y_logit = data['logit_share'].values
X_exog_logit = np.column_stack([
    np.ones(len(data)),
    data[['hppw', 'FuelEfficiency', 'size']].values
])
X_endog_logit = data[['price']].values
Z_iv_GH_logit = data[['iv_GH_own_hppw', 'iv_GH_own_FuelEfficiency', 'iv_GH_own_size',
                        'iv_GH_other_hppw', 'iv_GH_other_FuelEfficiency',
                        'iv_GH_other_size']].values

res_logit_GH = tsls_estimate(y_logit, X_exog_logit, X_endog_logit, Z_iv_GH_logit)

# beta order: [const, hppw, FuelEfficiency, size, price]
logit_var_names = ['const', 'hppw', 'FuelEfficiency', 'size', 'price']
print("GH-IV Logit 推定結果:")
print(f"{'変数':<20} {'推定値':>10} {'標準誤差':>10}")
print("-" * 42)
for name, b, s in zip(logit_var_names, res_logit_GH['beta'], res_logit_GH['se']):
    print(f"{name:<20} {b:>10.4f} {s:>10.4f}")
"""))

cells.append(nbf.v4.new_code_cell("""\
# 日評自動車の2016年データで弾力性行列を作成
data_NIPPYO = data[data['Nippyo'] == 1].copy()

dt2016 = data_NIPPYO[data_NIPPYO['year'] == 2016][['price', 'share', 'NameID', 'Name']].copy()
dt2016 = dt2016.sort_values('NameID').reset_index(drop=True)

price_vec = dt2016['price'].values
share_vec = dt2016['share'].values
nameID_vec = dt2016['NameID'].values
J = len(price_vec)

alpha_logit = res_logit_GH['beta'][4]  # price coefficient

# 自己弾力性: alpha * p_j * (1 - s_j)
own_elas = alpha_logit * price_vec * (1 - share_vec)

# 交差弾力性: -alpha * p_k * s_k (列kの価格変化が行jに与える影響)
cross_elas = (-1) * alpha_logit * price_vec * share_vec

# 弾力性行列: 各列kについて cross_elas[k] を全行に入れる
elas_mat = np.tile(cross_elas, (J, 1))  # (J, J) 行列
np.fill_diagonal(elas_mat, own_elas)

# 4車種の抽出
betard = dt2016.loc[dt2016['Name'] == 'アルファード', 'NameID'].values[0]
sedan = dt2016.loc[dt2016['Name'] == 'カローラ', 'NameID'].values[0]
suv = dt2016.loc[dt2016['Name'] == 'ジューク', 'NameID'].values[0]
kei = dt2016.loc[dt2016['Name'] == 'タント', 'NameID'].values[0]

target_ids = [betard, sedan, suv, kei]
target_names = ['ベータード', 'セダン(A)', 'SUV(B)', '軽自動車(C)']

# 行・列のインデックスを取得
idx = [np.where(nameID_vec == tid)[0][0] for tid in target_ids]

elas_mat_restricted = elas_mat[np.ix_(idx, idx)]
elas_df = pd.DataFrame(elas_mat_restricted, index=target_names, columns=target_names)

print("表 4.1: ロジットモデルにおける弾力性行列 (2016年)")
print(elas_df.round(4).to_string())
elas_df.round(4).to_csv(output_dir / 'tab4_1_elas_mat_restricted.txt', sep='\\t')
"""))

# ============================================================
# Part 2: Nested Logit (Tab 4.2, 4.3)
# ============================================================
cells.append(nbf.v4.new_markdown_cell(
    "## 4. 入れ子型ロジットモデルの推定 (表 4.2, 4.3)"
))

cells.append(nbf.v4.new_code_cell("""\
# 入れ子型ロジット推定のためのBLP操作変数の作成
# グループ = year + Maker + Type で定義

# マーケット・企業・Type レベルにおける各製品属性の和と二乗和
for var in ['hppw', 'FuelEfficiency', 'size']:
    data[f'{var}_sum_own_nest'] = data.groupby(['year', 'Maker', 'Type'])[var].transform('sum')
    data[f'{var}_sqr_sum_own_nest'] = data.groupby(['year', 'Maker', 'Type'])[var].transform(
        lambda x: (x**2).sum()
    )

data['group_n_nest'] = data.groupby(['year', 'Maker', 'Type'])['Sales'].transform('count')

# マーケット (year + Type) レベルでの各製品属性の和と二乗和
for var in ['hppw', 'FuelEfficiency', 'size']:
    data[f'{var}_sum_mkt_nest'] = data.groupby(['year', 'Type'])[var].transform('sum')
    data[f'{var}_sqr_sum_mkt_nest'] = data.groupby(['year', 'Type'])[var].transform(
        lambda x: (x**2).sum()
    )

data['mkt_n_nest'] = data.groupby(['year', 'Type'])['Sales'].transform('count')

# BLP操作変数 (nested)
for var in ['hppw', 'FuelEfficiency', 'size']:
    data[f'iv_BLP_own_{var}_nest'] = data[f'{var}_sum_own_nest'] - data[var]
    data[f'iv_BLP_other_{var}_nest'] = data[f'{var}_sum_mkt_nest'] - data[f'{var}_sum_own_nest']

data['iv_BLP_own_num_nest'] = data['group_n_nest'] - 1
data['iv_BLP_other_num_nest'] = data['mkt_n_nest'] - data['group_n_nest']

# Differentiation IV (nested)
for var in ['hppw', 'FuelEfficiency', 'size']:
    data[f'iv_GH_own_{var}_nest'] = (
        (data['group_n_nest'] - 1) * data[var]**2 +
        (data[f'{var}_sqr_sum_own_nest'] - data[var]**2) -
        2 * data[var] * (data[f'{var}_sum_own_nest'] - data[var])
    )
    data[f'iv_GH_other_{var}_nest'] = (
        (data['mkt_n_nest'] - data['group_n_nest']) * data[var]**2 +
        (data[f'{var}_sqr_sum_mkt_nest'] - data[f'{var}_sqr_sum_own_nest']) -
        2 * data[var] * (data[f'{var}_sum_mkt_nest'] - data[f'{var}_sum_own_nest'])
    )

# 不要な列を削除
drop_cols = [c for c in data.columns if c.endswith('_sum_own_nest') or c.endswith('_sum_mkt_nest')
             or c.endswith('_sqr_sum_own_nest') or c.endswith('_sqr_sum_mkt_nest')]
data = data.drop(columns=drop_cols + ['mkt_n_nest', 'group_n_nest'])

print("入れ子型ロジットの操作変数を構築しました。")
print(f"nested IV列: {[c for c in data.columns if 'nest' in c]}")
"""))

cells.append(nbf.v4.new_code_cell("""\
# インサイドシェアを定義
data['sum_year_body'] = data.groupby(['year', 'Type'])['Sales'].transform('sum')
data['inside_share'] = data['Sales'] / data['sum_year_body']
data['log_inside_share'] = np.log(data['inside_share'])
data = data.drop(columns=['sum_year_body'])

# logit_share を再計算（念のため更新）
data['logit_share'] = np.log(data['share']) - np.log(data['share0'])

# データの保存 (第5章で使用)
data.to_csv(intermediate_dir / 'data_for_estimation.csv', index=False)
print("data_for_estimation.csv を保存しました。")
print(f"データ形状: {data.shape}")
"""))

cells.append(nbf.v4.new_code_cell("""\
# --- OLS ---
y_nest = data['logit_share'].values
X_ols_nest = np.column_stack([
    np.ones(len(data)),
    data['price'].values,
    data['log_inside_share'].values,
    data['hppw'].values,
    data['FuelEfficiency'].values,
    data['size'].values
])
# OLS order: [const, price, log_inside_share, hppw, FuelEfficiency, size]
ols_nest_names = ['const', 'price', 'log_inside_share', 'hppw', 'FuelEfficiency', 'size']

beta_ols_nest = np.linalg.lstsq(X_ols_nest, y_nest, rcond=None)[0]
resid_ols_nest = y_nest - X_ols_nest @ beta_ols_nest
N_ols = len(y_nest)
k_ols = X_ols_nest.shape[1]

# HC1 robust SE for OLS
meat_ols = np.zeros((k_ols, k_ols))
for i in range(N_ols):
    xi = X_ols_nest[i:i+1, :].T
    meat_ols += (resid_ols_nest[i] ** 2) * (xi @ xi.T)
meat_ols *= N_ols / (N_ols - k_ols)
bread_ols = np.linalg.inv(X_ols_nest.T @ X_ols_nest)
V_ols = bread_ols @ meat_ols @ bread_ols
se_ols_nest = np.sqrt(np.diag(V_ols))

print("OLS推定結果:")
print(f"{'変数':<25} {'推定値':>10} {'標準誤差':>10}")
print("-" * 47)
for name, b, s in zip(ols_nest_names, beta_ols_nest, se_ols_nest):
    print(f"{name:<25} {b:>10.4f} {s:>10.4f}")

# --- IV推定: BLP IV (入れ子型) ---
# 内生変数: price, log_inside_share
# 外生変数: const, hppw, FuelEfficiency, size
X_exog_nest = np.column_stack([
    np.ones(len(data)),
    data[['hppw', 'FuelEfficiency', 'size']].values
])
X_endog_nest = data[['price', 'log_inside_share']].values

Z_iv_BLP_nest = data[['iv_BLP_own_hppw_nest', 'iv_BLP_own_FuelEfficiency_nest',
                        'iv_BLP_own_size_nest', 'iv_BLP_other_hppw_nest',
                        'iv_BLP_other_FuelEfficiency_nest', 'iv_BLP_other_size_nest',
                        'iv_BLP_own_num_nest', 'iv_BLP_other_num_nest']].values

res_nest_blp = tsls_estimate(y_nest, X_exog_nest, X_endog_nest, Z_iv_BLP_nest)
# beta order: [const, hppw, FuelEfficiency, size, price, log_inside_share]
nest_iv_names = ['const', 'hppw', 'FuelEfficiency', 'size', 'price', 'log_inside_share']

print("\\nIV (BLP IV) 推定結果:")
print(f"{'変数':<25} {'推定値':>10} {'標準誤差':>10}")
print("-" * 47)
for name, b, s in zip(nest_iv_names, res_nest_blp['beta'], res_nest_blp['se']):
    print(f"{name:<25} {b:>10.4f} {s:>10.4f}")

# 推定結果の表示 (R版ではOLSとIV_BLPの2列のみ表示)
print("\\n表 4.2: 入れ子型ロジットモデルの推定結果 (まとめ)")
# OLS beta dict
ols_dict = dict(zip(ols_nest_names, beta_ols_nest))
ols_se_dict = dict(zip(ols_nest_names, se_ols_nest))
iv_dict = dict(zip(nest_iv_names, res_nest_blp['beta']))
iv_se_dict = dict(zip(nest_iv_names, res_nest_blp['se']))

display_vars = ['const', 'price', 'log_inside_share', 'hppw', 'FuelEfficiency', 'size']
print(f"{'変数':<25} {'OLS':>10}  {'IV_BLP':>10}")
print("-" * 49)
for name in display_vars:
    print(f"{name:<25} {ols_dict[name]:>10.4f}  {iv_dict[name]:>10.4f}")
    print(f"{'':25s} ({ols_se_dict[name]:>8.4f})  ({iv_se_dict[name]:>8.4f})")

# Save
with open(output_dir / 'tab_4_2_nest_iv.txt', 'w', encoding='utf-8') as f:
    f.write(f"{'変数':<25} {'OLS':>12} {'IV_BLP':>12}\\n")
    for name in display_vars:
        f.write(f"{name:<25} {ols_dict[name]:>10.4f}   {iv_dict[name]:>10.4f}\\n")
        f.write(f"{'':25s} ({ols_se_dict[name]:>8.4f})   ({iv_se_dict[name]:>8.4f})\\n")
"""))

cells.append(nbf.v4.new_markdown_cell(
    "### 入れ子型ロジットモデルにおける弾力性 (表 4.3)"
))

cells.append(nbf.v4.new_code_cell("""\
# 入れ子型ロジットの弾力性
# IV BLP の結果から alpha2 (price), sigma2 (log_inside_share) を取得
# res_nest_blp beta order: [const, hppw, FuelEfficiency, size, price, log_inside_share]
alpha2 = res_nest_blp['beta'][4]   # price
sigma2 = res_nest_blp['beta'][5]   # log_inside_share

# OLS
# beta_ols_nest order: [const, price, log_inside_share, hppw, FuelEfficiency, size]
alpha1_ols = beta_ols_nest[1]   # price
sigma1_ols = beta_ols_nest[2]   # log_inside_share

# 全データについての自己弾力性の記述統計
data['own_elas_ols_nest'] = (alpha1_ols * data['price'] *
    (1 - sigma1_ols * data['inside_share'] -
     (1 - sigma1_ols) * data['share']) /
    (1 - sigma1_ols))

data['own_elas_ivblp_nest'] = (alpha2 * data['price'] *
    (1 - sigma2 * data['inside_share'] -
     (1 - sigma2) * data['share']) / (1 - sigma2))

tbl_nest_own_elas = data[['own_elas_ols_nest', 'own_elas_ivblp_nest']].describe().T[
    ['mean', 'std', '50%', 'min', 'max']
]
tbl_nest_own_elas.columns = ['Mean', 'Std.Dev.', 'Median', 'Min', 'Max']

print("入れ子型ロジットの自己弾力性 記述統計:")
print(tbl_nest_own_elas.round(4).to_string())
tbl_nest_own_elas.round(4).to_csv(output_dir / 'tbl_nest_own_elas.txt', sep='\\t')
"""))

cells.append(nbf.v4.new_code_cell("""\
# 日評自動車 2016年の入れ子型ロジット弾力性行列を作成
data_NIPPYO_nest = data[data['Nippyo'] == 1].copy()

dt2016_nest = data_NIPPYO_nest[data_NIPPYO_nest['year'] == 2016][
    ['price', 'Type', 'share', 'inside_share', 'NameID']
].copy().sort_values('NameID').reset_index(drop=True)

price_n = dt2016_nest['price'].values
share_n = dt2016_nest['share'].values
inside_share_n = dt2016_nest['inside_share'].values
nameID_n = dt2016_nest['NameID'].values
group_n = dt2016_nest['Type'].values
J_n = len(price_n)

# 自己弾力性
own_elas_nl = alpha2 * price_n * (1 - sigma2 * inside_share_n - (1 - sigma2) * share_n) / (1 - sigma2)

# グループ外の交差弾力性: -alpha * p_k * s_k
cross_elas_othergroup = (-1) * alpha2 * price_n * share_n
cross_elas_othergroup_mat = np.tile(cross_elas_othergroup, (J_n, 1))

# グループ内の交差弾力性
# cross_elas_samegroup[j, k] = -alpha * p_k * (sigma * inside_share_k + (1-sigma) * share_k) / (1-sigma)
price_l_mat = np.tile(price_n, (J_n, 1))
share_l_mat = np.tile(share_n, (J_n, 1))
insideshare_l_mat = np.tile(inside_share_n, (J_n, 1))
cross_elas_samegroup = (-1) * alpha2 * price_l_mat * (
    sigma2 * insideshare_l_mat + (1 - sigma2) * share_l_mat
) / (1 - sigma2)

# グループ一致のインジケータ
temp_mat1 = np.tile(group_n, (J_n, 1))
temp_mat2 = temp_mat1.T
ind_same_group = (temp_mat1 == temp_mat2).astype(float)
ind_other_group = (temp_mat1 != temp_mat2).astype(float)

# 弾力性行列を構築
elas_mat_nl = cross_elas_samegroup * ind_same_group + cross_elas_othergroup_mat * ind_other_group
np.fill_diagonal(elas_mat_nl, own_elas_nl)

# 4車種を抽出
idx_nl = [np.where(nameID_n == tid)[0][0] for tid in target_ids]
elas_mat_nl_restricted = elas_mat_nl[np.ix_(idx_nl, idx_nl)]
elas_nl_df = pd.DataFrame(elas_mat_nl_restricted, index=target_names, columns=target_names)

print("表 4.3: 入れ子型ロジットモデルにおける弾力性行列 (2016年)")
print(elas_nl_df.round(4).to_string())
elas_nl_df.round(4).to_csv(output_dir / 'tab4_3_elas_mat_nl_restricted.txt', sep='\\t')
"""))

# ============================================================
# Part 3: BLP Random Coefficients Setup (Step 1)
# ============================================================
cells.append(nbf.v4.new_markdown_cell(
    "## 5. ランダム係数ロジットモデル (BLP) の推定\n\n"
    "Berry, Levinsohn, and Pakes (1995) のランダム係数ロジットモデルを推定する。\n\n"
    "**推定のOverview**\n"
    "- Step 1: データに関する下準備\n"
    "- Step 2: 市場シェアを計算する関数\n"
    "- Step 3: 縮小写像によるBerryインバージョン\n"
    "- Step 4: GMMの目的関数\n"
    "- Step 5: 標準誤差の計算\n\n"
    "### Step 1: データの下準備"
))

cells.append(nbf.v4.new_code_cell("""\
# データのソート (マーケット順、モデル順)
data = data.sort_values(['year', 'NameID']).reset_index(drop=True)

# マーケットとモデルの情報
marketindex = data['year'].values
N = len(marketindex)
T_mkt = len(np.unique(marketindex))

print(f"N (total obs): {N}")
print(f"T (markets): {T_mkt}")

# X1: 平均効用に入る変数 [const, price, FuelEfficiency, hppw, size]
X1 = np.column_stack([
    np.ones(N),
    data['price'].values,
    data['FuelEfficiency'].values,
    data['hppw'].values,
    data['size'].values
])
X1_names = ['cons', 'price', 'FuelEfficiency', 'hppw', 'size']

# X2: ランダム係数とinteractする変数 [price, cons, size]
X2 = np.column_stack([
    data['price'].values,
    np.ones(N),
    data['size'].values
])
X2_names = ['price', 'cons', 'size']

# Z: 操作変数行列 (外生変数 + GH IV、ただしnest用を除く)
iv_GH_cols = [c for c in data.columns if c.startswith('iv_GH') and not c.endswith('nest')]
Z = np.column_stack([
    np.ones(N),
    data['FuelEfficiency'].values,
    data['hppw'].values,
    data['size'].values,
    data[iv_GH_cols].values
])

# 市場シェア
ShareVec = data['share'].values.reshape(-1, 1)

# logit share (contraction mapping 初期値)
logitshare = data['logit_share'].values

print(f"X1 shape: {X1.shape}")
print(f"X2 shape: {X2.shape}")
print(f"Z shape: {Z.shape}")
print(f"IV columns used: {iv_GH_cols}")
"""))

cells.append(nbf.v4.new_code_cell("""\
# 乱数の生成
np.random.seed(111)
Nsim = 500
K2 = X2.shape[1]  # = 3

draw_vec = np.random.randn(Nsim * K2)
draw_matrix = draw_vec.reshape(K2, Nsim)

print(f"Nsim: {Nsim}, K2: {K2}")
print(f"draw_matrix shape: {draw_matrix.shape}")

# マーケットインジケータ行列 (tempmat): N x T_mkt
unique_markets = np.sort(np.unique(marketindex))
tempmat = np.zeros((N, T_mkt), dtype=float)
for t_idx, yr in enumerate(unique_markets):
    tempmat[:, t_idx] = (marketindex == yr).astype(float)

print(f"tempmat shape: {tempmat.shape}")
"""))

# ============================================================
# BLP Core Functions (Steps 2-4)
# ============================================================
cells.append(nbf.v4.new_markdown_cell(
    "### Step 2-4: BLPの核となる関数の定義"
))

cells.append(nbf.v4.new_code_cell("""\
def calculate_mktshare(X2, theta2, delta, draw_matrix, tempmat, T_mkt, Nsim):
    \"\"\"
    シミュレーションにより市場シェアを計算する。

    Parameters
    ----------
    X2 : (N, K2) array - ランダム係数とinteractする変数
    theta2 : (K2,) array - ランダム係数の標準偏差パラメータ
    delta : (N,) or (N,1) array - 平均効用
    draw_matrix : (K2, Nsim) array - 標準正規乱数
    tempmat : (N, T_mkt) array - マーケットインジケータ
    T_mkt : int - マーケット数
    Nsim : int - シミュレーション数

    Returns
    -------
    s_jt : (N,) array - 予測市場シェア
    \"\"\"
    delta = delta.flatten()
    K2 = len(theta2)

    # 非線形の要素: mu = X2 @ diag(theta2) @ draw_matrix  -> (N, Nsim)
    mu = X2 @ np.diag(theta2) @ draw_matrix

    # delta + mu -> (N, Nsim)
    delta_mu = delta.reshape(-1, 1) + mu
    exp_delta_mu = np.exp(delta_mu)

    # outside option の分母 (=1 for each market-simulation)
    denom_outside = np.ones((T_mkt, Nsim))

    # マーケットごとの exp(delta+mu) の和
    # R code: denom_temp <- t(t(exp_delta_mu) %*% tempmat)
    #   t(exp_delta_mu) is (Nsim x N), tempmat is (N x T_mkt)
    #   result: (Nsim x T_mkt) -> transpose -> (T_mkt x Nsim)
    denom_temp = (exp_delta_mu.T @ tempmat).T + denom_outside  # (T_mkt, Nsim)

    # 各製品に対応するマーケットの分母: tempmat @ denom_temp -> (N, Nsim)
    denom = tempmat @ denom_temp

    # 選択確率 (N, Nsim)
    s_jt_i = exp_delta_mu / denom

    # シミュレーション平均 -> (N,)
    s_jt = s_jt_i.mean(axis=1)

    return s_jt


def calculate_avg_utility_by_Berry_inversion(X2, theta2, draw_matrix, tempmat,
                                              T_mkt, Nsim, ShareVec, delta_ini):
    \"\"\"
    縮小写像 (contraction mapping) によるBerryインバージョン。

    exp_delta_new = exp_delta_old * share_obs / pred_share

    Parameters
    ----------
    ShareVec : (N,1) or (N,) array - 観測市場シェア
    delta_ini : (N,) array - 初期値

    Returns
    -------
    delta : (N,) array - 収束した平均効用
    \"\"\"
    tol = 1e-11
    share_obs = ShareVec.flatten()

    delta_old = delta_ini.copy().flatten()
    exp_delta_old = np.exp(delta_old)

    for iteration in range(1000):
        pred_share = calculate_mktshare(X2, theta2, delta_old, draw_matrix,
                                         tempmat, T_mkt, Nsim)

        exp_delta_new = exp_delta_old * share_obs / pred_share

        norm_val = np.max(np.abs(exp_delta_new - exp_delta_old))

        exp_delta_old = exp_delta_new
        delta_old = np.log(exp_delta_new)

        if norm_val < tol:
            break

    return np.log(exp_delta_new)


def GMM_obj(theta2, X1, X2, Z, ShareVec, draw_matrix, tempmat, T_mkt, Nsim,
            delta_ini, option=0):
    \"\"\"
    GMM目的関数。

    Parameters
    ----------
    theta2 : (K2,) array - ランダム係数の標準偏差 [sigma_price, sigma_cons, sigma_size]
    option : int
        0: 最適化用 (スカラー返却)
        1: 推定値返却用 (dict返却)

    Returns
    -------
    scalar (option=0) or dict (option=1)
    \"\"\"
    # 縮小写像
    delta = calculate_avg_utility_by_Berry_inversion(
        X2, theta2, draw_matrix, tempmat, T_mkt, Nsim, ShareVec, delta_ini
    )

    delta_col = delta.reshape(-1, 1)

    # 重み行列 W = (Z'Z)^{-1}  (2SLS)
    W = np.linalg.inv(Z.T @ Z)

    # 線形パラメータ beta_hat = (X1'Z W Z'X1)^{-1} X1'Z W Z' delta
    beta_hat = (np.linalg.inv(X1.T @ Z @ W @ Z.T @ X1) @
                X1.T @ Z @ W @ Z.T @ delta_col)

    # 構造誤差項
    Xi = delta_col - X1 @ beta_hat

    # GMM目的関数値
    obj_val = float(Xi.T @ Z @ W @ Z.T @ Xi)

    if option == 0:
        return obj_val
    elif option == 1:
        return {
            'obj_val': obj_val,
            'beta_hat': beta_hat,
            'delta': delta
        }


def calculate_standard_error(theta2, X1, X2, Z, ShareVec, draw_matrix, tempmat,
                              T_mkt, Nsim, N, delta_ini):
    \"\"\"
    数値微分による漸近標準誤差の計算。

    Returns
    -------
    Ase : (K1 + K2,) array - 漸近標準誤差
    \"\"\"
    # 縮小写像
    delta = calculate_avg_utility_by_Berry_inversion(
        X2, theta2, draw_matrix, tempmat, T_mkt, Nsim, ShareVec, delta_ini
    )
    delta_col = delta.reshape(-1, 1)

    W = np.linalg.inv(Z.T @ Z)

    beta_hat = (np.linalg.inv(X1.T @ Z @ W @ Z.T @ X1) @
                X1.T @ Z @ W @ Z.T @ delta_col)

    Xi = delta_col - X1 @ beta_hat

    # Omega_hat
    n_z = Z.shape[1]
    Omega_hat = np.zeros((n_z, n_z))
    for ii in range(N):
        z_i = Z[ii, :].reshape(-1, 1)
        Omega_hat += (z_i @ z_i.T) * (Xi[ii, 0] ** 2) / N

    # delta の theta2 に関する数値微分
    K2_dim = len(theta2)
    Ddelta = np.zeros((N, K2_dim))
    h = 1e-6

    for k in range(K2_dim):
        theta2_plus = theta2.copy()
        theta2_plus[k] += h

        delta_plus = calculate_avg_utility_by_Berry_inversion(
            X2, theta2_plus, draw_matrix, tempmat, T_mkt, Nsim, ShareVec, delta_ini
        )

        Ddelta[:, k] = (delta_plus - delta) / h

    # G = N^{-1} * Z' @ [-X1, Ddelta]
    G = (1.0 / N) * Z.T @ np.hstack([-X1, Ddelta])

    # 漸近分散共分散行列
    GWG_inv = np.linalg.inv(G.T @ W @ G)
    AsyVarMat = GWG_inv @ G.T @ W @ Omega_hat @ W @ G @ GWG_inv

    # 漸近標準誤差
    Ase = np.sqrt(np.diag(AsyVarMat) / N)

    return Ase

print("BLPの関数を定義しました。")
print("  - calculate_mktshare()")
print("  - calculate_avg_utility_by_Berry_inversion()")
print("  - GMM_obj()")
print("  - calculate_standard_error()")
"""))

# ============================================================
# BLP Elasticity Function
# ============================================================
cells.append(nbf.v4.new_code_cell("""\
def calculate_elasticity(X1, X2, theta2, delta, draw_matrix, tempmat,
                          marketindex, T_mkt, Nsim, beta_hat, beta_names):
    \"\"\"
    シミュレーションベースの価格弾力性行列を年ごとに計算する。

    Parameters
    ----------
    beta_names : dict - パラメータ名 -> 値 のマッピング

    Returns
    -------
    elaslist : dict of {year: elasticity_matrix}
    \"\"\"
    K2 = len(theta2)

    # mu の計算 (calculate_mktshare と同じ)
    mu = X2 @ np.diag(theta2) @ draw_matrix
    delta_flat = delta.flatten()
    delta_mu = delta_flat.reshape(-1, 1) + mu
    exp_delta_mu = np.exp(delta_mu)

    denom_outside = np.ones((T_mkt, Nsim))
    denom_temp = (exp_delta_mu.T @ tempmat).T + denom_outside
    denom = tempmat @ denom_temp

    # 選択確率 (N, Nsim)
    s_jt_i = exp_delta_mu / denom

    # 価格列
    price_all = X1[:, 1]  # X1の2列目がprice

    # alpha_i = beta_price_mean + sigma_price * draw_i
    alpha_mean = beta_names['price_mean']
    sigma_price = beta_names['sigma_price']
    draw_for_price = draw_matrix[0, :]  # 1行目がprice用乱数
    alpha_i = alpha_mean + sigma_price * draw_for_price  # (Nsim,)

    unique_years = np.sort(np.unique(marketindex))
    elaslist = {}

    for t_idx, yr in enumerate(unique_years):
        mask = (marketindex == yr)
        J_t = mask.sum()

        # 製品ごと・個人ごとの選択確率 (J_t, Nsim)
        ag_model_s_i = s_jt_i[mask, :]

        # 製品ごとの選択確率 (J_t,)
        ag_model_s = ag_model_s_i.mean(axis=1)

        # 価格 (J_t,)
        price_t = price_all[mask]

        # ベクトル化された弾力性計算
        weighted_s = ag_model_s_i * alpha_i[None, :]  # (J_t, Nsim)
        cross_term = weighted_s @ ag_model_s_i.T / Nsim  # (J_t, J_t)

        elasmat = -price_t[:, None] / ag_model_s[None, :] * cross_term

        # 対角要素（自己価格弾力性）を上書き
        own_term = np.mean(alpha_i[None, :] * ag_model_s_i * (1 - ag_model_s_i), axis=1)
        np.fill_diagonal(elasmat, price_t / ag_model_s * own_term)

        elaslist[yr] = elasmat

    return elaslist

print("弾力性計算関数を定義しました。")
"""))

# ============================================================
# GMM Optimization
# ============================================================
cells.append(nbf.v4.new_markdown_cell("### GMM最適化"))

cells.append(nbf.v4.new_code_cell("""\
# GMM最適化
theta2_init = np.array([0.2, 10.0, 0.1])

print("GMM最適化を開始します...")
start_time = time.time()

result = optimize.minimize(
    GMM_obj,
    theta2_init,
    args=(X1, X2, Z, ShareVec, draw_matrix, tempmat, T_mkt, Nsim, logitshare, 0),
    method='L-BFGS-B',
    bounds=[(0, None), (0, None), (0, None)]
)

elapsed = time.time() - start_time
print(f"最適化完了: {elapsed:.1f} 秒")
print(f"theta2 = {result.x}")
print(f"GMM obj = {result.fun:.6f}")
print(f"収束: {result.success}")
"""))

# ============================================================
# Results (Tab 4.4)
# ============================================================
cells.append(nbf.v4.new_markdown_cell(
    "### Step 5: 推定結果 (表 4.4)"
))

cells.append(nbf.v4.new_code_cell("""\
# 推定結果の取得
result2 = GMM_obj(result.x, X1, X2, Z, ShareVec, draw_matrix, tempmat,
                   T_mkt, Nsim, logitshare, option=1)

delta_est = result2['delta']
beta_hat_linear = result2['beta_hat'].flatten()

print("線形パラメータ (beta_hat):")
for name, val in zip(X1_names, beta_hat_linear):
    print(f"  {name}: {val:.6f}")

# 標準誤差の計算
print("\\n標準誤差を計算中（数値微分のため時間がかかります）...")
start_time = time.time()
se = calculate_standard_error(result.x, X1, X2, Z, ShareVec, draw_matrix,
                               tempmat, T_mkt, Nsim, N, logitshare)
elapsed = time.time() - start_time
print(f"標準誤差の計算完了: {elapsed:.1f} 秒")
"""))

cells.append(nbf.v4.new_code_cell("""\
# パラメータをまとめる
# beta_hat_linear: [cons, price, FuelEfficiency, hppw, size]
# theta2 (result.x): [sigma_price, sigma_cons, sigma_size]
all_params = np.concatenate([beta_hat_linear, result.x])

param_names_jp = ['定数項：平均', '価格：平均', '燃費', '馬力', 'サイズ：平均',
                   '価格：標準偏差', '定数項：標準偏差', 'サイズ：標準偏差']

# 表示順に並べ替え（R版に合わせる）
display_order = ['定数項：平均', '定数項：標準偏差',
                 '価格：平均', '価格：標準偏差',
                 'サイズ：平均', 'サイズ：標準偏差',
                 '燃費', '馬力']

param_dict = dict(zip(param_names_jp, all_params))
se_dict = dict(zip(param_names_jp, se))

result_table = pd.DataFrame({
    '推定値': [param_dict[n] for n in display_order],
    '標準誤差': [se_dict[n] for n in display_order]
}, index=display_order)

print("表 4.4: ランダム係数ロジットモデルの推定結果")
print(result_table.round(4).to_string())
result_table.round(4).to_csv(output_dir / 'tab4_4_rand_coef_logit_result.txt', sep='\\t')
"""))

# ============================================================
# BLP Elasticity (Tab 4.5)
# ============================================================
cells.append(nbf.v4.new_markdown_cell(
    "### 価格弾力性行列 (表 4.5)"
))

cells.append(nbf.v4.new_code_cell("""\
# beta_names を構築 (弾力性計算用)
beta_names = {
    'price_mean': param_dict['価格：平均'],
    'sigma_price': param_dict['価格：標準偏差'],
}

elasticity = calculate_elasticity(
    X1, X2, result.x, delta_est, draw_matrix, tempmat,
    marketindex, T_mkt, Nsim, beta_hat_linear, beta_names
)

# 2016年の弾力性行列を取り出す
NameID2016 = data.loc[data['year'] == 2016, 'NameID'].values
elasmat_2016 = elasticity[2016]

# 4車種を抽出
idx_blp = [np.where(NameID2016 == tid)[0][0] for tid in target_ids]
elas_mat_blp_restricted = elasmat_2016[np.ix_(idx_blp, idx_blp)]
elas_blp_df = pd.DataFrame(elas_mat_blp_restricted, index=target_names, columns=target_names)

print("表 4.5: BLPモデルにおける弾力性行列 (2016年)")
print(elas_blp_df.round(4).to_string())
elas_blp_df.round(4).to_csv(output_dir / 'tab4_5_elas_mat_blp_restricted.txt', sep='\\t')
"""))

# ============================================================
# Part 4: Pricing Application
# ============================================================
cells.append(nbf.v4.new_markdown_cell(
    "## 6. 応用：プライシング\n\n"
    "ベータード（アルファード）の価格を変化させた場合の収入・利潤を計算し、\n"
    "最適価格を求める。"
))

cells.append(nbf.v4.new_code_cell("""\
# 日評自動車のNameID一覧
NIPPYOautoIDvec = data.loc[data['Nippyo'] == 1, 'NameID'].unique()

# ベータードの限界費用 (ラーナーの公式より簡便的に定義)
mc_betado = 3.198 * (1 - 1 / abs(-2.16720791))

NameID_target = betard

def calculate_revenue(price_cand, data_df, X1_orig, X2_orig, beta_hat_est, delta_est,
                       draw_matrix, tempmat, T_mkt, Nsim, theta2,
                       betard_id, nippyo_ids, option='ownpi'):
    \"\"\"
    ベータードの価格を変えた場合の収入・利潤を計算する。

    Parameters
    ----------
    price_cand : float - ベータードの候補価格
    option : str - 'own': ベータード収入, 'total': 日評全体収入,
                   'ownpi': ベータード利潤, 'totalpi': ベータード利潤+他収入
    \"\"\"
    mc_betado = 3.198 * (1 - 1 / abs(-2.16720791))

    # 価格をコピーして差し替え
    tempprice = data_df['price'].values.copy()
    mask_betard_2016 = ((data_df['NameID'].values == betard_id) &
                         (data_df['year'].values == 2016))
    tempprice[mask_betard_2016] = price_cand

    # X1, X2のpriceを差し替え
    X1_temp = X1_orig.copy()
    X2_temp = X2_orig.copy()
    X1_temp[:, 1] = tempprice   # X1のprice列
    X2_temp[:, 0] = tempprice   # X2のprice列

    # 新たな価格で平均効用を計算し直す
    org_xi = delta_est.reshape(-1, 1) - X1_orig @ beta_hat_est.reshape(-1, 1)
    new_delta = (X1_temp @ beta_hat_est.reshape(-1, 1) + org_xi).flatten()

    # 市場シェアを計算
    mktshare = calculate_mktshare(X2_temp, theta2, new_delta, draw_matrix,
                                   tempmat, T_mkt, Nsim)

    quant = mktshare * data_df['HH'].values
    revenue = tempprice * quant

    # ベータードのみの収入
    rev_betard = revenue[mask_betard_2016].sum()

    # 日評全体の収入
    mask_nippyo_2016 = (np.isin(data_df['NameID'].values, nippyo_ids) &
                         (data_df['year'].values == 2016))
    rev_total = revenue[mask_nippyo_2016].sum()

    # ベータード利潤
    quant_betard = quant[mask_betard_2016].sum()
    pi_betard = rev_betard - mc_betado * quant_betard

    # ベータード利潤 + 他収入
    pi_total = rev_total - mc_betado * quant_betard

    if option == 'own':
        return rev_betard
    elif option == 'total':
        return rev_total
    elif option == 'ownpi':
        return pi_betard
    elif option == 'totalpi':
        return pi_total

print(f"ベータード NameID: {betard}")
print(f"限界費用 mc: {mc_betado:.4f}")
print(f"日評自動車 車種数: {len(NIPPYOautoIDvec)}")
"""))

cells.append(nbf.v4.new_code_cell("""\
# 価格グリッドで利潤曲線を計算
pricevec = np.arange(1.73, 4.01, 0.01)
pivec = np.zeros(len(pricevec))
pivec2 = np.zeros(len(pricevec))

print("利潤曲線を計算中...")

for i, p in enumerate(pricevec):
    pivec[i] = calculate_revenue(
        p, data, X1, X2, beta_hat_linear, delta_est,
        draw_matrix, tempmat, T_mkt, Nsim, result.x,
        betard, NIPPYOautoIDvec, option='ownpi'
    )
    pivec2[i] = calculate_revenue(
        p, data, X1, X2, beta_hat_linear, delta_est,
        draw_matrix, tempmat, T_mkt, Nsim, result.x,
        betard, NIPPYOautoIDvec, option='totalpi'
    )
    if (i + 1) % 50 == 0:
        print(f"  {i + 1}/{len(pricevec)} 完了")

print("計算完了")
"""))

cells.append(nbf.v4.new_code_cell("""\
# ベータードのみの利潤曲線
fig, ax = plt.subplots(figsize=(9, 6))
ax.plot(pricevec * 100, pivec * 100 / 10000, 'o', markersize=2)
ax.set_xlabel('価格(万円)')
ax.set_ylabel('収入(億円)')
ax.set_title('ベータードのみの利潤')
ax.set_xticks(np.arange(150, 450, 50))
ax.grid(axis='y', linestyle='dotted', alpha=0.7)
plt.tight_layout()
plt.savefig(output_dir / 'fig4_2_revenue_betard.png', dpi=150)
plt.show()
"""))

cells.append(nbf.v4.new_code_cell("""\
# ベータード利潤 + 日評自動車他車種の収入
fig, ax = plt.subplots(figsize=(9, 6))
ax.plot(pricevec * 100, pivec2 * 100 / 10000, 'o', markersize=2)
ax.set_xlabel('価格(万円)')
ax.set_ylabel('収入(億円)')
ax.set_title('ベータード利潤 + 日評自動車他車種の収入')
ax.set_xticks(np.arange(150, 450, 50))
ax.grid(axis='y', linestyle='dotted', alpha=0.7)
plt.tight_layout()
plt.savefig(output_dir / 'fig4_2_revenue_all.png', dpi=150)
plt.show()
"""))

cells.append(nbf.v4.new_code_cell("""\
# まとめたテーブルを保存
dt_pricing = pd.DataFrame({
    'price': pricevec * 100,
    'pi1': pivec * 100 / 10000,
    'pi2': pivec2 * 100 / 10000
})
dt_pricing.to_csv(output_dir / 'dt_prof_Betard.csv', index=False)
print("dt_prof_Betard.csv を保存しました。")
"""))

# ============================================================
# Optimal Pricing
# ============================================================
cells.append(nbf.v4.new_markdown_cell("### 収入最大化価格の計算"))

cells.append(nbf.v4.new_code_cell("""\
# ベータードのみの利潤を最大化
def neg_revenue_own(price_cand):
    return -calculate_revenue(
        price_cand, data, X1, X2, beta_hat_linear, delta_est,
        draw_matrix, tempmat, T_mkt, Nsim, result.x,
        betard, NIPPYOautoIDvec, option='ownpi'
    )

optim_own = optimize.minimize_scalar(neg_revenue_own, bounds=(0.3, 5.0), method='bounded')

print("ベータードの利潤のみ最適化:")
print(f"  最適価格: {optim_own.x * 100:.2f} 万円")
print(f"  最大利潤: {-optim_own.fun:.6f}")

# ベータードのみの利潤を最大にする価格での、全体収入
total_rev_at_own_opt = calculate_revenue(
    optim_own.x, data, X1, X2, beta_hat_linear, delta_est,
    draw_matrix, tempmat, T_mkt, Nsim, result.x,
    betard, NIPPYOautoIDvec, option='totalpi'
)
print(f"  その価格でのベータード利潤+他収入: {total_rev_at_own_opt:.6f}")
"""))

cells.append(nbf.v4.new_code_cell("""\
# ベータード利潤 + 日評全体の収入を最大化
def neg_revenue_total(price_cand):
    return -calculate_revenue(
        price_cand, data, X1, X2, beta_hat_linear, delta_est,
        draw_matrix, tempmat, T_mkt, Nsim, result.x,
        betard, NIPPYOautoIDvec, option='totalpi'
    )

optim_total = optimize.minimize_scalar(neg_revenue_total, bounds=(0.3, 5.0), method='bounded')

print("ベータード利潤+他の収入最大化:")
print(f"  最適価格: {optim_total.x * 100:.2f} 万円")
print(f"  最大値: {-optim_total.fun:.6f}")
"""))

cells.append(nbf.v4.new_code_cell("""\
# 最適価格のまとめを保存
with open(output_dir / 'Ch04_opt_price.txt', 'w', encoding='utf-8') as f:
    f.write("ベータードの利潤のみ最適化\\n")
    f.write(f"  最適価格: {optim_own.x:.6f} (100万円単位)\\n")
    f.write(f"  最大利潤: {-optim_own.fun:.6f}\\n")
    f.write(f"  その価格でのベータード利潤+他収入: {total_rev_at_own_opt:.6f}\\n")
    f.write("\\nベータード利潤+他の収入最大化\\n")
    f.write(f"  最適価格: {optim_total.x:.6f} (100万円単位)\\n")
    f.write(f"  最大値: {-optim_total.fun:.6f}\\n")

print("Ch04_opt_price.txt を保存しました。")
"""))

# ============================================================
# Summary
# ============================================================
cells.append(nbf.v4.new_markdown_cell("## まとめ"))

cells.append(nbf.v4.new_code_cell("""\
print("=" * 60)
print("第4章の分析完了")
print("=" * 60)
print("\\n出力ファイル:")
for f in sorted(output_dir.glob('*')):
    if f.name.startswith(('tab4_', 'tab_4_', 'fig4_', 'Ch04_', 'dt_prof_', 'tbl_nest')):
        print(f"  {f.name}")
print("\\n中間ファイル:")
print("  data_for_estimation.csv")
"""))

nb.cells = cells
nbf.write(nb, 'main_ch04.ipynb')
print("Generated: main_ch04.ipynb")
