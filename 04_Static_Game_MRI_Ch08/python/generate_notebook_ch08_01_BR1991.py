"""Generate main_ch08_01_BR1991.ipynb for Chapter 8: Bresnahan-Reiss (1991) Entry Model."""
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
# Cell 1: Title
# ============================================================
cells.append(nbf.v4.new_markdown_cell(
    "# 第8章 静学ゲームの推定：Bresnahan and Reiss (1991) 参入モデル\n"
    "\n"
    "MRIスキャナーの導入に関する参入モデルを推定する。\n"
    "Bresnahan and Reiss (1991) の手法を用い、参入閾値と競争効果を分析する。"
))

# ============================================================
# Cell 2: Setup
# ============================================================
cells.append(nbf.v4.new_code_cell("""\
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.stats import norm
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib
import warnings
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
output_dir = base_dir / 'output'
output_dir.mkdir(exist_ok=True)

print("Setup complete.")
"""))

# ============================================================
# Cell 3: Data Loading & Preparation
# ============================================================
cells.append(nbf.v4.new_markdown_cell("## データの読み込みと準備"))

cells.append(nbf.v4.new_code_cell("""\
# 病院レベルのデータを読み込み
data = pd.read_csv(base_dir / 'data' / 'MRIData_BR1991.csv', encoding='utf-8')

print(f"病院レベルのデータ: {len(data)} 行")
print(f"カラム: {list(data.columns)}")
data.head()
"""))

cells.append(nbf.v4.new_code_cell("""\
# 病院レベルのデータを市区町村レベルの集計データに変換
listCode = data['CityCode'].unique()

records = []
for code in listCode:
    subdata = data[data['CityCode'] == code]
    records.append({
        'Code': code,
        'NumHospital': len(subdata),
        'NumMRI': subdata['MRIOwnDum'].sum(),
        'Pop': subdata['Population'].unique()[0],
        'Menseki': subdata['Menseki'].unique()[0],
        'PopDen': subdata['PopDensity'].unique()[0],
        'Income': subdata['TaxableIncome'].unique()[0],
    })

dataset = pd.DataFrame(records)
print(f"市区町村数: {len(dataset)}")
dataset.head()
"""))

# ============================================================
# Cell 4: Cross-tabulation Table (Tab 8.3)
# ============================================================
cells.append(nbf.v4.new_markdown_cell("## 表8.3: 病院数とMRI保有病院数のクロス集計"))

cells.append(nbf.v4.new_code_cell("""\
# 10以上は10としてまとめる
dt_table = dataset[['NumHospital', 'NumMRI']].copy()
dt_table['NumHospital'] = dt_table['NumHospital'].clip(upper=10)
dt_table['NumMRI'] = dt_table['NumMRI'].clip(upper=10)

# クロス集計
tbl = pd.crosstab(dt_table['NumHospital'], dt_table['NumMRI'])
print("病院数 x MRI保有病院数のクロス集計:")
print(tbl)

# CSV保存
tbl.to_csv(output_dir / 'Tab8_3_hospital_mri_table.csv')
print("\\n保存: output/Tab8_3_hospital_mri_table.csv")
"""))

# ============================================================
# Cell 5: Define Objective Function
# ============================================================
cells.append(nbf.v4.new_markdown_cell(
    "## 目的関数の定義\n"
    "\n"
    "Bresnahan and Reiss (1991) モデルの対数尤度関数を定義する。"
))

cells.append(nbf.v4.new_code_cell("""\
def obj(params, dataset_NumMRI, dataset_Pop, M, N_max):
    \"\"\"
    Bresnahan-Reiss (1991) モデルの目的関数（正規化対数尤度）。

    Parameters
    ----------
    params : array, shape (N_max + 1,)
        パラメータ。params[0]=alpha_1, params[1:N_max]=-alpha_2,...,-alpha_{N_max}, params[N_max]=gamma
    dataset_NumMRI : array, shape (M,)
        各市区町村のMRI導入病院数
    dataset_Pop : array, shape (M,)
        各市区町村の人口（百万人単位）
    M : int
        市区町村数
    N_max : int
        最大参入企業数

    Returns
    -------
    float
        正規化対数尤度 (sum of log-likelihood / M)
    \"\"\"
    # パラメータの定義
    alpha1 = params[0]
    alpha2 = -params[1:N_max]  # alpha_2, ..., alpha_{N_max} にはマイナスを乗じる
    alpha = np.concatenate([[alpha1], alpha2])
    gamma = params[N_max]

    # 人口の行列 (M x N_max)
    pop = np.tile(dataset_Pop.reshape(-1, 1), (1, N_max))

    # 下三角行列 V を構築
    # V[i:, i] = alpha[i] (i = 0, ..., len(alpha)-1)
    V = np.zeros((N_max, N_max))
    for i in range(len(alpha)):
        V[i:, i] = alpha[i]

    # VV: 可変利潤部分 (M x N_max)
    VV = (V @ np.ones((N_max, M))).T

    # 固定費用の行列 (M x N_max)
    F = np.full((M, N_max), gamma)

    # 利潤行列
    pi_mat = pop * VV - F

    # 標準正規分布のCDFを適用
    phi = norm.cdf(pi_mat)

    # 各企業数になる確率を計算
    # n=0: 1 - phi[:,0]
    # n=k (0 < k < N_max): phi[:,k-1] - phi[:,k]
    # n=N_max: phi[:,N_max-1]
    mat = np.column_stack([
        1 - phi[:, 0],
        phi[:, :N_max - 1] - phi[:, 1:N_max],
        phi[:, N_max - 1]
    ])

    # 観測された企業数に対応する確率の対数を抽出
    ml = np.zeros(M)
    for i in range(N_max + 1):
        mask = dataset_NumMRI == i
        ml[mask] = np.log(np.maximum(mat[mask, i], 1e-300))

    # 対数尤度が定義できない場合の対処
    ml[np.isinf(ml)] = -10000

    # 正規化対数尤度
    val = np.sum(ml) / M
    return val


print("目的関数の定義完了")
"""))

# ============================================================
# Cell 6: Numerical Hessian Function
# ============================================================
cells.append(nbf.v4.new_markdown_cell("## 数値ヘシアンの関数定義"))

cells.append(nbf.v4.new_code_cell("""\
def numerical_hessian(func, x, eps=1e-5):
    \"\"\"
    有限差分法による数値ヘシアンの計算。

    Parameters
    ----------
    func : callable
        スカラー値を返す関数
    x : array
        パラメータの値
    eps : float
        差分の幅

    Returns
    -------
    H : array, shape (n, n)
        ヘシアン行列
    \"\"\"
    n = len(x)
    H = np.zeros((n, n))
    f0 = func(x)
    for i in range(n):
        for j in range(n):
            e_i = np.zeros(n)
            e_j = np.zeros(n)
            e_i[i] = eps
            e_j[j] = eps
            H[i, j] = (func(x + e_i + e_j)
                        - func(x + e_i)
                        - func(x + e_j)
                        + f0) / eps**2
    return H


print("数値ヘシアン関数の定義完了")
"""))

# ============================================================
# Cell 7: Estimation Loop
# ============================================================
cells.append(nbf.v4.new_markdown_cell(
    "## Bresnahan and Reiss (1991) モデルの推定\n"
    "\n"
    "N_max = 6, 7, 8, 9, 10 について推定を行い、参入閾値を計算する。"
))

cells.append(nbf.v4.new_code_cell("""\
# 変数に欠損がある市区町村を落とす
dataset_clean = dataset.dropna().copy()

# 人口を百万で除する
dataset_clean['Pop'] = dataset_clean['Pop'] / 1_000_000

# 推定で用いるサンプルサイズ（市区町村数）
M = len(dataset_clean)
print(f"推定サンプルサイズ（市区町村数）: {M}")

# オリジナルのデータセットを保持
dataset_org = dataset_clean.copy()

# N_maxを6から10まで試す
N_max_list = list(range(6, 11))

# 結果の保管場所
result_estimates = {}
result_thresholds = {}

for N_max in N_max_list:
    print(f"\\n{'='*60}")
    print(f"N_max: {N_max}")
    print(f"{'='*60}")

    # データセットを初期化
    ds = dataset_org.copy()

    # MRI導入病院数がN_maxより大きい場合、N_maxに置換
    ds.loc[ds['NumMRI'] > N_max, 'NumMRI'] = N_max

    # 推定に使う変数を抽出
    dataset_NumMRI = ds['NumMRI'].values.astype(float)
    dataset_Pop = ds['Pop'].values.astype(float)

    # パラメータの初期値を全て1に設定
    initial = np.ones(N_max + 1)

    # 最適化（最大化するので目的関数を符号反転）
    result = minimize(
        lambda p: -obj(p, dataset_NumMRI, dataset_Pop, M, N_max),
        initial,
        method='L-BFGS-B',
        bounds=[(0, None)] * (N_max + 1)
    )

    estimates = result.x
    print(f"最適化成功: {result.success}")
    print(f"目的関数値（正規化対数尤度）: {-result.fun:.6f}")

    # ヘシアンの計算
    def neg_obj_for_hessian(p):
        return -obj(p, dataset_NumMRI, dataset_Pop, M, N_max)

    H = numerical_hessian(neg_obj_for_hessian, estimates)

    # 標準誤差の計算
    try:
        se = np.sqrt(np.diag(np.linalg.inv(H) / M))
    except np.linalg.LinAlgError:
        print("Warning: ヘシアンの逆行列が計算できません")
        se = np.full(N_max + 1, np.nan)

    # 推定値と標準誤差を表示
    print("\\n推定結果:")
    param_names = [f"alpha_1"] + [f"-alpha_{i}" for i in range(2, N_max + 1)] + ["gamma"]
    for name, est, s in zip(param_names, estimates, se):
        print(f"  {name:>12s}: {est:12.6f}  (se: {s:.8f})")

    # 結果を保管
    result_estimates[f"N_max={N_max}"] = np.column_stack([estimates, se])

    # === Entry Threshold の計算 ===
    alpha = estimates[:N_max]
    # alpha[0] = alpha_1, estimates[1:N_max] に -1 を乗じたものが alpha_2,...
    # ただし obj 内では alpha = [alpha1, -params[1:N_max]] としているので
    # 実際の alpha ベクトルは alpha_1, -estimates[1], -estimates[2], ...
    alpha_vals = np.zeros(N_max)
    alpha_vals[0] = estimates[0]
    alpha_vals[1:] = -estimates[1:N_max]

    gamma_val = estimates[N_max]

    EntryThreshold = np.zeros((N_max, 3))  # S_N, s_N, ratio

    # S_1, s_1
    deno = alpha_vals[0]
    S_N = gamma_val / deno * 1e6
    EntryThreshold[0, 0] = int(S_N)
    EntryThreshold[0, 1] = int(S_N)

    # S_n, s_n (n >= 2)
    for i in range(1, N_max):
        deno = deno + alpha_vals[i]  # alpha_vals[i] は負なので足す = 引く
        S_N = gamma_val / deno * 1e6
        EntryThreshold[i, 0] = int(S_N)
        EntryThreshold[i, 1] = int(S_N / (i + 1))

    # 比率 s_{N+1}/s_N
    for j in range(N_max):
        if j < N_max - 1:
            EntryThreshold[j, 2] = EntryThreshold[j + 1, 1] / EntryThreshold[j, 1]
        else:
            EntryThreshold[j, 2] = np.nan

    print("\\nEntry Threshold:")
    print(f"  {'N':>3s}  {'S_N':>10s}  {'s_N=S_N/N':>12s}  {'s_{N+1}/s_N':>14s}")
    for i in range(N_max):
        ratio_str = f"{EntryThreshold[i, 2]:.6f}" if not np.isnan(EntryThreshold[i, 2]) else "NA"
        print(f"  {i+1:>3d}  {EntryThreshold[i, 0]:>10.0f}  {EntryThreshold[i, 1]:>12.0f}  {ratio_str:>14s}")

    result_thresholds[f"N_max={N_max}"] = EntryThreshold
"""))

# ============================================================
# Cell 8: Save Results
# ============================================================
cells.append(nbf.v4.new_markdown_cell("## 結果の保存"))

cells.append(nbf.v4.new_code_cell("""\
# 推定値の保存
with open(output_dir / 'Tab8_4_BR1991_Estimates.txt', 'w') as f:
    for name, est in result_estimates.items():
        f.write(f"\\n{name}\\n")
        f.write(f"{'estimates':>14s} {'se':>14s}\\n")
        for row in est:
            f.write(f"{row[0]:14.6f} {row[1]:14.8f}\\n")

print("保存: output/Tab8_4_BR1991_Estimates.txt")

# Entry Threshold の保存
with open(output_dir / 'Tab8_5_Entry_Thresholds.txt', 'w') as f:
    for name, thr in result_thresholds.items():
        f.write(f"\\n{name}\\n")
        f.write(f"{'S_N':>10s} {'s_N=S_N/N':>12s} {'s_{N+1}/s_N':>14s}\\n")
        for row in thr:
            ratio_str = f"{row[2]:14.6f}" if not np.isnan(row[2]) else f"{'NA':>14s}"
            f.write(f"{row[0]:10.0f} {row[1]:12.0f} {ratio_str}\\n")

print("保存: output/Tab8_5_Entry_Thresholds.txt")
"""))

# ============================================================
# Cell 9: Results Summary
# ============================================================
cells.append(nbf.v4.new_markdown_cell("## 結果のまとめ"))

cells.append(nbf.v4.new_code_cell("""\
# 全てのN_maxの推定結果をまとめて表示
print("=" * 80)
print("Bresnahan and Reiss (1991) 推定結果のまとめ")
print("=" * 80)

for name in result_estimates:
    est = result_estimates[name]
    thr = result_thresholds[name]
    N_max = int(name.split("=")[1])

    print(f"\\n--- {name} ---")
    print(f"  パラメータ推定値:")
    print(f"    {'param':>12s}  {'estimate':>12s}  {'se':>12s}")
    param_names = ["alpha_1"] + [f"-alpha_{i}" for i in range(2, N_max + 1)] + ["gamma"]
    for i, pname in enumerate(param_names):
        print(f"    {pname:>12s}  {est[i, 0]:12.6f}  {est[i, 1]:12.8f}")

    print(f"  Entry Threshold:")
    print(f"    {'N':>3s}  {'S_N':>10s}  {'s_N=S_N/N':>12s}  {'s_{N+1}/s_N':>14s}")
    for i in range(N_max):
        ratio_str = f"{thr[i, 2]:.6f}" if not np.isnan(thr[i, 2]) else "NA"
        print(f"    {i+1:>3d}  {thr[i, 0]:>10.0f}  {thr[i, 1]:>12.0f}  {ratio_str:>14s}")

print("\\n" + "=" * 80)
print("推定完了")
"""))

nb.cells = cells
nbf.write(nb, 'main_ch08_01_BR1991.ipynb')
print("Generated: main_ch08_01_BR1991.ipynb")
