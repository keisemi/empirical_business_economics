"""Generate main_ch06.ipynb for Chapter 6: Dynamic Discrete Choice (NFXP)."""
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
# Cell 1: Setup
# ============================================================
cells.append(nbf.v4.new_markdown_cell(
    "# 第6章 動的離散選択モデル：入れ子不動点アルゴリズム\n"
    "\n"
    "中古車の購入に関する動的離散選択モデルを推定する。\n"
    "Rust (1987) 型のモデルを用い、入れ子不動点 (NFXP) アルゴリズムで構造パラメータを推定する。"
))

cells.append(nbf.v4.new_code_cell("""\
import numpy as np
import pandas as pd
from scipy import optimize, stats
from pathlib import Path
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
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
intermediate_dir = base_dir / 'intermediate'
output_dir = base_dir / 'output'
intermediate_dir.mkdir(exist_ok=True)
output_dir.mkdir(exist_ok=True)

print("Setup complete.")
"""))

# ============================================================
# Cell 2: Parameters & State Space
# ============================================================
cells.append(nbf.v4.new_markdown_cell("## パラメータの設定とState空間の作成"))

cells.append(nbf.v4.new_code_cell("""\
# パラメータの設定
theta_true = np.array([15.0, 6.0])  # theta_c, theta_p
beta = 0.95                          # 時間割引率
Euler_const = -np.euler_gamma        # オイラー定数 (= digamma(1) in R)
num_choice = 2                       # 選択肢の数（購入する/しない）

# 状態空間の作成
price_states = np.round(np.arange(2.0, 2.6, 0.1), 1)     # 価格 (2.0 ~ 2.5): 6個
mileage_states = np.round(np.arange(0, 0.105, 0.005), 3)  # 走行距離 (0 ~ 0.1): 21個
# np.arangeの浮動小数点誤差で要素数が変わるのを防止
price_states = price_states[:6]
mileage_states = mileage_states[:21]

num_price_states = len(price_states)      # 6
num_mileage_states = len(mileage_states)  # 21
num_states = num_price_states * num_mileage_states  # 126

# 状態変数のデータフレーム
# 順番: (p,m) = (2.0,0), (2.1,0), ..., (2.5,0), (2.0,0.005), (2.1,0.005), ...
state_df = pd.DataFrame({
    'state_id': np.arange(1, num_states + 1),
    'price_id': np.tile(np.arange(1, num_price_states + 1), num_mileage_states),
    'mileage_id': np.repeat(np.arange(1, num_mileage_states + 1), num_price_states),
    'price': np.tile(price_states, num_mileage_states),
    'mileage': np.repeat(mileage_states, num_price_states)
})

print(f"状態空間: {num_price_states} prices x {num_mileage_states} mileage = {num_states} states")
state_df.tail(3)
"""))

# ============================================================
# Cell 3: Helper Functions
# ============================================================
cells.append(nbf.v4.new_markdown_cell("## 関数定義"))

cells.append(nbf.v4.new_code_cell("""\
def gen_mileage_trans(kappa):
    \"\"\"走行距離の遷移行列を生成する（購買なし/ありの2つ）\"\"\"
    kappa_1, kappa_2 = kappa

    # 購買しなかった場合の遷移行列
    mat_not_buy = np.zeros((num_mileage_states, num_mileage_states))
    for i in range(num_mileage_states):
        for j in range(num_mileage_states):
            if i == j:
                mat_not_buy[i, j] = 1 - kappa_1 - kappa_2
            elif i == j - 1:
                mat_not_buy[i, j] = kappa_1
            elif i == j - 2:
                mat_not_buy[i, j] = kappa_2
    # 境界条件
    mat_not_buy[num_mileage_states - 2, num_mileage_states - 1] = kappa_1 + kappa_2
    mat_not_buy[num_mileage_states - 1, num_mileage_states - 1] = 1

    # 購買した場合の遷移行列（購入後 m=0 にリセット）
    mat_buy = np.tile(mat_not_buy[0, :], (num_mileage_states, 1))

    # (num_mileage_states, num_mileage_states, 2) の3次元配列
    return np.stack([mat_not_buy, mat_buy], axis=2)


def gen_price_trans(lam):
    \"\"\"価格の遷移行列を生成する (Rコードに忠実な実装)\"\"\"
    # Rではrow-byrow(byrow=T)で行列を構築
    # 各行: [off-diag要素をlambdaの順番に入れ、対角要素は1-行の残り]
    # lam[0:5]  → 行0: diag, lam[0], lam[1], lam[2], lam[3], lam[4]
    # lam[5:10] → 行1: lam[5], diag, lam[6], lam[7], lam[8], lam[9]
    # ...
    mat = np.zeros((num_price_states, num_price_states))
    idx = 0
    for i in range(num_price_states):
        row_sum = 0
        for j in range(num_price_states):
            if i != j:
                mat[i, j] = lam[idx]
                row_sum += lam[idx]
                idx += 1
        mat[i, i] = 1 - row_sum
    return mat


def flow_utility(theta, state_df_local):
    \"\"\"状態変数・選択肢毎の今期の効用を返す\"\"\"
    theta_c, theta_p = theta
    U_not_buy = -theta_c * state_df_local['mileage'].values
    U_buy = -theta_p * state_df_local['price'].values
    return np.column_stack([U_not_buy, U_buy])


def contraction(theta, beta, trans_mat, state_df_local):
    \"\"\"価値関数反復法（縮小写像）\"\"\"
    n = len(state_df_local)
    V_old = np.zeros(n)
    U = flow_utility(theta, state_df_local)

    diff = 1000
    tol_level = 1.0e-12

    while diff > tol_level:
        EV = np.column_stack([
            trans_mat['not_buy'] @ V_old,
            trans_mat['buy'] @ V_old
        ])
        V_new = np.log(np.sum(np.exp(U + beta * EV), axis=1)) + Euler_const
        diff = np.max(np.abs(V_new - V_old))
        V_old = V_new.copy()

    return V_old


print("関数定義完了")
"""))

# ============================================================
# Cell 4: Transition Matrices
# ============================================================
cells.append(nbf.v4.new_markdown_cell("## 遷移行列の作成"))

cells.append(nbf.v4.new_code_cell("""\
# 走行距離の遷移行列パラメータ
kappa_true = np.array([0.4, 0.1])
mileage_trans_mat_true = gen_mileage_trans(kappa_true)

print("走行距離の遷移行列（購買なし、左上4x4）:")
print(mileage_trans_mat_true[:4, :4, 0])

# 価格の遷移行列パラメータ
lambda_true = np.array([
    0.1, 0.2, 0.2, 0.2, 0.2,
    0.1, 0.2, 0.2, 0.2, 0.2,
    0.1, 0.1, 0.2, 0.2, 0.1,
    0.1, 0.1, 0.2, 0.2, 0.1,
    0.05, 0.05, 0.1, 0.1, 0.2,
    0.05, 0.05, 0.1, 0.1, 0.2
])
price_trans_mat_true = gen_price_trans(lambda_true)

print("\\n価格の遷移行列:")
print(price_trans_mat_true)

# コントロール変数ごとの遷移行列（クロネッカー積）
trans_mat_true = {
    'not_buy': np.kron(mileage_trans_mat_true[:, :, 0], price_trans_mat_true),
    'buy': np.kron(mileage_trans_mat_true[:, :, 1], price_trans_mat_true)
}

# 定常状態での価格の分布を計算
eigenvalues, eigenvectors = np.linalg.eig(price_trans_mat_true.T)
# 固有値が1に最も近いものを選択
idx = np.argmin(np.abs(eigenvalues - 1))
price_dist_steady = np.real(eigenvectors[:, idx])
price_dist_steady = price_dist_steady / price_dist_steady.sum()

print("\\n価格の定常分布:")
print(price_dist_steady)
"""))

# ============================================================
# Cell 5: Value Function
# ============================================================
cells.append(nbf.v4.new_markdown_cell("## 価値関数の計算"))

cells.append(nbf.v4.new_code_cell("""\
import time

start_time = time.time()
V_true = contraction(theta_true, beta, trans_mat_true, state_df)
elapsed = time.time() - start_time
print(f"Runtime: {elapsed:.3f} sec")

# 選択毎の価値関数
U_true = flow_utility(theta_true, state_df)
EV_true = np.column_stack([
    trans_mat_true['not_buy'] @ V_true,
    trans_mat_true['buy'] @ V_true
])
V_CS_true = U_true + beta * EV_true

# logitによる理論上の条件付き購入確率
prob_buy_true = np.exp(V_CS_true[:, 1]) / np.sum(np.exp(V_CS_true), axis=1)
prob_buy_true_mat = prob_buy_true.reshape(num_mileage_states, num_price_states).T

print("\\n理論上の購入確率（price x mileage 行列）:")
print(prob_buy_true_mat)
"""))

# ============================================================
# Cell 6: CCP Visualization
# ============================================================
cells.append(nbf.v4.new_markdown_cell("## 条件付き購入確率（CCP）の可視化"))

cells.append(nbf.v4.new_code_cell("""\
# 3次元プロット
fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection='3d')

X, Y = np.meshgrid(mileage_states * 100, price_states * 100)
ax.plot_surface(X, Y, prob_buy_true_mat, alpha=0.7, color='grey', edgecolor='black', linewidth=0.3)
ax.set_xlabel('Mileage (10k km)')
ax.set_ylabel('Price (10k JPY)')
ax.set_zlabel('Purchase Prob.')
ax.set_zlim(0, 0.8)
ax.view_init(elev=10, azim=-60)

plt.tight_layout()
plt.savefig(output_dir / 'CCP_true_3D.png', dpi=150)
plt.show()
"""))

cells.append(nbf.v4.new_code_cell("""\
# 2次元プロット: 価格ごとの購入確率
fig, ax = plt.subplots(figsize=(8, 5))

for i, p in enumerate(price_states):
    ax.plot(mileage_states * 100, prob_buy_true_mat[i, :],
            marker='o', markersize=4,
            label=f'{int(p * 100)}万円',
            color=plt.cm.Greys(0.3 + 0.1 * i))

ax.set_xlabel('走行距離（万km）')
ax.set_ylabel('購買確率')
ax.set_xticks(np.arange(0, 11, 2))
ax.set_ylim(bottom=0)
ax.legend(title='価格')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig(output_dir / 'fig6_3_CCP_true_2D.png', dpi=150)
plt.show()
"""))

# ============================================================
# Cell 7: Data Simulation
# ============================================================
cells.append(nbf.v4.new_markdown_cell("## シミュレーション"))

cells.append(nbf.v4.new_code_cell("""\
# サンプルサイズ
num_consumer = 1000
num_period = 50
num_period_obs = 10
num_obs = num_consumer * num_period

# 累積遷移行列
trans_mat_cum = {
    'not_buy': np.cumsum(trans_mat_true['not_buy'], axis=1),
    'buy': np.cumsum(trans_mat_true['buy'], axis=1)
}

# 乱数生成（Rのset.seed(1)とは異なるため独自シード）
rng = np.random.RandomState(1)

# Gumbel(0,1) 乱数 = Type-I Extreme Value
eps_type1_not_buy = stats.gumbel_r.rvs(size=num_obs, random_state=rng)
eps_type1_buy = stats.gumbel_r.rvs(size=num_obs, random_state=rng)
eps_unif = rng.uniform(size=num_obs)
eps_price_state_unif = rng.uniform(size=num_obs)

# データ生成
consumer_ids = np.repeat(np.arange(1, num_consumer + 1), num_period)
period_ids = np.tile(np.arange(1, num_period + 1), num_consumer)
state_ids = np.zeros(num_obs, dtype=int)
actions = np.zeros(num_obs, dtype=int)

# 価格定常分布の累積値
price_dist_steady_cumsum = np.cumsum(price_dist_steady)

# state_dfの値をnumpyで高速にアクセスするための準備
state_mileage_id = state_df['mileage_id'].values
state_price_id = state_df['price_id'].values

print("データ生成中...")
start_time = time.time()

for c in range(num_consumer):
    base = c * num_period

    # 初期の価格状態を定常分布から決定
    price_id_consumer = np.searchsorted(price_dist_steady_cumsum,
                                         eps_price_state_unif[base]) + 1
    price_id_consumer = min(price_id_consumer, num_price_states)

    # 初期のstate_id（mileage_id=1）
    init_state = (price_id_consumer - 1) + (1 - 1) * num_price_states + 1
    state_ids[base] = init_state

    for t in range(num_period - 1):
        idx = base + t
        sid = state_ids[idx] - 1  # 0-indexed

        # 購入判定
        v_not_buy = V_CS_true[sid, 0] + eps_type1_not_buy[idx]
        v_buy = V_CS_true[sid, 1] + eps_type1_buy[idx]

        if v_not_buy > v_buy:
            actions[idx] = 0
            cum_row = trans_mat_cum['not_buy'][sid, :]
        else:
            actions[idx] = 1
            cum_row = trans_mat_cum['buy'][sid, :]

        # 次期の状態
        next_state = np.searchsorted(cum_row, eps_unif[idx]) + 1
        next_state = min(next_state, num_states)
        state_ids[base + t + 1] = next_state

elapsed = time.time() - start_time
print(f"データ生成完了: {elapsed:.1f} sec")

# DataFrame作成
data_gen = pd.DataFrame({
    'consumer': consumer_ids,
    'period': period_ids,
    'eps_type1_not_buy': eps_type1_not_buy,
    'eps_type1_buy': eps_type1_buy,
    'eps_unif': eps_unif,
    'eps_price_state_unif': eps_price_state_unif,
    'state_id': state_ids,
    'action': actions
})

# 最後の10年のみ観察
data_gen = data_gen[data_gen['period'] > (num_period - num_period_obs)].copy()

# 状態変数を結合
data_gen = data_gen.merge(state_df, on='state_id', how='left')

print(f"観察数: {len(data_gen)}")
data_gen.tail(3)
"""))

# ============================================================
# Cell 8: Save Data & Descriptive Stats
# ============================================================
cells.append(nbf.v4.new_markdown_cell("## データの保存と記述統計"))

cells.append(nbf.v4.new_code_cell("""\
# CSVとして保存
data_gen.to_csv(intermediate_dir / 'Chap6_data.csv', index=True)
print("データを保存しました: intermediate/Chap6_data.csv")

# 記述統計
desc_stats = data_gen[['price', 'mileage', 'action']].describe()
print("\\n記述統計:")
print(desc_stats.loc[['mean', 'std', 'min', 'max']].to_string())
"""))

# ============================================================
# Cell 9: Distribution Plots
# ============================================================
cells.append(nbf.v4.new_markdown_cell("## データの分布"))

cells.append(nbf.v4.new_code_cell("""\
fig, axes = plt.subplots(1, 2, figsize=(12, 4))

# 走行距離の分布
axes[0].hist(data_gen['mileage'] * 100, bins=21, color='grey', edgecolor='black')
axes[0].set_xlabel('走行距離（万km）')
axes[0].set_ylabel('頻度')
axes[0].set_xticks(np.arange(0, 11, 2))
axes[0].spines['top'].set_visible(False)
axes[0].spines['right'].set_visible(False)

# 価格の分布
axes[1].hist(data_gen['price'] * 100, bins=6, color='grey', edgecolor='black')
axes[1].set_xlabel('価格（万円）')
axes[1].set_ylabel('頻度')
axes[1].set_xticks(np.arange(200, 260, 10))
axes[1].spines['top'].set_visible(False)
axes[1].spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig(output_dir / 'dist_mile_price.png', dpi=150)
plt.show()
"""))

# ============================================================
# Cell 10: Observed CCP
# ============================================================
cells.append(nbf.v4.new_markdown_cell("## 観測された条件付き購入確率（CCP）"))

cells.append(nbf.v4.new_code_cell("""\
fig, axes = plt.subplots(1, 2, figsize=(12, 4))

# 走行距離ごとの購入確率
ccp_mile = data_gen.groupby('mileage').agg(
    num_state=('action', 'count'),
    sum_action=('action', 'sum')
).reset_index()
ccp_mile['prob_buy'] = ccp_mile['sum_action'] / ccp_mile['num_state']

axes[0].bar(ccp_mile['mileage'] * 100, ccp_mile['prob_buy'],
            color='grey', edgecolor='black', width=0.4)
axes[0].set_xlabel('走行距離（万km）')
axes[0].set_ylabel('購入確率')
axes[0].set_xticks(np.arange(0, 11, 2))
axes[0].spines['top'].set_visible(False)
axes[0].spines['right'].set_visible(False)

# 価格ごとの購入確率
ccp_price = data_gen.groupby('price').agg(
    num_state=('action', 'count'),
    sum_action=('action', 'sum')
).reset_index()
ccp_price['prob_buy'] = ccp_price['sum_action'] / ccp_price['num_state']

axes[1].bar(ccp_price['price'] * 100, ccp_price['prob_buy'],
            color='grey', edgecolor='black', width=8)
axes[1].set_xlabel('価格（万円）')
axes[1].set_ylabel('購入確率')
axes[1].set_xticks(np.arange(200, 260, 10))
axes[1].spines['top'].set_visible(False)
axes[1].spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig(output_dir / 'CCP_mile_price.png', dpi=150)
plt.show()
"""))

cells.append(nbf.v4.new_code_cell("""\
# state(p,m)ごとに観測された条件付き購入確率
prob_buy_obs = data_gen.groupby(['mileage', 'price']).agg(
    num_state=('action', 'count'),
    sum_action=('action', 'sum')
).reset_index()
prob_buy_obs['prob_buy'] = prob_buy_obs['sum_action'] / prob_buy_obs['num_state']
prob_buy_obs_mat = prob_buy_obs['prob_buy'].values.reshape(num_mileage_states, num_price_states).T

# 3Dプロット
fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection='3d')

X, Y = np.meshgrid(mileage_states * 100, price_states * 100)
ax.plot_surface(X, Y, prob_buy_obs_mat, alpha=0.7, color='grey', edgecolor='black', linewidth=0.3)
ax.set_xlabel('走行距離（万km）')
ax.set_ylabel('価格（万円）')
ax.set_zlabel('購入確率')
ax.set_zlim(0, 0.8)
ax.view_init(elev=10, azim=-60)

plt.tight_layout()
plt.savefig(output_dir / 'CCP_3D.png', dpi=150)
plt.show()
"""))

# ============================================================
# Cell 11: Transition Matrix Estimation
# ============================================================
cells.append(nbf.v4.new_markdown_cell("## 遷移行列の推定"))

cells.append(nbf.v4.new_code_cell("""\
# ラグ変数を追加
data_gen = data_gen.sort_values(['consumer', 'period']).copy()
data_gen['lag_price_id'] = data_gen.groupby('consumer')['price_id'].shift(1)
data_gen['lag_mileage_id'] = data_gen.groupby('consumer')['mileage_id'].shift(1)
data_gen['lag_action'] = data_gen.groupby('consumer')['action'].shift(1)

# === 走行距離の遷移行列の推定 ===
# 1期目を除外
df_est = data_gen.dropna(subset=['lag_action']).copy()
df_est['lag_price_id'] = df_est['lag_price_id'].astype(int)
df_est['lag_mileage_id'] = df_est['lag_mileage_id'].astype(int)
df_est['lag_action'] = df_est['lag_action'].astype(int)

# 確率カテゴリを判定
def classify_mileage_obs(row):
    lag_m = row['lag_mileage_id']
    m = row['mileage_id']
    lag_a = row['lag_action']

    # 1 - kappa_1 - kappa_2
    if ((lag_a == 0 and 1 <= lag_m <= 20 and lag_m == m) or
        (lag_a == 1 and m == 1)):
        return 'cond1'
    # kappa_1
    elif ((lag_a == 0 and 1 <= lag_m <= 19 and lag_m == m - 1) or
          (lag_a == 1 and m == 2)):
        return 'cond2'
    # kappa_2
    elif ((lag_a == 0 and 1 <= lag_m <= 19 and lag_m == m - 2) or
          (lag_a == 1 and m == 3)):
        return 'cond3'
    # kappa_1 + kappa_2
    elif (lag_a == 0 and lag_m == 20 and m == 21):
        return 'cond4'
    else:
        return 'other'

df_est['cond_obs_mileage'] = df_est.apply(classify_mileage_obs, axis=1)

num_cond_obs_mileage = (
    df_est[df_est['cond_obs_mileage'] != 'other']
    .groupby('cond_obs_mileage')
    .size()
    .sort_index()
    .values
    .astype(float)
)

print("走行距離の遷移カテゴリ別観察数:", num_cond_obs_mileage)

# 最尤法の解析解
n1, n2, n3, n4 = num_cond_obs_mileage
kappa_est = np.array([
    n2 * (n2 + n3 + n4) / ((n2 + n3) * (n1 + n2 + n3 + n4)),
    n3 * (n2 + n3 + n4) / ((n2 + n3) * (n1 + n2 + n3 + n4))
])
print(f"\\nkappa推定値: {kappa_est}")

# 標準誤差（フィッシャー情報量から）
I_mat = np.zeros((2, 2))
I_mat[0, 0] = n1 / (1 - kappa_est[0] - kappa_est[1])**2 + n2 / kappa_est[0]**2 + n4 / (kappa_est[0] + kappa_est[1])**2
I_mat[0, 1] = n1 / (1 - kappa_est[0] - kappa_est[1])**2 + n4 / (kappa_est[0] + kappa_est[1])**2
I_mat[1, 0] = I_mat[0, 1]
I_mat[1, 1] = n1 / (1 - kappa_est[0] - kappa_est[1])**2 + n3 / kappa_est[1]**2 + n4 / (kappa_est[0] + kappa_est[1])**2

kappa_se = np.sqrt(np.diag(np.linalg.inv(I_mat)))
print(f"kappa標準誤差: {kappa_se}")
print(f"kappa真値: {kappa_true}")
"""))

cells.append(nbf.v4.new_code_cell("""\
# === 価格の遷移行列の推定 ===
num_cond_obs_price = pd.crosstab(
    df_est['lag_price_id'], df_est['price_id']
).values.astype(float)

lambda_est_mat = num_cond_obs_price / num_cond_obs_price.sum(axis=1, keepdims=True)
print("価格の遷移行列推定値:")
print(lambda_est_mat)

# 対角要素を除いて1次元ベクトルに変換
lambda_est = []
for i in range(num_price_states):
    for j in range(num_price_states):
        if i != j:
            lambda_est.append(lambda_est_mat[i, j])
lambda_est = np.array(lambda_est)

# 標準誤差
lambda_se = []
for i in range(num_price_states):
    n_row = num_cond_obs_price[i, :]
    diag_part = np.diag(n_row)
    # i行i列を除いた (num_price_states-1) x (num_price_states-1) 行列
    idx = [j for j in range(num_price_states) if j != i]
    info_mat = (diag_part[np.ix_(idx, idx)] / lambda_est_mat[np.ix_([i], idx)].flatten()**2 +
                n_row[i] / lambda_est_mat[i, i]**2)
    se_i = np.sqrt(np.diag(np.linalg.inv(info_mat)))
    lambda_se.extend(se_i)

lambda_se = np.array(lambda_se)

print(f"\\nlambda推定値 (off-diag): {lambda_est[:5]} ...")
print(f"lambda標準誤差: {lambda_se[:5]} ...")
"""))

# ============================================================
# Cell 12: Static Logit
# ============================================================
cells.append(nbf.v4.new_markdown_cell("## 静学的ロジットによる推定"))

cells.append(nbf.v4.new_code_cell("""\
def logLH_stat(theta, state_df_local, df):
    \"\"\"静学的ロジットの対数尤度関数\"\"\"
    U = flow_utility(theta, state_df_local)
    prob_C = np.exp(U) / np.sum(np.exp(U), axis=1, keepdims=True)

    # 各観察の対数尤度
    state_idx = df['state_id'].values - 1  # 0-indexed
    action_idx = df['action'].values
    log_lik = np.sum(np.log(prob_C[state_idx, action_idx]))
    return log_lik

# 最適化（最大化 → fnscale=-1 に対応）
start_time = time.time()

result_stat = optimize.minimize(
    lambda theta: -logLH_stat(theta, state_df, data_gen),
    theta_true,
    method='Nelder-Mead'
)

elapsed = time.time() - start_time
theta_est_stat = result_stat.x
print(f"Runtime: {elapsed:.3f} sec")
print(f"静学ロジット推定値: theta_c={theta_est_stat[0]:.2f}, theta_p={theta_est_stat[1]:.2f}")

# 標準誤差（数値ヘシアン）
def neg_logLH_stat(theta):
    return -logLH_stat(theta, state_df, data_gen)

eps = 1e-5
n_params = len(theta_est_stat)
hessian_stat = np.zeros((n_params, n_params))
f0 = neg_logLH_stat(theta_est_stat)
for i in range(n_params):
    for j in range(n_params):
        e_i = np.zeros(n_params); e_i[i] = eps
        e_j = np.zeros(n_params); e_j[j] = eps
        hessian_stat[i, j] = (neg_logLH_stat(theta_est_stat + e_i + e_j)
                               - neg_logLH_stat(theta_est_stat + e_i)
                               - neg_logLH_stat(theta_est_stat + e_j)
                               + f0) / eps**2

theta_se_stat = np.sqrt(np.diag(np.linalg.inv(hessian_stat)))
print(f"静学ロジット標準誤差: theta_c={theta_se_stat[0]:.2f}, theta_p={theta_se_stat[1]:.2f}")
"""))

# ============================================================
# Cell 13: NFXP Estimation
# ============================================================
cells.append(nbf.v4.new_markdown_cell("## 入れ子不動点アルゴリズム（NFXP）による推定"))

cells.append(nbf.v4.new_code_cell("""\
def logLH_nfxp(theta, beta, trans_mat, state_df_local, df):
    \"\"\"NFXP対数尤度関数\"\"\"
    U = flow_utility(theta, state_df_local)
    V = contraction(theta, beta, trans_mat, state_df_local)
    EV = np.column_stack([
        trans_mat['not_buy'] @ V,
        trans_mat['buy'] @ V
    ])
    V_CS = U + beta * EV
    prob_C = np.exp(V_CS) / np.sum(np.exp(V_CS), axis=1, keepdims=True)

    state_idx = df['state_id'].values - 1
    action_idx = df['action'].values
    return np.sum(np.log(prob_C[state_idx, action_idx]))

# 推定された遷移行列を使用
trans_mat_hat = {
    'not_buy': np.kron(gen_mileage_trans(kappa_est)[:, :, 0], gen_price_trans(lambda_est)),
    'buy': np.kron(gen_mileage_trans(kappa_est)[:, :, 1], gen_price_trans(lambda_est))
}

start_time = time.time()

result_nfxp = optimize.minimize(
    lambda theta: -logLH_nfxp(theta, beta, trans_mat_hat, state_df, data_gen),
    theta_true,
    method='Nelder-Mead'
)

elapsed = time.time() - start_time
theta_est_nfxp = result_nfxp.x
print(f"Runtime: {elapsed:.3f} sec")
print(f"NFXP推定値: theta_c={theta_est_nfxp[0]:.2f}, theta_p={theta_est_nfxp[1]:.2f}")

# 標準誤差
def neg_logLH_nfxp(theta):
    return -logLH_nfxp(theta, beta, trans_mat_hat, state_df, data_gen)

hessian_nfxp = np.zeros((n_params, n_params))
f0 = neg_logLH_nfxp(theta_est_nfxp)
for i in range(n_params):
    for j in range(n_params):
        e_i = np.zeros(n_params); e_i[i] = eps
        e_j = np.zeros(n_params); e_j[j] = eps
        hessian_nfxp[i, j] = (neg_logLH_nfxp(theta_est_nfxp + e_i + e_j)
                               - neg_logLH_nfxp(theta_est_nfxp + e_i)
                               - neg_logLH_nfxp(theta_est_nfxp + e_j)
                               + f0) / eps**2

theta_se_nfxp = np.sqrt(np.diag(np.linalg.inv(hessian_nfxp)))
print(f"NFXP標準誤差: theta_c={theta_se_nfxp[0]:.2f}, theta_p={theta_se_nfxp[1]:.2f}")
"""))

# ============================================================
# Cell 14: Comparison Table
# ============================================================
cells.append(nbf.v4.new_markdown_cell("## 推定結果の比較"))

cells.append(nbf.v4.new_code_cell("""\
results_ch6 = pd.DataFrame({
    'algorithm': ['Static', 'NFXP', 'True'],
    'theta_c': [theta_est_stat[0], theta_est_nfxp[0], theta_true[0]],
    'theta_se_c': [theta_se_stat[0], theta_se_nfxp[0], np.nan],
    'theta_p': [theta_est_stat[1], theta_est_nfxp[1], theta_true[1]],
    'theta_se_p': [theta_se_stat[1], theta_se_nfxp[1], np.nan]
})

print("推定結果の比較:")
print(results_ch6.to_string(index=False))

# テキストファイルとして保存
with open(output_dir / 'tab6_3_compare_est.txt', 'w') as f:
    f.write(results_ch6.to_string(index=False))
"""))

nb.cells = cells
nbf.write(nb, 'main_ch06.ipynb')
print("Generated: main_ch06.ipynb")
