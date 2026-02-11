"""Generate main_ch10_equilibrium.ipynb for Chapter 10: Dynamic Game - Equilibrium Computation."""
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
    "# 第10章 動的ゲーム: 均衡CCPの計算と疑似データの生成\n\n"
    "2企業の参入・退出モデルにおけるマルコフ完全均衡（MPE）を計算し、\n"
    "均衡CCPに基づいて疑似データを生成する。"
))

# ============================================================
# Setup
# ============================================================
cells.append(nbf.v4.new_markdown_cell("## 1. Pythonに関する下準備"))

cells.append(nbf.v4.new_code_cell("""\
import numpy as np
import pandas as pd
from pathlib import Path
import warnings

warnings.filterwarnings('ignore')

# パス設定
base_dir = Path('..')
data_dir = base_dir / 'data_from_matlab'
output_dir = base_dir / 'output'
output_dir.mkdir(exist_ok=True)

print("Setup complete.")
"""))

# ============================================================
# Model Primitives
# ============================================================
cells.append(nbf.v4.new_markdown_cell(
    "## 2. モデルの基本設定\n\n"
    "### 状態空間\n"
    "状態変数は (景気, 企業1の出店状況, 企業2の出店状況) の組み合わせで8通り:\n\n"
    "| 状態 (0-indexed) | 景気 | n1 | n2 | R での状態番号 |\n"
    "|:---:|:---:|:---:|:---:|:---:|\n"
    "| 0 | Good (1) | 0 | 0 | 1 |\n"
    "| 1 | Good (1) | 0 | 1 | 2 |\n"
    "| 2 | Good (1) | 1 | 0 | 3 |\n"
    "| 3 | Good (1) | 1 | 1 | 4 |\n"
    "| 4 | Bad (2)  | 0 | 0 | 5 |\n"
    "| 5 | Bad (2)  | 0 | 1 | 6 |\n"
    "| 6 | Bad (2)  | 1 | 0 | 7 |\n"
    "| 7 | Bad (2)  | 1 | 1 | 8 |\n\n"
    "### 行動空間\n"
    "各企業は 3 つの行動を取りうる（列インデックスで管理）:\n"
    "- 列 0: a = -1 (退出)\n"
    "- 列 1: a = 0  (現状維持)\n"
    "- 列 2: a = +1 (参入/投資)\n\n"
    "ただし、出店していない企業は退出できず、出店済みの企業は参入できない。\n"
    "CCP Adjuster行列で実行不可能な行動をマスクする。"
))

cells.append(nbf.v4.new_code_cell("""\
# --- 基本パラメータ ---
beta = 0.8                # 割引因子
euler_gamma = 0.5772      # オイラー・マスケローニ定数

# 景気の遷移行列
TransitionMat = np.array([
    [0.7, 0.3],
    [0.4, 0.6]
])

# パラメータの設定
# Parameters = [企業1ベース利潤, 企業2ベース利潤, 顧客収奪効果,
#               好景気追加利潤, 退出コスト, 参入コスト]
Parameters = np.array([0.3, 0.2, -0.27, 0.45, -0.15, -2.10])

# パラメータを並べ替え (theta: 10 x 1)
# theta[0]: 企業1のベース利潤
# theta[1]: ライバルが企業1の利潤に与える影響
# theta[2]: 好景気の追加的利潤 (企業1)
# theta[3]: 企業1の退出コスト
# theta[4]: 企業1の参入コスト
# theta[5]: 企業2のベース利潤
# theta[6]: ライバルが企業2の利潤に与える影響
# theta[7]: 好景気の追加的利潤 (企業2)
# theta[8]: 企業2の退出コスト
# theta[9]: 企業2の参入コスト
TrueParameterValues = np.array([
    Parameters[0],   # 企業1のベース利潤
    Parameters[2],   # 顧客収奪効果 (企業1)
    Parameters[3],   # 好景気の追加利潤 (企業1)
    Parameters[4],   # 退出コスト (企業1)
    Parameters[5],   # 参入コスト (企業1)
    Parameters[1],   # 企業2のベース利潤
    Parameters[2],   # 顧客収奪効果 (企業2)
    Parameters[3],   # 好景気の追加利潤 (企業2)
    Parameters[4],   # 退出コスト (企業2)
    Parameters[5],   # 参入コスト (企業2)
])

print("TrueParameterValues:")
labels = ['theta_1 (Firm1 base profit)', 'theta_2 (Rival effect on Firm1)',
          'theta_3 (Good economy, Firm1)', 'theta_4 (Exit cost, Firm1)',
          'theta_5 (Entry cost, Firm1)', 'theta_6 (Firm2 base profit)',
          'theta_7 (Rival effect on Firm2)', 'theta_8 (Good economy, Firm2)',
          'theta_9 (Exit cost, Firm2)', 'theta_10 (Entry cost, Firm2)']
for i, (label, val) in enumerate(zip(labels, TrueParameterValues)):
    print(f"  [{i}] {label}: {val:.2f}")
"""))

# ============================================================
# CCP Adjuster
# ============================================================
cells.append(nbf.v4.new_markdown_cell("### CCP Adjuster 行列"))

cells.append(nbf.v4.new_code_cell("""\
# CCP Adjuster: 実行可能な行動に1、不可能な行動に0を割り振る (8 x 3)
# 列: [a=-1(退出), a=0(現状維持), a=+1(参入)]
#
# 企業1について:
#   状態0 (G,0,0): 出店なし -> 退出不可 -> [0, 1, 1]
#   状態1 (G,0,1): 出店なし -> 退出不可 -> [0, 1, 1]
#   状態2 (G,1,0): 出店あり -> 参入不可 -> [1, 1, 0]
#   状態3 (G,1,1): 出店あり -> 参入不可 -> [1, 1, 0]
#   (Bad景気でも同様のパターン)
CCP1Adjuster = np.array([
    [0, 1, 1],  # s=0: G,0,0
    [0, 1, 1],  # s=1: G,0,1
    [1, 1, 0],  # s=2: G,1,0
    [1, 1, 0],  # s=3: G,1,1
    [0, 1, 1],  # s=4: B,0,0
    [0, 1, 1],  # s=5: B,0,1
    [1, 1, 0],  # s=6: B,1,0
    [1, 1, 0],  # s=7: B,1,1
], dtype=float)

# 企業2について:
#   状態0 (G,0,0): 出店なし -> 退出不可 -> [0, 1, 1]
#   状態1 (G,0,1): 出店あり -> 参入不可 -> [1, 1, 0]
#   状態2 (G,1,0): 出店なし -> 退出不可 -> [0, 1, 1]
#   状態3 (G,1,1): 出店あり -> 参入不可 -> [1, 1, 0]
#   (Bad景気でも同様のパターン)
CCP2Adjuster = np.array([
    [0, 1, 1],  # s=0: G,0,0
    [1, 1, 0],  # s=1: G,0,1
    [0, 1, 1],  # s=2: G,1,0
    [1, 1, 0],  # s=3: G,1,1
    [0, 1, 1],  # s=4: B,0,0
    [1, 1, 0],  # s=5: B,0,1
    [0, 1, 1],  # s=6: B,1,0
    [1, 1, 0],  # s=7: B,1,1
], dtype=float)

print("CCP1Adjuster:")
print(CCP1Adjuster)
print("\\nCCP2Adjuster:")
print(CCP2Adjuster)
"""))

# ============================================================
# Profit Functions
# ============================================================
cells.append(nbf.v4.new_markdown_cell("### 利潤関数の定義"))

cells.append(nbf.v4.new_code_cell("""\
def pi1gen(theta):
    \"\"\"企業1の利潤行列を生成する。

    Parameters
    ----------
    theta : array (10,)
        パラメータベクトル

    Returns
    -------
    output : array (8, 3)
        状態(行) x 行動(列) の利潤行列
        列: [a=-1(退出), a=0(現状維持), a=+1(参入)]

    状態ごとのベース利潤 (8 x 1):
        s=0 (G,0,0): 0 (出店していないので利潤なし)
        s=1 (G,0,1): 0
        s=2 (G,1,0): theta[0] + theta[2]  (ベース + 好景気)
        s=3 (G,1,1): theta[0] + theta[1] + theta[2]  (ベース + 競合影響 + 好景気)
        s=4 (B,0,0): 0
        s=5 (B,0,1): 0
        s=6 (B,1,0): theta[0]  (ベースのみ)
        s=7 (B,1,1): theta[0] + theta[1]  (ベース + 競合影響)

    投資/退出コスト (1 x 3):
        [theta[3](退出コスト), 0(現状維持), theta[4](参入コスト)]
    \"\"\"
    base = np.array([
        0,
        0,
        theta[0] + theta[2],
        theta[0] + theta[1] + theta[2],
        0,
        0,
        theta[0],
        theta[0] + theta[1]
    ])  # (8,)

    invdiv = np.array([theta[3], 0, theta[4]])  # (3,)

    # base を (8, 3) に拡張し、invdiv を (8, 3) に拡張して加算
    output = base[:, np.newaxis] + invdiv[np.newaxis, :]

    return output


def pi2gen(theta):
    \"\"\"企業2の利潤行列を生成する。

    Parameters
    ----------
    theta : array (10,)
        パラメータベクトル

    Returns
    -------
    output : array (8, 3)
        状態(行) x 行動(列) の利潤行列
    \"\"\"
    base = np.array([
        0,
        theta[5] + theta[7],
        0,
        theta[5] + theta[6] + theta[7],
        0,
        theta[5],
        0,
        theta[5] + theta[6]
    ])  # (8,)

    invdiv = np.array([theta[8], 0, theta[9]])  # (3,)

    output = base[:, np.newaxis] + invdiv[np.newaxis, :]

    return output


# 真のパラメータの下で利潤行列を計算し、Adjusterを適用
pi1 = pi1gen(TrueParameterValues) * CCP1Adjuster
pi2 = pi2gen(TrueParameterValues) * CCP2Adjuster

print("pi1 (企業1の利潤行列, Adjuster適用後):")
print(np.round(pi1, 4))
print("\\npi2 (企業2の利潤行列, Adjuster適用後):")
print(np.round(pi2, 4))
"""))

# ============================================================
# CCP Transform Functions
# ============================================================
cells.append(nbf.v4.new_markdown_cell("### CCP変換関数の定義"))

cells.append(nbf.v4.new_code_cell("""\
def CCP1Transform(x):
    \"\"\"CCPベクトル (8,) からCCP行列 (8, 3) への変換（企業1）。

    x[i] は状態 i において企業1が「現状維持 (a=0)」を選ぶ確率。
    出店していない状態 (n1=0) では: [0, x[i], 1-x[i]]  (退出不可)
    出店している状態 (n1=1) では:   [1-x[i], x[i], 0]  (参入不可)
    \"\"\"
    return np.array([
        [0,        x[0], 1 - x[0]],   # s=0: G,0,0
        [0,        x[1], 1 - x[1]],   # s=1: G,0,1
        [1 - x[2], x[2], 0       ],   # s=2: G,1,0
        [1 - x[3], x[3], 0       ],   # s=3: G,1,1
        [0,        x[4], 1 - x[4]],   # s=4: B,0,0
        [0,        x[5], 1 - x[5]],   # s=5: B,0,1
        [1 - x[6], x[6], 0       ],   # s=6: B,1,0
        [1 - x[7], x[7], 0       ],   # s=7: B,1,1
    ])


def CCP2Transform(x):
    \"\"\"CCPベクトル (8,) からCCP行列 (8, 3) への変換（企業2）。

    x[i] は状態 i において企業2が「現状維持 (a=0)」を選ぶ確率。
    出店していない状態 (n2=0) では: [0, x[i], 1-x[i]]  (退出不可)
    出店している状態 (n2=1) では:   [1-x[i], x[i], 0]  (参入不可)
    \"\"\"
    return np.array([
        [0,        x[0], 1 - x[0]],   # s=0: G,0,0
        [1 - x[1], x[1], 0       ],   # s=1: G,0,1
        [0,        x[2], 1 - x[2]],   # s=2: G,1,0
        [1 - x[3], x[3], 0       ],   # s=3: G,1,1
        [0,        x[4], 1 - x[4]],   # s=4: B,0,0
        [1 - x[5], x[5], 0       ],   # s=5: B,0,1
        [0,        x[6], 1 - x[6]],   # s=6: B,1,0
        [1 - x[7], x[7], 0       ],   # s=7: B,1,1
    ])


def CCP1LogTransform(x):
    \"\"\"CCPベクトル (8,) から対数CCP行列 (8, 3) への変換（企業1）。\"\"\"
    return np.array([
        [0,             np.log(x[0]), np.log(1 - x[0])],
        [0,             np.log(x[1]), np.log(1 - x[1])],
        [np.log(1 - x[2]), np.log(x[2]), 0             ],
        [np.log(1 - x[3]), np.log(x[3]), 0             ],
        [0,             np.log(x[4]), np.log(1 - x[4])],
        [0,             np.log(x[5]), np.log(1 - x[5])],
        [np.log(1 - x[6]), np.log(x[6]), 0             ],
        [np.log(1 - x[7]), np.log(x[7]), 0             ],
    ])


def CCP2LogTransform(x):
    \"\"\"CCPベクトル (8,) から対数CCP行列 (8, 3) への変換（企業2）。\"\"\"
    return np.array([
        [0,             np.log(x[0]), np.log(1 - x[0])],
        [np.log(1 - x[1]), np.log(x[1]), 0             ],
        [0,             np.log(x[2]), np.log(1 - x[2])],
        [np.log(1 - x[3]), np.log(x[3]), 0             ],
        [0,             np.log(x[4]), np.log(1 - x[4])],
        [np.log(1 - x[5]), np.log(x[5]), 0             ],
        [0,             np.log(x[6]), np.log(1 - x[6])],
        [np.log(1 - x[7]), np.log(x[7]), 0             ],
    ])


# 動作確認
test_ccp = np.full(8, 0.5)
print("CCP1Transform(0.5):")
print(CCP1Transform(test_ccp))
print("\\nCCP2Transform(0.5):")
print(CCP2Transform(test_ccp))
"""))

# ============================================================
# Transition Matrix Functions
# ============================================================
cells.append(nbf.v4.new_markdown_cell("### 遷移行列の構築"))

cells.append(nbf.v4.new_code_cell("""\
def fP(trans_mat, ccp1_vec, ccp2_vec):
    \"\"\"CCPが与えられたときの状態遷移行列 F^{P,sigma} を構築する (8 x 8)。

    3つの行列のアダマール積（要素ごとの積）として構成される:
    1. TempMat0: 外生状態（景気）の遷移確率
    2. TempMat1: 企業1のCCPに基づく遷移確率
    3. TempMat2: 企業2のCCPに基づく遷移確率
    \"\"\"
    # 外生状態の遷移確率を 8x8 に拡張
    TempMat0 = np.kron(trans_mat, np.ones((4, 4)))

    # 企業1の行動による遷移
    # ccp1_vec[i] = 状態iで企業1がa=0を選ぶ確率
    # n1=0のとき: 確率 ccp1_vec[i] で n1=0維持, 確率 1-ccp1_vec[i] で n1=1へ
    # n1=1のとき: 確率 1-ccp1_vec[i] で n1=0へ, 確率 ccp1_vec[i] で n1=1維持
    v = ccp1_vec
    rows1 = np.array([
        [v[0], 1 - v[0]],
        [v[1], 1 - v[1]],
        [1 - v[2], v[2]],
        [1 - v[3], v[3]],
        [v[4], 1 - v[4]],
        [v[5], 1 - v[5]],
        [1 - v[6], v[6]],
        [1 - v[7], v[7]],
    ])  # (8, 2)
    TempMat1 = np.kron(rows1, np.array([[1, 1]]))  # (8, 4)
    TempMat1 = np.hstack([TempMat1, TempMat1])     # (8, 8)

    # 企業2の行動による遷移
    v2 = ccp2_vec
    rows2 = np.array([
        [v2[0], 1 - v2[0]],
        [1 - v2[1], v2[1]],
        [v2[2], 1 - v2[2]],
        [1 - v2[3], v2[3]],
        [v2[4], 1 - v2[4]],
        [1 - v2[5], v2[5]],
        [v2[6], 1 - v2[6]],
        [1 - v2[7], v2[7]],
    ])  # (8, 2)
    TempMat2 = np.kron(np.array([[1, 1, 1, 1]]), rows2)  # (8, 8)

    output = TempMat0 * TempMat1 * TempMat2
    return output


def fP_a1given(trans_mat, ccp2_vec):
    \"\"\"企業1の行動が与えられたときの遷移行列を構築する。

    Returns
    -------
    list of 3 arrays (8, 8)
        [0]: a1 = -1 (退出) の時の遷移行列
        [1]: a1 = 0  (現状維持) の時の遷移行列
        [2]: a1 = +1 (参入) の時の遷移行列
    \"\"\"
    # 外生状態の遷移
    TempMat0 = np.kron(trans_mat, np.ones((4, 4)))

    # 企業2のCCPによる遷移
    v2 = ccp2_vec
    rows2 = np.array([
        [v2[0], 1 - v2[0], v2[0], 1 - v2[0]],
        [1 - v2[1], v2[1], 1 - v2[1], v2[1]],
        [v2[2], 1 - v2[2], v2[2], 1 - v2[2]],
        [1 - v2[3], v2[3], 1 - v2[3], v2[3]],
        [v2[4], 1 - v2[4], v2[4], 1 - v2[4]],
        [1 - v2[5], v2[5], 1 - v2[5], v2[5]],
        [v2[6], 1 - v2[6], v2[6], 1 - v2[6]],
        [1 - v2[7], v2[7], 1 - v2[7], v2[7]],
    ])  # (8, 4)
    TempMat2 = np.hstack([rows2, rows2])  # (8, 8)

    # a1 = -1 (退出): n1 は 1->0 に固定
    # R: MatAdjustMinus1 <- matrix(rep(rep(c(0,1), each=16), 2), nrow=8, byrow=TRUE)
    # vec = [0]*16 + [1]*16 + [0]*16 + [1]*16 (length 64)
    # matrix(vec, nrow=8, byrow=TRUE) => rows 0-1: all 0, rows 2-3: all 1, rows 4-5: all 0, rows 6-7: all 1
    MatAdjustMinus1 = np.array(([0]*16 + [1]*16 + [0]*16 + [1]*16), dtype=float).reshape(8, 8)

    # R: MatAdjustMinus2 <- matrix(rep(rep(c(1,0), each=16), 2), nrow=8)
    # vec = [1]*16 + [0]*16 + [1]*16 + [0]*16 (length 64)
    # matrix(vec, nrow=8) column-major => col 0-1: 1, col 2-3: 0, col 4-5: 1, col 6-7: 0
    MatAdjustMinus2 = np.array(([1]*16 + [0]*16 + [1]*16 + [0]*16), dtype=float).reshape(8, 8, order='F')

    output1 = TempMat0 * TempMat2 * MatAdjustMinus1 * MatAdjustMinus2

    # a1 = 0 (現状維持): n1 は変わらない
    ForZero = np.array([[1, 0], [0, 1]])
    MatAdjustZero = np.kron(ForZero, np.ones((2, 2)))  # (4, 4)
    MatAdjustZero = np.hstack([MatAdjustZero, MatAdjustZero])  # (4, 8)
    MatAdjustZero = np.vstack([MatAdjustZero, MatAdjustZero])  # (8, 8)
    output2 = TempMat0 * TempMat2 * MatAdjustZero

    # a1 = +1 (参入): n1 は 0->1 に固定
    # R: MatAdjustPlus1 <- matrix(rep(rep(c(1,0), each=16), 2), nrow=8, byrow=TRUE)
    # vec = [1]*16 + [0]*16 + [1]*16 + [0]*16 (length 64)
    # rows 0-1: all 1, rows 2-3: all 0, rows 4-5: all 1, rows 6-7: all 0
    MatAdjustPlus1 = np.array(([1]*16 + [0]*16 + [1]*16 + [0]*16), dtype=float).reshape(8, 8)

    # R: MatAdjustPlus2 <- matrix(rep(rep(c(0,1), each=16), 2), nrow=8)
    # vec = [0]*16 + [1]*16 + [0]*16 + [1]*16 (length 64)
    # column-major => col 0-1: 0, col 2-3: 1, col 4-5: 0, col 6-7: 1
    MatAdjustPlus2 = np.array(([0]*16 + [1]*16 + [0]*16 + [1]*16), dtype=float).reshape(8, 8, order='F')

    output3 = TempMat0 * TempMat2 * MatAdjustPlus1 * MatAdjustPlus2

    return [output1, output2, output3]


def fP_a2given(trans_mat, ccp1_vec):
    \"\"\"企業2の行動が与えられたときの遷移行列を構築する。

    Returns
    -------
    list of 3 arrays (8, 8)
        [0]: a2 = -1 (退出) の時の遷移行列
        [1]: a2 = 0  (現状維持) の時の遷移行列
        [2]: a2 = +1 (参入) の時の遷移行列
    \"\"\"
    # 外生状態の遷移
    TempMat0 = np.kron(trans_mat, np.ones((4, 4)))

    # 企業1のCCPによる遷移
    v = ccp1_vec
    rows1 = np.array([
        [v[0], v[0], 1 - v[0], 1 - v[0]],
        [v[1], v[1], 1 - v[1], 1 - v[1]],
        [1 - v[2], 1 - v[2], v[2], v[2]],
        [1 - v[3], 1 - v[3], v[3], v[3]],
        [v[4], v[4], 1 - v[4], 1 - v[4]],
        [v[5], v[5], 1 - v[5], 1 - v[5]],
        [1 - v[6], 1 - v[6], v[6], v[6]],
        [1 - v[7], 1 - v[7], v[7], v[7]],
    ])  # (8, 4)
    TempMat1 = np.kron(np.ones((1, 2)), rows1)  # (8, 8)

    # a2 = -1 (退出): n2 は 1->0 に固定
    # R: MatAdjustMinus1 <- matrix(rep(rep(c(0,1), each=8), 4), nrow=8, byrow=TRUE)
    # vec = ([0]*8 + [1]*8) * 4 (length 64)
    # matrix(vec, nrow=8, byrow=TRUE) => rows 0,2,4,6: all 0; rows 1,3,5,7: all 1
    MatAdjustMinus1 = np.array((([0]*8 + [1]*8) * 4), dtype=float).reshape(8, 8)

    # R: MatAdjustMinus2 <- matrix(rep(rep(c(1,0), each=8), 4), nrow=8)
    # vec = ([1]*8 + [0]*8) * 4 (length 64)
    # matrix(vec, nrow=8) column-major => col 0: 1, col 1: 0, col 2: 1, ...
    MatAdjustMinus2 = np.array((([1]*8 + [0]*8) * 4), dtype=float).reshape(8, 8, order='F')

    output1 = TempMat0 * TempMat1 * MatAdjustMinus1 * MatAdjustMinus2

    # a2 = 0 (現状維持): n2 は変わらない
    # R: ForZero <- c(1,0,1,0,1,0,1,0,0,1,0,1,0,1,0,1)
    # MatAdjustZero <- matrix(rep(ForZero, 4), ncol=8, byrow=TRUE)
    ForZero_vec = [1,0,1,0,1,0,1,0,0,1,0,1,0,1,0,1] * 4  # length 64
    MatAdjustZero = np.array(ForZero_vec, dtype=float).reshape(8, 8)
    output2 = TempMat0 * TempMat1 * MatAdjustZero

    # a2 = +1 (参入): n2 は 0->1 に固定
    # R: MatAdjustPlus1 <- matrix(rep(rep(c(1,0), each=8), 4), nrow=8, byrow=TRUE)
    # vec = ([1]*8 + [0]*8) * 4 (length 64)
    # rows 0,2,4,6: all 1; rows 1,3,5,7: all 0
    MatAdjustPlus1 = np.array((([1]*8 + [0]*8) * 4), dtype=float).reshape(8, 8)

    # R: MatAdjustPlus2 <- matrix(rep(rep(c(0,1), each=8), 4), nrow=8)
    # vec = ([0]*8 + [1]*8) * 4 (length 64)
    # column-major => col 0: 0, col 1: 1, col 2: 0, col 3: 1, ...
    MatAdjustPlus2 = np.array((([0]*8 + [1]*8) * 4), dtype=float).reshape(8, 8, order='F')

    output3 = TempMat0 * TempMat1 * MatAdjustPlus1 * MatAdjustPlus2

    return [output1, output2, output3]


# 動作確認: 遷移行列の行和がCCPと整合的かチェック
test_ccp1 = np.full(8, 0.5)
test_ccp2 = np.full(8, 0.5)
fPsigma_test = fP(TransitionMat, test_ccp1, test_ccp2)
print("fP(TransitionMat, 0.5, 0.5) の行和:")
print(fPsigma_test.sum(axis=1))
print("\\n全ての行和が1.0であることを確認。")
"""))

# ============================================================
# pi_Psigma functions
# ============================================================
cells.append(nbf.v4.new_markdown_cell("### 均衡利潤関数"))

cells.append(nbf.v4.new_code_cell("""\
def pi1PsigmaGen(pi1, CCP2Mat):
    \"\"\"企業1の均衡期待利潤を計算する (8 x 3)。

    各行動に対して、企業2のCCPで重み付けした期待利潤を計算。
    \"\"\"
    pi1_dec = pi1[:, 0]   # (8,)
    pi1_0   = pi1[:, 1]   # (8,)
    pi1_inc = pi1[:, 2]   # (8,)

    # (pi1_act * CCP2Mat) の行和 = 企業2の行動の期待値
    pi1_dec_p = (pi1_dec[:, np.newaxis] * CCP2Mat).sum(axis=1)  # (8,)
    pi1_0_p   = (pi1_0[:, np.newaxis]   * CCP2Mat).sum(axis=1)
    pi1_inc_p = (pi1_inc[:, np.newaxis] * CCP2Mat).sum(axis=1)

    return np.column_stack([pi1_dec_p, pi1_0_p, pi1_inc_p])


def pi2PsigmaGen(pi2, CCP1Mat):
    \"\"\"企業2の均衡期待利潤を計算する (8 x 3)。\"\"\"
    pi2_dec = pi2[:, 0]
    pi2_0   = pi2[:, 1]
    pi2_inc = pi2[:, 2]

    pi2_dec_p = (pi2_dec[:, np.newaxis] * CCP1Mat).sum(axis=1)
    pi2_0_p   = (pi2_0[:, np.newaxis]   * CCP1Mat).sum(axis=1)
    pi2_inc_p = (pi2_inc[:, np.newaxis] * CCP1Mat).sum(axis=1)

    return np.column_stack([pi2_dec_p, pi2_0_p, pi2_inc_p])


print("関数定義完了: pi1PsigmaGen, pi2PsigmaGen")
"""))

# ============================================================
# MPE Solver
# ============================================================
cells.append(nbf.v4.new_markdown_cell(
    "## 3. マルコフ完全均衡 (MPE) ソルバー\n\n"
    "以下のアルゴリズムで均衡CCPを計算する:\n\n"
    "1. CCPの初期値を設定\n"
    "2. 現在のCCPに基づいて事前価値関数を計算\n"
    "3. 事前価値関数に基づいてCCPを更新\n"
    "4. 更新されたCCPに基づいて事前価値関数を再計算\n"
    "5. 収束するまで3-4を繰り返す（収束基準: 価値関数の変化の二乗和 < 1e-12）"
))

cells.append(nbf.v4.new_code_cell("""\
def f_MPE(trans_mat, pi1, pi2, beta, ccp1_init=None, ccp2_init=None):
    \"\"\"マルコフ完全均衡のCCPを求める。

    Parameters
    ----------
    trans_mat : array (2, 2)
        景気の遷移行列
    pi1 : array (8, 3)
        企業1の利潤行列 (Adjuster適用済み)
    pi2 : array (8, 3)
        企業2の利潤行列 (Adjuster適用済み)
    beta : float
        割引因子
    ccp1_init : array (8,) or None
        企業1のCCP初期値。Noneの場合は全て0.5
    ccp2_init : array (8,) or None
        企業2のCCP初期値。Noneの場合は全て0.5

    Returns
    -------
    CCP1UpdatedMat : array (8, 3)
        均衡における企業1のCCP行列
    CCP2UpdatedMat : array (8, 3)
        均衡における企業2のCCP行列
    ExanteV1 : array (8,)
        均衡における企業1の事前価値関数
    ExanteV2 : array (8,)
        均衡における企業2の事前価値関数
    \"\"\"
    # Step 1: CCPの初期値
    CCP1 = np.full(8, 0.5) if ccp1_init is None else ccp1_init.copy()
    CCP2 = np.full(8, 0.5) if ccp2_init is None else ccp2_init.copy()

    CCP1Mat = CCP1Transform(CCP1)
    CCP2Mat = CCP2Transform(CCP2)

    # Step 2: 事前の価値関数を計算
    fPsigma = fP(trans_mat, CCP1, CCP2)

    pi1Psigma = pi1PsigmaGen(pi1, CCP2Mat)
    pi2Psigma = pi2PsigmaGen(pi2, CCP1Mat)

    eP1 = euler_gamma - CCP1LogTransform(CCP1)
    eP2 = euler_gamma - CCP2LogTransform(CCP2)

    # (4) 式: ExanteV = (I - beta * fPsigma)^{-1} * sum_a[ sigma(a|s) * (pi + eP) ]
    ones3 = np.ones(3)
    ExanteV1 = np.linalg.solve(
        np.eye(8) - beta * fPsigma,
        (CCP1Mat * (pi1Psigma + eP1)) @ ones3
    )
    ExanteV2 = np.linalg.solve(
        np.eye(8) - beta * fPsigma,
        (CCP2Mat * (pi2Psigma + eP2)) @ ones3
    )

    # Step 3: CCPを更新
    fP_a1 = fP_a1given(trans_mat, CCP2)
    fP_a2 = fP_a2given(trans_mat, CCP1)

    # 企業1
    NewSigmaSeed1 = (pi1Psigma + beta * np.column_stack([
        fP_a1[0] @ ExanteV1,
        fP_a1[1] @ ExanteV1,
        fP_a1[2] @ ExanteV1
    ])) * CCP1Adjuster
    NewSigmaDeno1 = np.exp(NewSigmaSeed1).sum(axis=1) - 1.0
    CCP1UpdatedMat = (np.exp(NewSigmaSeed1) / NewSigmaDeno1[:, np.newaxis]) * CCP1Adjuster
    CCP1Updated = CCP1UpdatedMat[:, 1]

    # 企業2
    NewSigmaSeed2 = (pi2Psigma + beta * np.column_stack([
        fP_a2[0] @ ExanteV2,
        fP_a2[1] @ ExanteV2,
        fP_a2[2] @ ExanteV2
    ])) * CCP2Adjuster
    NewSigmaDeno2 = np.exp(NewSigmaSeed2).sum(axis=1) - 1.0
    CCP2UpdatedMat = (np.exp(NewSigmaSeed2) / NewSigmaDeno2[:, np.newaxis]) * CCP2Adjuster
    CCP2Updated = CCP2UpdatedMat[:, 1]

    # Step 4: 更新されたCCPで事前価値関数を再計算
    fPsigma = fP(trans_mat, CCP1Updated, CCP2Updated)
    pi1Psigma = pi1PsigmaGen(pi1, CCP2UpdatedMat)
    pi2Psigma = pi2PsigmaGen(pi2, CCP1UpdatedMat)
    eP1 = euler_gamma - CCP1LogTransform(CCP1Updated)
    eP2 = euler_gamma - CCP2LogTransform(CCP2Updated)

    ExanteV1Updated = np.linalg.solve(
        np.eye(8) - beta * fPsigma,
        (CCP1UpdatedMat * (pi1Psigma + eP1)).sum(axis=1)
    )
    ExanteV2Updated = np.linalg.solve(
        np.eye(8) - beta * fPsigma,
        (CCP2UpdatedMat * (pi2Psigma + eP2)).sum(axis=1)
    )

    # Step 5: 収束するまで繰り返す
    DiffExanteV = np.sum((ExanteV1Updated - ExanteV1)**2 + (ExanteV2Updated - ExanteV2)**2)

    iteration = 0
    while DiffExanteV > 1.0e-12:
        iteration += 1
        CCP1 = CCP1Updated.copy()
        CCP2 = CCP2Updated.copy()
        ExanteV1 = ExanteV1Updated.copy()
        ExanteV2 = ExanteV2Updated.copy()

        # Step 3 再実行
        fP_a1 = fP_a1given(trans_mat, CCP2)
        fP_a2 = fP_a2given(trans_mat, CCP1)

        NewSigmaSeed1 = (pi1Psigma + beta * np.column_stack([
            fP_a1[0] @ ExanteV1,
            fP_a1[1] @ ExanteV1,
            fP_a1[2] @ ExanteV1
        ])) * CCP1Adjuster
        NewSigmaDeno1 = np.exp(NewSigmaSeed1).sum(axis=1) - 1.0
        CCP1UpdatedMat = (np.exp(NewSigmaSeed1) / NewSigmaDeno1[:, np.newaxis]) * CCP1Adjuster
        CCP1Updated = CCP1UpdatedMat[:, 1]

        NewSigmaSeed2 = (pi2Psigma + beta * np.column_stack([
            fP_a2[0] @ ExanteV2,
            fP_a2[1] @ ExanteV2,
            fP_a2[2] @ ExanteV2
        ])) * CCP2Adjuster
        NewSigmaDeno2 = np.exp(NewSigmaSeed2).sum(axis=1) - 1.0
        CCP2UpdatedMat = (np.exp(NewSigmaSeed2) / NewSigmaDeno2[:, np.newaxis]) * CCP2Adjuster
        CCP2Updated = CCP2UpdatedMat[:, 1]

        # Step 4 再実行
        fPsigma = fP(trans_mat, CCP1Updated, CCP2Updated)
        pi1Psigma = pi1PsigmaGen(pi1, CCP2UpdatedMat)
        pi2Psigma = pi2PsigmaGen(pi2, CCP1UpdatedMat)
        eP1 = euler_gamma - CCP1LogTransform(CCP1Updated)
        eP2 = euler_gamma - CCP2LogTransform(CCP2Updated)

        ExanteV1Updated = np.linalg.solve(
            np.eye(8) - beta * fPsigma,
            (CCP1UpdatedMat * (pi1Psigma + eP1)).sum(axis=1)
        )
        ExanteV2Updated = np.linalg.solve(
            np.eye(8) - beta * fPsigma,
            (CCP2UpdatedMat * (pi2Psigma + eP2)).sum(axis=1)
        )

        DiffExanteV = np.sum((ExanteV1Updated - ExanteV1)**2 + (ExanteV2Updated - ExanteV2)**2)

    print(f"  収束: {iteration} 回の反復")
    return CCP1UpdatedMat, CCP2UpdatedMat, ExanteV1Updated, ExanteV2Updated
"""))

# ============================================================
# Compute Equilibrium
# ============================================================
cells.append(nbf.v4.new_markdown_cell("## 4. 均衡の計算"))

cells.append(nbf.v4.new_code_cell("""\
# デフォルトの初期値 (全て0.5) で均衡を計算
print("均衡CCPの計算（初期値: 全て0.5）...")
CCP1Mat_eq, CCP2Mat_eq, ExanteV1_eq, ExanteV2_eq = f_MPE(
    TransitionMat, pi1, pi2, beta
)

print("\\n均衡CCP行列 [企業1 | 企業2]:")
state_labels = ['G,0,0', 'G,0,1', 'G,1,0', 'G,1,1',
                'B,0,0', 'B,0,1', 'B,1,0', 'B,1,1']
col_labels = ['a1=-1', 'a1=0', 'a1=+1', 'a2=-1', 'a2=0', 'a2=+1']

eq_mat = np.hstack([CCP1Mat_eq, CCP2Mat_eq])
df_eq = pd.DataFrame(eq_mat, index=state_labels, columns=col_labels)
print(df_eq.round(6).to_string())

print("\\n事前価値関数:")
df_v = pd.DataFrame({
    'ExanteV1': ExanteV1_eq,
    'ExanteV2': ExanteV2_eq
}, index=state_labels)
print(df_v.round(6).to_string())
"""))

# ============================================================
# Multiple Equilibria Search
# ============================================================
cells.append(nbf.v4.new_markdown_cell(
    "### 複数均衡の探索\n\n"
    "異なる初期値から均衡計算を行い、複数のMPEが存在するか確認する。"
))

cells.append(nbf.v4.new_code_cell("""\
# 複数の初期値から均衡を探索
np.random.seed(42)
num_trials = 20
equilibria = []
tol_eq = 1e-6  # 均衡の一致判定のトレランス

for trial in range(num_trials):
    if trial == 0:
        ccp1_init = np.full(8, 0.5)
        ccp2_init = np.full(8, 0.5)
    else:
        ccp1_init = np.random.uniform(0.1, 0.9, 8)
        ccp2_init = np.random.uniform(0.1, 0.9, 8)

    try:
        CCP1Mat_t, CCP2Mat_t, V1_t, V2_t = f_MPE(
            TransitionMat, pi1, pi2, beta,
            ccp1_init=ccp1_init, ccp2_init=ccp2_init
        )

        # 既知の均衡と比較
        is_new = True
        for eq in equilibria:
            diff = np.max(np.abs(eq[0] - CCP1Mat_t)) + np.max(np.abs(eq[1] - CCP2Mat_t))
            if diff < tol_eq:
                is_new = False
                break

        if is_new:
            equilibria.append((CCP1Mat_t, CCP2Mat_t, V1_t, V2_t))
            print(f"  新しい均衡を発見 (試行 {trial})")
    except Exception as e:
        print(f"  試行 {trial}: 収束しませんでした - {e}")

print(f"\\n発見された均衡の数: {len(equilibria)}")

# 最初の均衡を使用（デフォルト初期値からの結果）
CCP1UpdatedMat = equilibria[0][0]
CCP2UpdatedMat = equilibria[0][1]

print("\\n使用する均衡CCP (企業1の現状維持確率):")
print(CCP1UpdatedMat[:, 1].round(6))
print("\\n使用する均衡CCP (企業2の現状維持確率):")
print(CCP2UpdatedMat[:, 1].round(6))
"""))

# ============================================================
# Display Results
# ============================================================
cells.append(nbf.v4.new_markdown_cell("### 均衡CCPの表示と保存"))

cells.append(nbf.v4.new_code_cell("""\
# 均衡CCPの表を作成して表示
print("=" * 70)
print("均衡における条件付き選択確率 (CCP)")
print("=" * 70)

state_labels = ['G,0,0', 'G,0,1', 'G,1,0', 'G,1,1',
                'B,0,0', 'B,0,1', 'B,1,0', 'B,1,1']

# 企業1
print("\\n企業1:")
df_ccp1 = pd.DataFrame(
    CCP1UpdatedMat,
    index=state_labels,
    columns=['a1=-1 (退出)', 'a1=0 (現状維持)', 'a1=+1 (参入)']
)
print(df_ccp1.round(6).to_string())

# 企業2
print("\\n企業2:")
df_ccp2 = pd.DataFrame(
    CCP2UpdatedMat,
    index=state_labels,
    columns=['a2=-1 (退出)', 'a2=0 (現状維持)', 'a2=+1 (参入)']
)
print(df_ccp2.round(6).to_string())

# CSVに保存 (R版と同じフォーマット)
eq_mat_save = np.hstack([CCP1UpdatedMat, CCP2UpdatedMat])
df_save = pd.DataFrame(
    eq_mat_save,
    columns=['a1=-1', 'a1=0', 'a1=+1', 'a2=-1', 'a2=0', 'a2=+1']
)
df_save.index = range(1, 9)  # R と同じ 1-indexed
df_save.to_csv(output_dir / 'Sec_10_3_2_Equilibrium_CCP_python.csv')
print("\\n保存: Sec_10_3_2_Equilibrium_CCP_python.csv")
"""))

# ============================================================
# Compare with R results
# ============================================================
cells.append(nbf.v4.new_markdown_cell("### Rの結果との比較"))

cells.append(nbf.v4.new_code_cell("""\
# Rの結果と比較
r_csv_path = output_dir / 'Sec_10_3_2_Equilibrium_CCP.csv'
if r_csv_path.exists():
    df_r = pd.read_csv(r_csv_path, index_col=0)
    r_values = df_r.values
    py_values = np.hstack([CCP1UpdatedMat, CCP2UpdatedMat])
    max_diff = np.max(np.abs(r_values - py_values))
    print(f"R版との最大差: {max_diff:.2e}")
    if max_diff < 1e-8:
        print("結果はR版と一致しています。")
    else:
        print("注意: R版との差が大きい可能性があります。")
        print("\\nPython版:")
        print(np.round(py_values, 8))
        print("\\nR版:")
        print(np.round(r_values, 8))
else:
    print("R版の結果ファイルが見つかりません。比較をスキップします。")
"""))

# ============================================================
# Data Generation
# ============================================================
cells.append(nbf.v4.new_markdown_cell(
    "## 5. 疑似データの生成\n\n"
    "均衡CCPに基づいて、500市場 x 50期間の疑似データを生成する。\n"
    "各期間で、一様乱数とCCPを比較して企業の行動を決定する。"
))

cells.append(nbf.v4.new_code_cell("""\
# シミュレーション設定
np.random.seed(2023)  # R版と同じシード

NumSimMarkets = 500   # 市場の数
NumSimPeriods = 50    # 期間
NumSimFirms = 2       # 企業の数

# 各市場の初期状態変数 (1~8 の一様乱数, R の 1-indexed に合わせる)
InitialState = np.random.randint(1, 9, size=NumSimMarkets)  # 1~8

# 乱数の生成 (市場 x 期間 x (企業数+1))
RandomNumbers = np.random.uniform(0, 1, size=(NumSimMarkets, NumSimPeriods, NumSimFirms + 1))

# 疑似データの格納用行列
# 列: [MarketID, Period, State, Economy, n1, n2, a1, a2]
FakeData = np.zeros((NumSimMarkets * NumSimPeriods, 8))

# 状態変数 -> (景気, n1, n2) の対応表 (1-indexed の状態番号)
# state 1: (1, 0, 0) = Good, n1=0, n2=0
# state 2: (1, 0, 1) = Good, n1=0, n2=1
# state 3: (1, 1, 0) = Good, n1=1, n2=0
# state 4: (1, 1, 1) = Good, n1=1, n2=1
# state 5: (2, 0, 0) = Bad,  n1=0, n2=0
# state 6: (2, 0, 1) = Bad,  n1=0, n2=1
# state 7: (2, 1, 0) = Bad,  n1=1, n2=0
# state 8: (2, 1, 1) = Bad,  n1=1, n2=1
state_to_info = {
    1: (1, 0, 0), 2: (1, 0, 1), 3: (1, 1, 0), 4: (1, 1, 1),
    5: (2, 0, 0), 6: (2, 0, 1), 7: (2, 1, 0), 8: (2, 1, 1)
}

# 行動の決定ルール (状態 s に応じて、乱数 > CCP(s,1) の場合にどの行動を取るか)
# n1=0 の状態: a1 は 0 (維持) or +1 (参入)。乱数 > CCP1(s,col=1) なら a1=+1
# n1=1 の状態: a1 は 0 (維持) or -1 (退出)。乱数 > CCP1(s,col=1) なら a1=-1
# n2=0 の状態: a2 は 0 (維持) or +1 (参入)。乱数 > CCP2(s,col=1) なら a2=+1
# n2=1 の状態: a2 は 0 (維持) or -1 (退出)。乱数 > CCP2(s,col=1) なら a2=-1

# 状態番号(1-indexed) -> 企業1の「乱数がCCPを超えた場合」の行動
action_if_exceed_ccp1 = {
    1: 1, 2: 1, 3: -1, 4: -1,
    5: 1, 6: 1, 7: -1, 8: -1
}
action_if_exceed_ccp2 = {
    1: 1, 2: -1, 3: 1, 4: -1,
    5: 1, 6: -1, 7: 1, 8: -1
}

for m in range(NumSimMarkets):
    for t in range(NumSimPeriods):
        row = m * NumSimPeriods + t

        # 市場IDと期間
        FakeData[row, 0] = m + 1       # 1-indexed Market ID
        FakeData[row, 1] = t + 1       # 1-indexed Period

        if t == 0:
            # 初期状態
            s = InitialState[m]
            FakeData[row, 2] = s

            econ, n1_val, n2_val = state_to_info[s]
            FakeData[row, 3] = econ
        else:
            # 前期の情報
            prev_row = row - 1
            sprev = int(FakeData[prev_row, 2])
            a1prev = int(FakeData[prev_row, 6])
            a2prev = int(FakeData[prev_row, 7])

            # 景気の遷移
            if sprev >= 1 and sprev <= 4:
                # 前期は好景気
                if RandomNumbers[m, t, 2] < TransitionMat[0, 0]:
                    snow = 1  # 好景気維持
                else:
                    snow = 2  # 不景気へ
            else:
                # 前期は不景気
                if RandomNumbers[m, t, 2] < TransitionMat[1, 1]:
                    snow = 2  # 不景気維持
                else:
                    snow = 1  # 好景気へ

            FakeData[row, 3] = snow

            # 出店状況の更新
            n1t = int(FakeData[prev_row, 4]) + a1prev
            n2t = int(FakeData[prev_row, 5]) + a2prev

            # 状態変数の決定
            s = (snow - 1) * 4 + n1t * 2 + n2t + 1
            FakeData[row, 2] = s

        # 状態情報の設定
        s = int(FakeData[row, 2])
        econ, n1_val, n2_val = state_to_info[s]
        FakeData[row, 3:6] = [econ, n1_val, n2_val]

        # 企業1の行動: 0-indexed でアクセス (s-1 で Python の行インデックス)
        if RandomNumbers[m, t, 0] > CCP1UpdatedMat[s - 1, 1]:
            FakeData[row, 6] = action_if_exceed_ccp1[s]
        # else: FakeData[row, 6] = 0 (既にゼロ初期化済み)

        # 企業2の行動
        if RandomNumbers[m, t, 1] > CCP2UpdatedMat[s - 1, 1]:
            FakeData[row, 7] = action_if_exceed_ccp2[s]

print(f"疑似データ生成完了: {FakeData.shape[0]} 行 x {FakeData.shape[1]} 列")
print(f"  市場数: {NumSimMarkets}")
print(f"  期間数: {NumSimPeriods}")
"""))

# ============================================================
# Display Generated Data
# ============================================================
cells.append(nbf.v4.new_markdown_cell("### 生成されたデータの確認"))

cells.append(nbf.v4.new_code_cell("""\
# データフレームに変換
col_names = ['MarketID', 'Period', 'State', 'Economy', 'n1', 'n2', 'a1', 'a2']
df_fake = pd.DataFrame(FakeData, columns=col_names)
df_fake = df_fake.astype({
    'MarketID': int, 'Period': int, 'State': int,
    'Economy': int, 'n1': int, 'n2': int, 'a1': int, 'a2': int
})

print("先頭20行:")
print(df_fake.head(20).to_string(index=False))

print("\\n記述統計:")
print(df_fake.describe().round(3).to_string())

print(f"\\n状態変数の分布:")
print(df_fake['State'].value_counts().sort_index())

print(f"\\n企業1の行動の分布:")
print(df_fake['a1'].value_counts().sort_index())

print(f"\\n企業2の行動の分布:")
print(df_fake['a2'].value_counts().sort_index())
"""))

# ============================================================
# Save Data
# ============================================================
cells.append(nbf.v4.new_markdown_cell("### データの保存"))

cells.append(nbf.v4.new_code_cell("""\
# CSVに保存 (R版と同じフォーマット: ヘッダーなし)
save_path = output_dir / 'FakeData_Python.csv'
df_fake.to_csv(save_path, index=False, header=False)
print(f"保存: {save_path}")

# R版のデータとの比較
r_data_path = output_dir / 'FakeData_R.csv'
if r_data_path.exists():
    df_r_data = pd.read_csv(r_data_path, index_col=0, header=0)
    df_r_data.columns = col_names

    # CCPが一致しているかを確認（ランダムシードが異なるので完全一致はしない）
    print(f"\\nR版のデータ形状: {df_r_data.shape}")
    print(f"Python版のデータ形状: {df_fake.shape}")

    # 状態分布の比較
    print("\\n状態分布の比較:")
    r_state_dist = df_r_data['State'].value_counts().sort_index()
    py_state_dist = df_fake['State'].value_counts().sort_index()
    comparison = pd.DataFrame({
        'R': r_state_dist,
        'Python': py_state_dist
    })
    print(comparison.to_string())
else:
    print("R版のデータファイルが見つかりません。比較をスキップします。")
"""))

# ============================================================
# Compare with Matlab Data
# ============================================================
cells.append(nbf.v4.new_markdown_cell("### Matlabデータとの比較"))

cells.append(nbf.v4.new_code_cell("""\
# Matlabで生成されたデータの読み込み
matlab_path = data_dir / 'FakeData_Matlab.csv'
if matlab_path.exists():
    df_matlab = pd.read_csv(matlab_path, header=None,
                            names=col_names)
    print(f"Matlabデータ形状: {df_matlab.shape}")
    print("\\n先頭10行:")
    print(df_matlab.head(10).to_string(index=False))

    print("\\nMatlabデータの状態分布:")
    print(df_matlab['State'].value_counts().sort_index())

    print("\\nMatlabデータの行動分布 (企業1):")
    print(df_matlab['a1'].value_counts().sort_index())
else:
    print("Matlabデータが見つかりません。")
"""))

# ============================================================
# Summary
# ============================================================
cells.append(nbf.v4.new_markdown_cell("## まとめ"))

cells.append(nbf.v4.new_code_cell("""\
print("=" * 60)
print("第10章の分析完了")
print("=" * 60)
print("\\n実行内容:")
print("  1. モデルの基本パラメータを設定")
print("  2. 利潤関数と遷移行列を構築")
print("  3. マルコフ完全均衡のCCPを計算")
print("  4. 均衡CCPに基づいて疑似データを生成")
print(f"  5. 発見された均衡の数: {len(equilibria)}")
print(f"  6. 疑似データ: {NumSimMarkets} 市場 x {NumSimPeriods} 期間 = {FakeData.shape[0]} 観測")
print("\\n出力ファイル:")
for f in sorted(output_dir.glob('*')):
    if 'python' in f.name.lower() or 'Python' in f.name:
        print(f"  {f.name}")
"""))

nb.cells = cells
nbf.write(nb, 'main_ch10_equilibrium.ipynb')
print("Generated: main_ch10_equilibrium.ipynb")
