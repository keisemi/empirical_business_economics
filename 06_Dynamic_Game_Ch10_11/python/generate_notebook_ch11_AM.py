"""Generate main_ch11_AM.ipynb for Chapter 11: AM (Aguirregabiria-Mira 2007) Estimation."""
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
    "# 第11章 動的ゲームの推定：Aguirregabiria-Mira (2007) 推定量\n"
    "\n"
    "2企業の参入退出の動的ゲームにおいて、\n"
    "Aguirregabiria and Mira (2007) のCCP反復推定量を用いてパラメータを推定する。\n"
    "\n"
    "手順:\n"
    "1. 均衡CCPの計算と疑似データの生成 (sub_1_prepare + sub_2_DGP 相当)\n"
    "2. Matlab生成のFakeDataを読み込み\n"
    "3. データからCCPと遷移確率を推定\n"
    "4. AM推定量によるパラメータ推定 (CCP反復)\n"
    "5. Bootstrapによる標準誤差の計算\n"
    "6. 結果の表示と保存"
))

# ============================================================
# Cell 2: Setup
# ============================================================
cells.append(nbf.v4.new_code_cell("""\
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from pathlib import Path
import time
import warnings
warnings.filterwarnings('ignore')

# パス設定
base_dir = Path('..')
output_dir = base_dir / 'output'
data_dir = base_dir / 'data_from_matlab'
output_dir.mkdir(exist_ok=True)

print("Setup complete.")
"""))

# ============================================================
# Cell 3: Parameters (sub_1_prepare)
# ============================================================
cells.append(nbf.v4.new_markdown_cell(
    "## パラメータの設定 (sub_1_prepare 相当)\n"
    "\n"
    "割引因子、遷移行列、利潤パラメータ等の基本設定を行う。"
))

cells.append(nbf.v4.new_code_cell("""\
# 割引因子
beta = 0.8

# オイラー定数
eulergamma = 0.5772

# 景気の遷移行列: [[P(G|G), P(B|G)], [P(G|B), P(B|B)]]
TransitionMat = np.array([[0.7, 0.3],
                          [0.4, 0.6]])

# パラメータの設定
# Parameters[0]: 企業1のベース利潤
# Parameters[1]: 企業2のベース利潤
# Parameters[2]: 顧客収奪効果
# Parameters[3]: 景気が良い時の追加的利潤
# Parameters[4]: 退出のためのコスト
# Parameters[5]: 参入のためのコスト
Parameters = np.array([0.3, 0.2, -0.27, 0.45, -0.15, -2.10])

# パラメータを一般化しやすいよう並べ替える (10x1)
# theta[0]: 企業1のベース利潤
# theta[1]: ライバルの店舗数が企業1の利潤に与える影響
# theta[2]: 景気が良い時の企業1への追加的利潤
# theta[3]: 企業1の退出のためのコスト
# theta[4]: 企業1の出店のためのコスト
# theta[5]-theta[9]: 企業2について同様
TrueParameterValues = np.array([
    Parameters[0],  # 企業1のベース利潤
    Parameters[2],  # ライバルの店舗数が企業1の利潤に与える影響
    Parameters[3],  # 景気が良い時の企業1への追加的利潤
    Parameters[4],  # 企業1の退出のためのコスト
    Parameters[5],  # 企業1の出店のためのコスト
    Parameters[1],  # 企業2のベース利潤
    Parameters[2],  # ライバルの店舗数が企業2の利潤に与える影響
    Parameters[3],  # 景気が良い時の企業2への追加的利潤
    Parameters[4],  # 企業2の退出のためのコスト
    Parameters[5],  # 企業2の出店のためのコスト
])

# CCP Adjuster: 各状態で選択可能な行動に1、不可能な行動に0を割り当てる行列 (8x3)
# 列は a_i = -1, 0, 1 に対応
# 状態: G00, G01, G10, G11, B00, B01, B10, B11
# 企業1について:
#   n1=0 のとき: a1=-1 は不可, a1=0 は可, a1=1 は可 -> [0,1,1]
#   n1=1 のとき: a1=-1 は可, a1=0 は可, a1=1 は不可 -> [1,1,0]
CCP1Adjuster = np.array([
    [0, 1, 1],  # G00: n1=0
    [0, 1, 1],  # G01: n1=0
    [1, 1, 0],  # G10: n1=1
    [1, 1, 0],  # G11: n1=1
    [0, 1, 1],  # B00: n1=0
    [0, 1, 1],  # B01: n1=0
    [1, 1, 0],  # B10: n1=1
    [1, 1, 0],  # B11: n1=1
], dtype=float)

# 企業2について:
#   n2=0 のとき: [0,1,1]
#   n2=1 のとき: [1,1,0]
CCP2Adjuster = np.array([
    [0, 1, 1],  # G00: n2=0
    [1, 1, 0],  # G01: n2=1
    [0, 1, 1],  # G10: n2=0
    [1, 1, 0],  # G11: n2=1
    [0, 1, 1],  # B00: n2=0
    [1, 1, 0],  # B01: n2=1
    [0, 1, 1],  # B10: n2=0
    [1, 1, 0],  # B11: n2=1
], dtype=float)

print("パラメータ設定完了")
print(f"beta = {beta}")
print(f"遷移行列:\\n{TransitionMat}")
print(f"Parameters = {Parameters}")
print(f"TrueParameterValues = {TrueParameterValues}")
"""))

# ============================================================
# Cell 4: Helper Functions - Profit
# ============================================================
cells.append(nbf.v4.new_markdown_cell("## ヘルパー関数の定義: 利潤関数"))

cells.append(nbf.v4.new_code_cell("""\
def pi1gen(theta):
    \"\"\"企業1の利潤行列を生成する (8x3)
    行: 状態 (G00, G01, G10, G11, B00, B01, B10, B11)
    列: 行動 (a1=-1, a1=0, a1=1)
    \"\"\"
    # base: 各状態での基本利潤 (8x1)
    # theta[0]: 企業1のベース利潤, theta[1]: ライバル効果, theta[2]: 好景気効果
    base = np.array([
        0,                              # G00: n1=0
        0,                              # G01: n1=0
        theta[0] + theta[2],            # G10: n1=1, 好景気
        theta[0] + theta[1] + theta[2], # G11: n1=1, 好景気, ライバルあり
        0,                              # B00: n1=0
        0,                              # B01: n1=0
        theta[0],                       # B10: n1=1
        theta[0] + theta[1],            # B11: n1=1, ライバルあり
    ])

    # invdiv: 参入退出コスト (1x3)
    # theta[3]: 退出コスト, theta[4]: 参入コスト
    invdiv = np.array([theta[3], 0, theta[4]])

    # 8x3 行列: base を3列に拡張 + invdiv を8行に拡張
    output = np.tile(base.reshape(-1, 1), (1, 3)) + np.tile(invdiv.reshape(1, -1), (8, 1))
    return output


def pi2gen(theta):
    \"\"\"企業2の利潤行列を生成する (8x3)\"\"\"
    # theta[5]: 企業2のベース利潤, theta[6]: ライバル効果, theta[7]: 好景気効果
    base = np.array([
        0,                              # G00: n2=0
        theta[5] + theta[7],            # G01: n2=1, 好景気
        0,                              # G10: n2=0
        theta[5] + theta[6] + theta[7], # G11: n2=1, 好景気, ライバルあり
        0,                              # B00: n2=0
        theta[5],                       # B01: n2=1
        0,                              # B10: n2=0
        theta[5] + theta[6],            # B11: n2=1, ライバルあり
    ])

    # theta[8]: 退出コスト, theta[9]: 参入コスト
    invdiv = np.array([theta[8], 0, theta[9]])

    output = np.tile(base.reshape(-1, 1), (1, 3)) + np.tile(invdiv.reshape(1, -1), (8, 1))
    return output


# 真のパラメータで利潤行列を計算
pi1 = pi1gen(TrueParameterValues) * CCP1Adjuster
pi2 = pi2gen(TrueParameterValues) * CCP2Adjuster

print("企業1の利潤行列 pi1 (8x3):")
print(pi1)
print("\\n企業2の利潤行列 pi2 (8x3):")
print(pi2)
"""))

# ============================================================
# Cell 5: Helper Functions - CCP Transform
# ============================================================
cells.append(nbf.v4.new_markdown_cell("## ヘルパー関数の定義: CCP変換"))

cells.append(nbf.v4.new_code_cell("""\
def CCP1Transform(x):
    \"\"\"CCP1ベクトル (8,) を CCP1行列 (8x3) に変換する
    x[s]: 状態sにおける現状維持(a1=0)の確率
    \"\"\"
    return np.array([
        [0,        x[0], 1-x[0]],  # G00: n1=0
        [0,        x[1], 1-x[1]],  # G01: n1=0
        [1-x[2],   x[2], 0      ],  # G10: n1=1
        [1-x[3],   x[3], 0      ],  # G11: n1=1
        [0,        x[4], 1-x[4]],  # B00: n1=0
        [0,        x[5], 1-x[5]],  # B01: n1=0
        [1-x[6],   x[6], 0      ],  # B10: n1=1
        [1-x[7],   x[7], 0      ],  # B11: n1=1
    ])


def CCP2Transform(x):
    \"\"\"CCP2ベクトル (8,) を CCP2行列 (8x3) に変換する\"\"\"
    return np.array([
        [0,        x[0], 1-x[0]],  # G00: n2=0
        [1-x[1],   x[1], 0      ],  # G01: n2=1
        [0,        x[2], 1-x[2]],  # G10: n2=0
        [1-x[3],   x[3], 0      ],  # G11: n2=1
        [0,        x[4], 1-x[4]],  # B00: n2=0
        [1-x[5],   x[5], 0      ],  # B01: n2=1
        [0,        x[6], 1-x[6]],  # B10: n2=0
        [1-x[7],   x[7], 0      ],  # B11: n2=1
    ])


def CCP1LogTransform(x):
    \"\"\"CCP1ベクトル (8,) を log CCP1行列 (8x3) に変換する\"\"\"
    return np.array([
        [0,              np.log(x[0]), np.log(1-x[0])],
        [0,              np.log(x[1]), np.log(1-x[1])],
        [np.log(1-x[2]), np.log(x[2]), 0             ],
        [np.log(1-x[3]), np.log(x[3]), 0             ],
        [0,              np.log(x[4]), np.log(1-x[4])],
        [0,              np.log(x[5]), np.log(1-x[5])],
        [np.log(1-x[6]), np.log(x[6]), 0             ],
        [np.log(1-x[7]), np.log(x[7]), 0             ],
    ])


def CCP2LogTransform(x):
    \"\"\"CCP2ベクトル (8,) を log CCP2行列 (8x3) に変換する\"\"\"
    return np.array([
        [0,              np.log(x[0]), np.log(1-x[0])],
        [np.log(1-x[1]), np.log(x[1]), 0             ],
        [0,              np.log(x[2]), np.log(1-x[2])],
        [np.log(1-x[3]), np.log(x[3]), 0             ],
        [0,              np.log(x[4]), np.log(1-x[4])],
        [np.log(1-x[5]), np.log(x[5]), 0             ],
        [0,              np.log(x[6]), np.log(1-x[6])],
        [np.log(1-x[7]), np.log(x[7]), 0             ],
    ])


print("CCP変換関数定義完了")
"""))

# ============================================================
# Cell 6: Helper Functions - Transition Matrices
# ============================================================
cells.append(nbf.v4.new_markdown_cell("## ヘルパー関数の定義: 遷移行列"))

cells.append(nbf.v4.new_code_cell("""\
def fP(Matrix1, Vec1, Vec2):
    \"\"\"CCPの下での状態遷移行列 F^{P,sigma} を構築する (8x8)

    Parameters
    ----------
    Matrix1 : (2,2) 景気の遷移行列
    Vec1 : (8,) 企業1のCCPベクトル (P(stay))
    Vec2 : (8,) 企業2のCCPベクトル (P(stay))
    \"\"\"
    # 景気の遷移確率を8x8に拡張
    TempMat0 = np.kron(Matrix1, np.ones((4, 4)))

    # 企業1の行動による状態遷移
    rows1 = np.array([
        [Vec1[0], 1-Vec1[0]],
        [Vec1[1], 1-Vec1[1]],
        [1-Vec1[2], Vec1[2]],
        [1-Vec1[3], Vec1[3]],
        [Vec1[4], 1-Vec1[4]],
        [Vec1[5], 1-Vec1[5]],
        [1-Vec1[6], Vec1[6]],
        [1-Vec1[7], Vec1[7]],
    ])
    TempMat1 = np.kron(rows1, np.ones((1, 2)))
    TempMat1 = np.hstack([TempMat1, TempMat1])

    # 企業2の行動による状態遷移
    rows2 = np.array([
        [Vec2[0], 1-Vec2[0]],
        [1-Vec2[1], Vec2[1]],
        [Vec2[2], 1-Vec2[2]],
        [1-Vec2[3], Vec2[3]],
        [Vec2[4], 1-Vec2[4]],
        [1-Vec2[5], Vec2[5]],
        [Vec2[6], 1-Vec2[6]],
        [1-Vec2[7], Vec2[7]],
    ])
    TempMat2 = np.kron(np.ones((1, 4)), rows2)

    output = TempMat0 * TempMat1 * TempMat2
    return output


def fP_a1given(Matrix1, Vec2):
    \"\"\"企業1の行動を条件付けた状態遷移行列のリスト [a1=-1, a1=0, a1=1] (各8x8)\"\"\"
    TempMat0 = np.kron(Matrix1, np.ones((4, 4)))

    rows2 = np.array([
        [Vec2[0], 1-Vec2[0], Vec2[0], 1-Vec2[0]],
        [1-Vec2[1], Vec2[1], 1-Vec2[1], Vec2[1]],
        [Vec2[2], 1-Vec2[2], Vec2[2], 1-Vec2[2]],
        [1-Vec2[3], Vec2[3], 1-Vec2[3], Vec2[3]],
        [Vec2[4], 1-Vec2[4], Vec2[4], 1-Vec2[4]],
        [1-Vec2[5], Vec2[5], 1-Vec2[5], Vec2[5]],
        [Vec2[6], 1-Vec2[6], Vec2[6], 1-Vec2[6]],
        [1-Vec2[7], Vec2[7], 1-Vec2[7], Vec2[7]],
    ])
    TempMat2 = np.hstack([rows2, rows2])

    # a1 = -1
    vec_m1 = np.concatenate([np.zeros(16), np.ones(16), np.zeros(16), np.ones(16)])
    MatAdjustMinus1 = vec_m1.reshape(8, 8, order='C')
    vec_m2 = np.concatenate([np.ones(16), np.zeros(16), np.ones(16), np.zeros(16)])
    MatAdjustMinus2 = vec_m2.reshape(8, 8, order='F')
    output1 = TempMat0 * TempMat2 * MatAdjustMinus1 * MatAdjustMinus2

    # a1 = 0
    ForZero = np.array([[1, 0], [0, 1]])
    block = np.kron(ForZero, np.ones((2, 2)))
    MatAdjustZero = np.hstack([block, block])
    MatAdjustZero = np.vstack([MatAdjustZero, MatAdjustZero])
    output2 = TempMat0 * TempMat2 * MatAdjustZero

    # a1 = 1
    vec_p1 = np.concatenate([np.ones(16), np.zeros(16), np.ones(16), np.zeros(16)])
    MatAdjustPlus1 = vec_p1.reshape(8, 8, order='C')
    vec_p2 = np.concatenate([np.zeros(16), np.ones(16), np.zeros(16), np.ones(16)])
    MatAdjustPlus2 = vec_p2.reshape(8, 8, order='F')
    output3 = TempMat0 * TempMat2 * MatAdjustPlus1 * MatAdjustPlus2

    return [output1, output2, output3]


def fP_a2given(Matrix1, Vec1):
    \"\"\"企業2の行動を条件付けた状態遷移行列のリスト [a2=-1, a2=0, a2=1] (各8x8)\"\"\"
    TempMat0 = np.kron(Matrix1, np.ones((4, 4)))

    rows1 = np.array([
        [Vec1[0], Vec1[0], 1-Vec1[0], 1-Vec1[0]],
        [Vec1[1], Vec1[1], 1-Vec1[1], 1-Vec1[1]],
        [1-Vec1[2], 1-Vec1[2], Vec1[2], Vec1[2]],
        [1-Vec1[3], 1-Vec1[3], Vec1[3], Vec1[3]],
        [Vec1[4], Vec1[4], 1-Vec1[4], 1-Vec1[4]],
        [Vec1[5], Vec1[5], 1-Vec1[5], 1-Vec1[5]],
        [1-Vec1[6], 1-Vec1[6], Vec1[6], Vec1[6]],
        [1-Vec1[7], 1-Vec1[7], Vec1[7], Vec1[7]],
    ])
    TempMat1 = np.kron(np.ones((1, 2)), rows1)

    # a2 = -1
    vec_m1 = np.tile(np.concatenate([np.zeros(8), np.ones(8)]), 4)
    MatAdjustMinus1 = vec_m1.reshape(8, 8, order='C')
    vec_m2 = np.tile(np.concatenate([np.ones(8), np.zeros(8)]), 4)
    MatAdjustMinus2 = vec_m2.reshape(8, 8, order='F')
    output1 = TempMat0 * TempMat1 * MatAdjustMinus1 * MatAdjustMinus2

    # a2 = 0
    ForZero = np.array([1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1], dtype=float)
    vec_z = np.tile(ForZero, 4)
    MatAdjustZero = vec_z.reshape(8, 8, order='C')
    output2 = TempMat0 * TempMat1 * MatAdjustZero

    # a2 = 1
    vec_p1 = np.tile(np.concatenate([np.ones(8), np.zeros(8)]), 4)
    MatAdjustPlus1 = vec_p1.reshape(8, 8, order='C')
    vec_p2 = np.tile(np.concatenate([np.zeros(8), np.ones(8)]), 4)
    MatAdjustPlus2 = vec_p2.reshape(8, 8, order='F')
    output3 = TempMat0 * TempMat1 * MatAdjustPlus1 * MatAdjustPlus2

    return [output1, output2, output3]


print("遷移行列関数定義完了")
"""))

# ============================================================
# Cell 7: Helper Functions - piPsigmaGen
# ============================================================
cells.append(nbf.v4.new_markdown_cell("## ヘルパー関数の定義: 期待利潤・収束判定"))

cells.append(nbf.v4.new_code_cell("""\
def pi1PsigmaGen(pi1, Mat2):
    \"\"\"企業2のCCP行列 Mat2 の下での企業1の期待利潤 (8x3)\"\"\"
    ones3 = np.ones(3)
    pi1_dec = (pi1[:, 0].reshape(-1, 1) * Mat2) @ ones3
    pi1_0   = (pi1[:, 1].reshape(-1, 1) * Mat2) @ ones3
    pi1_inc = (pi1[:, 2].reshape(-1, 1) * Mat2) @ ones3
    return np.column_stack([pi1_dec, pi1_0, pi1_inc])


def pi2PsigmaGen(pi2, Mat1):
    \"\"\"企業1のCCP行列 Mat1 の下での企業2の期待利潤 (8x3)\"\"\"
    ones3 = np.ones(3)
    pi2_dec = (pi2[:, 0].reshape(-1, 1) * Mat1) @ ones3
    pi2_0   = (pi2[:, 1].reshape(-1, 1) * Mat1) @ ones3
    pi2_inc = (pi2[:, 2].reshape(-1, 1) * Mat1) @ ones3
    return np.column_stack([pi2_dec, pi2_0, pi2_inc])


def check_convergence(oldCCP, newCCP, tol=1e-12):
    \"\"\"CCPが収束しているかどうか確認する
    oldCCP, newCCP: (8, 2) の行列 [ccp1, ccp2]
    Returns True if converged
    \"\"\"
    diff1 = np.dot(oldCCP[:, 0] - newCCP[:, 0], oldCCP[:, 0] - newCCP[:, 0])
    diff2 = np.dot(oldCCP[:, 1] - newCCP[:, 1], oldCCP[:, 1] - newCCP[:, 1])
    return (diff1 < tol) and (diff2 < tol)


print("期待利潤関数・収束判定関数定義完了")
"""))

# ============================================================
# Cell 8: CCP_to_Value_to_prediction
# ============================================================
cells.append(nbf.v4.new_markdown_cell(
    "## CCP_to_Value_to_prediction 関数\n"
    "\n"
    "パラメータとCCPから、新しいCCPを予測する関数。\n"
    "AM推定の中核をなす「Value inversion + CCP更新」を1ステップ実行する。"
))

cells.append(nbf.v4.new_code_cell("""\
def CCP_to_Value_to_prediction(theta, CCP1, CCP2, TransitionMat, beta):
    \"\"\"CCPから事前の価値をinvertし、その価値からCCPを予測する

    Parameters
    ----------
    theta : (10,) パラメータベクトル
    CCP1, CCP2 : (8,) CCP本質ベクトル
    TransitionMat : (2,2) 外生状態遷移行列
    beta : float 割引因子

    Returns
    -------
    output : (8, 2) 更新されたCCP [CCP1Updated, CCP2Updated]
    \"\"\"
    # Step 1: setup
    CCP1Adj = np.array([
        [0, 1, 1], [0, 1, 1], [1, 1, 0], [1, 1, 0],
        [0, 1, 1], [0, 1, 1], [1, 1, 0], [1, 1, 0],
    ], dtype=float)
    CCP2Adj = np.array([
        [0, 1, 1], [1, 1, 0], [0, 1, 1], [1, 1, 0],
        [0, 1, 1], [1, 1, 0], [0, 1, 1], [1, 1, 0],
    ], dtype=float)

    # 利潤行列の計算
    pi1_local = pi1gen(theta) * CCP1Adj
    pi2_local = pi2gen(theta) * CCP2Adj

    # CCPベクトルを行列に変換
    CCP1Mat = CCP1Transform(CCP1)
    CCP2Mat = CCP2Transform(CCP2)

    # Step 2: 事前の価値関数を計算 (Value inversion)
    fPsigma = fP(TransitionMat, CCP1, CCP2)
    pi1Psigma = pi1PsigmaGen(pi1_local, CCP2Mat)
    pi2Psigma = pi2PsigmaGen(pi2_local, CCP1Mat)
    eP1 = eulergamma - CCP1LogTransform(CCP1)
    eP2 = eulergamma - CCP2LogTransform(CCP2)

    inv_mat = np.linalg.inv(np.eye(8) - beta * fPsigma)
    ExanteV1 = inv_mat @ np.sum(CCP1Mat * (pi1Psigma + eP1), axis=1)
    ExanteV2 = inv_mat @ np.sum(CCP2Mat * (pi2Psigma + eP2), axis=1)

    # Step 3: 事前の価値関数からCCPを更新
    fP_a1_list = fP_a1given(TransitionMat, CCP2)
    fP_a2_list = fP_a2given(TransitionMat, CCP1)

    # 企業1
    future1 = np.column_stack([fP_a1_list[k] @ ExanteV1 for k in range(3)])
    NewSigmaSeed1 = (pi1Psigma + beta * future1) * CCP1Adj
    NewSigmaDeno1 = np.sum(np.exp(NewSigmaSeed1), axis=1) - 1.0
    NewSigma1 = np.exp(NewSigmaSeed1) / NewSigmaDeno1.reshape(-1, 1)
    CCP1UpdatedMat = NewSigma1 * CCP1Adj
    CCP1Updated = CCP1UpdatedMat[:, 1]

    # 企業2
    future2 = np.column_stack([fP_a2_list[k] @ ExanteV2 for k in range(3)])
    NewSigmaSeed2 = (pi2Psigma + beta * future2) * CCP2Adj
    NewSigmaDeno2 = np.sum(np.exp(NewSigmaSeed2), axis=1) - 1.0
    NewSigma2 = np.exp(NewSigmaSeed2) / NewSigmaDeno2.reshape(-1, 1)
    CCP2UpdatedMat = NewSigma2 * CCP2Adj
    CCP2Updated = CCP2UpdatedMat[:, 1]

    return np.column_stack([CCP1Updated, CCP2Updated])


print("CCP_to_Value_to_prediction 関数定義完了")
"""))

# ============================================================
# Cell 9: obj_lik function
# ============================================================
cells.append(nbf.v4.new_markdown_cell(
    "## 擬似尤度関数 obj_lik\n"
    "\n"
    "AM推定で用いる擬似対数尤度関数。\n"
    "パラメータとCCPから予測CCPを計算し、データとの尤度を求める。"
))

cells.append(nbf.v4.new_code_cell("""\
def obj_lik(param, GivenCCP1, GivenCCP2, TransitionProb, beta, Data):
    \"\"\"擬似対数尤度関数

    Parameters
    ----------
    param : (10,) パラメータベクトル
    GivenCCP1, GivenCCP2 : (8,) 与えられたCCP
    TransitionProb : (2,2) 遷移確率行列
    beta : float 割引因子
    Data : (N, 3) データ [action1, action2, state(1-indexed)]

    Returns
    -------
    ll : float 対数尤度 (最大化する)
    \"\"\"
    # パラメータとCCPから新たなCCPを計算
    CCPs = CCP_to_Value_to_prediction(param, GivenCCP1, GivenCCP2, TransitionProb, beta)
    # CCPs: (8, 2), 列0=CCP1, 列1=CCP2
    # CCPs[s, i] = P(a_i=0 | state s+1)  (0-indexed)

    ll = 0.0
    for s in range(8):
        s_1indexed = s + 1  # R は 1-indexed
        mask_s = (Data[:, 2] == s_1indexed)

        # 企業1: action==0 のとき log(CCP), action!=0 のとき log(1-CCP)
        mask_a1_zero = (Data[:, 0] == 0) & mask_s
        mask_a1_nonzero = (Data[:, 0] != 0) & mask_s
        ll += np.sum(mask_a1_zero) * np.log(CCPs[s, 0])
        ll += np.sum(mask_a1_nonzero) * np.log(1 - CCPs[s, 0])

        # 企業2: action==0 のとき log(CCP), action!=0 のとき log(1-CCP)
        mask_a2_zero = (Data[:, 1] == 0) & mask_s
        mask_a2_nonzero = (Data[:, 1] != 0) & mask_s
        ll += np.sum(mask_a2_zero) * np.log(CCPs[s, 1])
        ll += np.sum(mask_a2_nonzero) * np.log(1 - CCPs[s, 1])

    return ll


print("obj_lik 関数定義完了")
"""))

# ============================================================
# Cell 10: MPE Solver (for computing true equilibrium CCP)
# ============================================================
cells.append(nbf.v4.new_markdown_cell(
    "## MPEソルバー\n"
    "\n"
    "真のパラメータの下で均衡CCPを計算し、疑似データを生成する (sub_2_DGP 相当)。"
))

cells.append(nbf.v4.new_code_cell("""\
def f_MPE(TransitionMat, pi1, pi2, beta, eulergamma, CCP1Adjuster, CCP2Adjuster, tol=1e-12):
    \"\"\"マルコフ完全均衡（MPE）を固定点反復で計算する\"\"\"
    CCP1 = np.full(8, 0.5)
    CCP2 = np.full(8, 0.5)

    CCP1Mat = CCP1Transform(CCP1)
    CCP2Mat = CCP2Transform(CCP2)

    # Step 2
    fPsigma = fP(TransitionMat, CCP1, CCP2)
    pi1Psigma = pi1PsigmaGen(pi1, CCP2Mat)
    pi2Psigma = pi2PsigmaGen(pi2, CCP1Mat)
    eP1 = eulergamma - CCP1LogTransform(CCP1)
    eP2 = eulergamma - CCP2LogTransform(CCP2)

    inv_mat = np.linalg.inv(np.eye(8) - beta * fPsigma)
    ExanteV1 = inv_mat @ np.sum(CCP1Mat * (pi1Psigma + eP1), axis=1)
    ExanteV2 = inv_mat @ np.sum(CCP2Mat * (pi2Psigma + eP2), axis=1)

    # Step 3
    fP_a1 = fP_a1given(TransitionMat, CCP2)
    fP_a2 = fP_a2given(TransitionMat, CCP1)

    future1 = np.column_stack([fP_a1[k] @ ExanteV1 for k in range(3)])
    NewSigmaSeed1 = (pi1Psigma + beta * future1) * CCP1Adjuster
    NewSigmaDeno1 = np.sum(np.exp(NewSigmaSeed1), axis=1) - 1.0
    NewSigma1 = np.exp(NewSigmaSeed1) / NewSigmaDeno1.reshape(-1, 1)
    CCP1UpdatedMat = NewSigma1 * CCP1Adjuster
    CCP1Updated = CCP1UpdatedMat[:, 1]

    future2 = np.column_stack([fP_a2[k] @ ExanteV2 for k in range(3)])
    NewSigmaSeed2 = (pi2Psigma + beta * future2) * CCP2Adjuster
    NewSigmaDeno2 = np.sum(np.exp(NewSigmaSeed2), axis=1) - 1.0
    NewSigma2 = np.exp(NewSigmaSeed2) / NewSigmaDeno2.reshape(-1, 1)
    CCP2UpdatedMat = NewSigma2 * CCP2Adjuster
    CCP2Updated = CCP2UpdatedMat[:, 1]

    # Step 4
    fPsigma = fP(TransitionMat, CCP1Updated, CCP2Updated)
    pi1Psigma = pi1PsigmaGen(pi1, CCP2UpdatedMat)
    pi2Psigma = pi2PsigmaGen(pi2, CCP1UpdatedMat)
    eP1 = eulergamma - CCP1LogTransform(CCP1Updated)
    eP2 = eulergamma - CCP2LogTransform(CCP2Updated)
    inv_mat = np.linalg.inv(np.eye(8) - beta * fPsigma)
    ExanteV1Updated = inv_mat @ np.sum(CCP1UpdatedMat * (pi1Psigma + eP1), axis=1)
    ExanteV2Updated = inv_mat @ np.sum(CCP2UpdatedMat * (pi2Psigma + eP2), axis=1)

    # Step 5: iterate until convergence
    DiffExanteV = np.sum((ExanteV1Updated - ExanteV1)**2 + (ExanteV2Updated - ExanteV2)**2)

    iteration = 0
    while DiffExanteV > tol:
        iteration += 1
        CCP1 = CCP1Updated.copy()
        CCP2 = CCP2Updated.copy()
        ExanteV1 = ExanteV1Updated.copy()
        ExanteV2 = ExanteV2Updated.copy()

        fP_a1 = fP_a1given(TransitionMat, CCP2)
        fP_a2 = fP_a2given(TransitionMat, CCP1)

        future1 = np.column_stack([fP_a1[k] @ ExanteV1 for k in range(3)])
        NewSigmaSeed1 = (pi1Psigma + beta * future1) * CCP1Adjuster
        NewSigmaDeno1 = np.sum(np.exp(NewSigmaSeed1), axis=1) - 1.0
        NewSigma1 = np.exp(NewSigmaSeed1) / NewSigmaDeno1.reshape(-1, 1)
        CCP1UpdatedMat = NewSigma1 * CCP1Adjuster
        CCP1Updated = CCP1UpdatedMat[:, 1]

        future2 = np.column_stack([fP_a2[k] @ ExanteV2 for k in range(3)])
        NewSigmaSeed2 = (pi2Psigma + beta * future2) * CCP2Adjuster
        NewSigmaDeno2 = np.sum(np.exp(NewSigmaSeed2), axis=1) - 1.0
        NewSigma2 = np.exp(NewSigmaSeed2) / NewSigmaDeno2.reshape(-1, 1)
        CCP2UpdatedMat = NewSigma2 * CCP2Adjuster
        CCP2Updated = CCP2UpdatedMat[:, 1]

        fPsigma = fP(TransitionMat, CCP1Updated, CCP2Updated)
        pi1Psigma = pi1PsigmaGen(pi1, CCP2UpdatedMat)
        pi2Psigma = pi2PsigmaGen(pi2, CCP1UpdatedMat)
        eP1 = eulergamma - CCP1LogTransform(CCP1Updated)
        eP2 = eulergamma - CCP2LogTransform(CCP2Updated)
        inv_mat = np.linalg.inv(np.eye(8) - beta * fPsigma)
        ExanteV1Updated = inv_mat @ np.sum(CCP1UpdatedMat * (pi1Psigma + eP1), axis=1)
        ExanteV2Updated = inv_mat @ np.sum(CCP2UpdatedMat * (pi2Psigma + eP2), axis=1)

        DiffExanteV = np.sum((ExanteV1Updated - ExanteV1)**2 + (ExanteV2Updated - ExanteV2)**2)

    print(f"MPE収束: {iteration} 回の反復, 差分 = {DiffExanteV:.2e}")
    return CCP1UpdatedMat, CCP2UpdatedMat, ExanteV1Updated, ExanteV2Updated


print("MPEソルバー定義完了")
"""))

# ============================================================
# Cell 11: Compute Equilibrium CCP (sub_2_DGP part 1)
# ============================================================
cells.append(nbf.v4.new_markdown_cell(
    "## 均衡CCPの計算\n"
    "\n"
    "真のパラメータの下で均衡CCPを求める。"
))

cells.append(nbf.v4.new_code_cell("""\
start_time = time.time()
CCP1UpdatedMat, CCP2UpdatedMat, ExanteV1, ExanteV2 = f_MPE(
    TransitionMat, pi1, pi2, beta, eulergamma, CCP1Adjuster, CCP2Adjuster
)
elapsed = time.time() - start_time
print(f"計算時間: {elapsed:.3f} 秒")

# 均衡CCP行列の表示
state_labels = ['G00', 'G01', 'G10', 'G11', 'B00', 'B01', 'B10', 'B11']
eq_ccp_df = pd.DataFrame(
    np.hstack([CCP1UpdatedMat, CCP2UpdatedMat]),
    columns=['a1=-1', 'a1=0', 'a1=1', 'a2=-1', 'a2=0', 'a2=1'],
    index=state_labels
)
print("\\n均衡CCP:")
print(eq_ccp_df.to_string())
"""))

# ============================================================
# Cell 12: Simulate Synthetic Data (sub_2_DGP part 2)
# ============================================================
cells.append(nbf.v4.new_markdown_cell(
    "## 疑似データの生成\n"
    "\n"
    "均衡CCPに基づいて、500市場 x 50期間の疑似データを生成する。"
))

cells.append(nbf.v4.new_code_cell("""\
np.random.seed(2023)

NumSimMarkets = 500
NumSimPeriods = 50
NumSimFirms = 2

InitialState = np.random.randint(1, 9, size=NumSimMarkets)
RandomNumbers = np.random.uniform(0, 1, size=(NumSimMarkets, NumSimPeriods, NumSimFirms + 1))

total_rows = NumSimMarkets * NumSimPeriods
FakeData = np.zeros((total_rows, 8))

state_to_demand_n1_n2 = {
    1: (1, 0, 0), 2: (1, 0, 1), 3: (1, 1, 0), 4: (1, 1, 1),
    5: (2, 0, 0), 6: (2, 0, 1), 7: (2, 1, 0), 8: (2, 1, 1),
}
demand_n1_n2_to_state = {v: k for k, v in state_to_demand_n1_n2.items()}

print("疑似データ生成中...")
start_time = time.time()

for m in range(NumSimMarkets):
    for t in range(NumSimPeriods):
        row_idx = m * NumSimPeriods + t
        FakeData[row_idx, 0] = m + 1
        FakeData[row_idx, 1] = t + 1

        if t == 0:
            s = InitialState[m]
            FakeData[row_idx, 2] = s
            demand, n1, n2 = state_to_demand_n1_n2[s]
            FakeData[row_idx, 3] = demand
            FakeData[row_idx, 4] = n1
            FakeData[row_idx, 5] = n2
        else:
            prev_idx = row_idx - 1
            sprev = int(FakeData[prev_idx, 2])
            a1prev = int(FakeData[prev_idx, 6])
            a2prev = int(FakeData[prev_idx, 7])

            if sprev <= 4:
                if RandomNumbers[m, t, 2] < TransitionMat[0, 0]:
                    demand_now = 1
                else:
                    demand_now = 2
            else:
                if RandomNumbers[m, t, 2] < TransitionMat[1, 1]:
                    demand_now = 2
                else:
                    demand_now = 1

            n1_now = int(FakeData[prev_idx, 4]) + a1prev
            n2_now = int(FakeData[prev_idx, 5]) + a2prev
            s = demand_n1_n2_to_state[(demand_now, n1_now, n2_now)]
            FakeData[row_idx, 2] = s
            FakeData[row_idx, 3] = demand_now
            FakeData[row_idx, 4] = n1_now
            FakeData[row_idx, 5] = n2_now

        s = int(FakeData[row_idx, 2])
        s_idx = s - 1
        n1_cur = int(FakeData[row_idx, 4])
        n2_cur = int(FakeData[row_idx, 5])

        if RandomNumbers[m, t, 0] > CCP1UpdatedMat[s_idx, 1]:
            if n1_cur == 0:
                FakeData[row_idx, 6] = 1
            else:
                FakeData[row_idx, 6] = -1

        if RandomNumbers[m, t, 1] > CCP2UpdatedMat[s_idx, 1]:
            if n2_cur == 0:
                FakeData[row_idx, 7] = 1
            else:
                FakeData[row_idx, 7] = -1

elapsed = time.time() - start_time
print(f"疑似データ生成完了: {elapsed:.1f} 秒, shape = {FakeData.shape}")
"""))

# ============================================================
# Cell 13: Load Matlab FakeData
# ============================================================
cells.append(nbf.v4.new_markdown_cell(
    "## Matlab生成のFakeDataを読み込み\n"
    "\n"
    "紙面における点推定値を同じものを得るために、Matlabで生成したFakeDataを利用する。"
))

cells.append(nbf.v4.new_code_cell("""\
isUseMatlabData = 1

if isUseMatlabData == 1:
    matlab_data_path = data_dir / 'FakeData_Matlab.csv'
    try:
        FakeDataMatlab = pd.read_csv(matlab_data_path, header=None).values
        print(f"FakeData_Matlab.csv 読み込み成功: shape = {FakeDataMatlab.shape}")
        FakeData = FakeDataMatlab
    except FileNotFoundError:
        print("FakeData_Matlab.csv が見つからないため、Python生成データを使用")

print(f"使用するデータ: shape = {FakeData.shape}")
print(f"先頭5行:\\n{FakeData[:5]}")
"""))

# ============================================================
# Cell 14: Step 1 - Estimate CCP and Transition from Data
# ============================================================
cells.append(nbf.v4.new_markdown_cell(
    "## Step 1: データからCCPと遷移確率を推定"
))

cells.append(nbf.v4.new_code_cell("""\
# --- CCP の推定 ---
EstimatedCCP1 = np.zeros(8)
EstimatedCCP2 = np.zeros(8)

for s in range(1, 9):  # 1-indexed states
    sub = FakeData[FakeData[:, 2] == s]
    n_total = len(sub)
    n_a1_zero = np.sum(sub[:, 6] == 0)
    n_a2_zero = np.sum(sub[:, 7] == 0)
    EstimatedCCP1[s-1] = n_a1_zero / n_total
    EstimatedCCP2[s-1] = n_a2_zero / n_total

print("推定CCP1 (P(a1=0|s)):")
for s in range(8):
    print(f"  State {s+1}: {EstimatedCCP1[s]:.6f}")

print("\\n推定CCP2 (P(a2=0|s)):")
for s in range(8):
    print(f"  State {s+1}: {EstimatedCCP2[s]:.6f}")

# --- 遷移確率の推定 ---
EstimatedTransition = np.zeros((2, 2))

n_rows = len(FakeData)
lag_demand = np.zeros(n_rows)
lag_demand[1:] = FakeData[:-1, 3]  # 1期ラグの景気状態

# t=1 のデータを除外
mask = FakeData[:, 1] != 1
data_with_lag = np.column_stack([FakeData[mask], lag_demand[mask]])

for z in range(1, 3):
    sub_z = data_with_lag[data_with_lag[:, 3] == z]     # 今期がz
    sub_z_prev = sub_z[sub_z[:, 8] == z]                # 前期もz
    EstimatedTransition[z-1, z-1] = len(sub_z_prev) / len(sub_z)
    EstimatedTransition[z-1, 2-z] = 1 - EstimatedTransition[z-1, z-1]

print("\\n推定遷移確率行列:")
print(EstimatedTransition)
print(f"\\n真の遷移行列との差:")
print(f"  P(G|G): 真={0.7:.4f}, 推定={EstimatedTransition[0,0]:.6f}")
print(f"  P(B|B): 真={0.6:.4f}, 推定={EstimatedTransition[1,1]:.6f}")
"""))

# ============================================================
# Cell 15: Step 2 - AM Estimation (main estimation)
# ============================================================
cells.append(nbf.v4.new_markdown_cell(
    "## Step 2: AM推定量によるパラメータ推定\n"
    "\n"
    "Aguirregabiria-Mira (2007) のCCP反復アルゴリズム:\n"
    "1. 推定CCPを初期値とする\n"
    "2. CCPを所与として擬似尤度を最大化 -> パラメータを推定\n"
    "3. 推定パラメータからCCPを更新\n"
    "4. CCPが収束するまで2-3を繰り返す"
))

cells.append(nbf.v4.new_code_cell("""\
# 初期パラメータの設定
# InitialParameters1: 真のパラメータに近い初期値
InitialParameters = np.array([
    0.3,    # 企業1のベース利潤
    0.3,    # 企業2のベース利潤
    -0.25,  # 顧客収奪効果
    0.4,    # 景気が良い時の追加的利潤
    -0.12,  # 退出のためのコスト
    -2.2    # 参入のためのコスト
])

# CCPの初期値 (推定されたCCPをそのまま使う)
ccp1 = EstimatedCCP1.copy()
ccp2 = EstimatedCCP2.copy()

# データから観察された行動 (N, 3): [action1, action2, state(1-indexed)]
Actions = FakeData[:, [6, 7, 2]]

# 5パラメータの初期値: [theta1, theta2, rival, z_good, entry_cost]
# (退出コストは0に正規化)
initial = np.array([InitialParameters[0], InitialParameters[1],
                    InitialParameters[2], InitialParameters[3],
                    InitialParameters[5]])

print("AM推定を開始...")
print(f"初期パラメータ (5-dim): {initial}")
print(f"初期CCP1: {np.round(ccp1, 4)}")
print(f"初期CCP2: {np.round(ccp2, 4)}")

start_time = time.time()

for i in range(1, 10001):
    # 目的関数: 5パラメータ -> 10パラメータに展開して擬似尤度を計算
    def obj(x):
        param10 = np.array([x[0], x[2], x[3], 0, x[4],
                            x[1], x[2], x[3], 0, x[4]])
        return obj_lik(param10, ccp1, ccp2, EstimatedTransition, beta, Actions)

    # 擬似尤度を最大化 (minimizeに負号を付ける)
    sol = minimize(lambda x: -obj(x), initial, method='Nelder-Mead',
                   options={'maxiter': 10000, 'xatol': 1e-8, 'fatol': 1e-8})

    # 推定パラメータから10次元パラメータを構築
    param10 = np.array([sol.x[0], sol.x[2], sol.x[3], 0, sol.x[4],
                        sol.x[1], sol.x[2], sol.x[3], 0, sol.x[4]])

    # 推定パラメータのもとでCCPを更新
    newCCP = CCP_to_Value_to_prediction(param10, ccp1, ccp2, EstimatedTransition, beta)

    # 収束判定
    diff1 = np.dot(newCCP[:, 0] - ccp1, newCCP[:, 0] - ccp1)
    diff2 = np.dot(newCCP[:, 1] - ccp2, newCCP[:, 1] - ccp2)
    print(f"  Iteration {i}: LL={-sol.fun:.4f}, CCP diff1={diff1:.2e}, diff2={diff2:.2e}")

    if check_convergence(newCCP, np.column_stack([ccp1, ccp2]), tol=1e-6):
        print(f"\\nCCP収束! ({i} 回の反復)")
        break

    # CCPを更新
    ccp1 = newCCP[:, 0].copy()
    ccp2 = newCCP[:, 1].copy()

    # 次のループの初期値を更新
    initial = sol.x.copy()

elapsed = time.time() - start_time
print(f"\\nAM推定完了: {elapsed:.1f} 秒")
print(f"\\n推定されたパラメータ (5-dim):")
print(f"  theta1 (firm1 base profit) = {sol.x[0]:.6f}")
print(f"  theta2 (firm2 base profit) = {sol.x[1]:.6f}")
print(f"  rival effect               = {sol.x[2]:.6f}")
print(f"  good economy bonus         = {sol.x[3]:.6f}")
print(f"  entry cost                 = {sol.x[4]:.6f}")
print(f"\\n真のパラメータ:")
print(f"  theta1={Parameters[0]}, theta2={Parameters[1]}, rival={Parameters[2]}, z_good={Parameters[3]}, entry={Parameters[5]}")

print(f"\\n収束時のCCP:")
print(f"  CCP1 = {np.round(newCCP[:, 0], 6)}")
print(f"  CCP2 = {np.round(newCCP[:, 1], 6)}")
"""))

# ============================================================
# Cell 16: Estimation_AM_bootstrap function
# ============================================================
cells.append(nbf.v4.new_markdown_cell(
    "## Bootstrap用のAM推定関数\n"
    "\n"
    "BootstrapのリサンプルごとにAM推定を実行する関数を定義する。"
))

cells.append(nbf.v4.new_code_cell("""\
def Estimation_AM_bootstrap(FakeData_boot, beta):
    \"\"\"Aguirregabiria and Mira の推定をBootstrapサンプルに対して実行する

    Parameters
    ----------
    FakeData_boot : (N, 8) ブートストラップサンプル
    beta : float 割引因子

    Returns
    -------
    (EstCCP, EstTransition, est_params) or None (収束しない場合)
    EstCCP : (8, 2) [EstimatedCCP1, EstimatedCCP2]
    EstTransition : (2, 2) 推定遷移確率行列
    est_params : (5,) 推定パラメータ
    \"\"\"
    # Step 1: CCP と遷移確率の推定
    EstCCP1 = np.zeros(8)
    EstCCP2 = np.zeros(8)

    for s in range(1, 9):
        sub = FakeData_boot[FakeData_boot[:, 2] == s]
        n_total = len(sub)
        if n_total == 0:
            EstCCP1[s-1] = 0.5
            EstCCP2[s-1] = 0.5
            continue
        EstCCP1[s-1] = np.sum(sub[:, 6] == 0) / n_total
        EstCCP2[s-1] = np.sum(sub[:, 7] == 0) / n_total

    EstTransition = np.zeros((2, 2))
    n_rows = len(FakeData_boot)
    lag_demand = np.zeros(n_rows)
    lag_demand[1:] = FakeData_boot[:-1, 3]
    mask = FakeData_boot[:, 1] != 1
    data_with_lag = np.column_stack([FakeData_boot[mask], lag_demand[mask]])

    for z in range(1, 3):
        sub_z = data_with_lag[data_with_lag[:, 3] == z]
        sub_zz = sub_z[sub_z[:, 8] == z]
        if len(sub_z) > 0:
            EstTransition[z-1, z-1] = len(sub_zz) / len(sub_z)
        else:
            EstTransition[z-1, z-1] = 0.5
        EstTransition[z-1, 2-z] = 1 - EstTransition[z-1, z-1]

    # Step 2: 初期パラメータ
    InitParams = np.array([0.3, 0.3, -0.25, 0.4, -0.12, -2.2])
    init_x = np.array([InitParams[0], InitParams[1],
                       InitParams[2], InitParams[3], InitParams[5]])

    ccp1_b = EstCCP1.copy()
    ccp2_b = EstCCP2.copy()

    Actions_b = FakeData_boot[:, [6, 7, 2]]

    # Stuck防止用の乱数
    rng_perturb = np.random.RandomState(123456)
    random_draw = rng_perturb.uniform(0.75, 1.25, size=8 * 2 * 100)

    # Step 3: AM estimator loop
    for i in range(1, 101):
        if i == 100:
            return None  # 収束しなかった場合

        def obj_boot(x):
            p10 = np.array([x[0], x[2], x[3], 0, x[4],
                            x[1], x[2], x[3], 0, x[4]])
            return obj_lik(p10, ccp1_b, ccp2_b, EstTransition, beta, Actions_b)

        sol_b = minimize(lambda x: -obj_boot(x), init_x, method='Nelder-Mead',
                         options={'maxiter': 10000, 'xatol': 1e-8, 'fatol': 1e-8})

        p10 = np.array([sol_b.x[0], sol_b.x[2], sol_b.x[3], 0, sol_b.x[4],
                        sol_b.x[1], sol_b.x[2], sol_b.x[3], 0, sol_b.x[4]])
        newCCP_b = CCP_to_Value_to_prediction(p10, ccp1_b, ccp2_b, EstTransition, beta)

        if check_convergence(newCCP_b, np.column_stack([ccp1_b, ccp2_b]), tol=1e-6):
            break

        ccp1_b = newCCP_b[:, 0].copy()
        ccp2_b = newCCP_b[:, 1].copy()
        init_x = sol_b.x.copy()

        # 10の倍数のとき、CCPをランダムに揺らがせる (Stuck防止)
        if i % 10 == 0:
            idx_base = (i // 10 - 1) * 16
            ccp1_b = ccp1_b * random_draw[idx_base:idx_base + 8]
            ccp2_b = ccp2_b * random_draw[idx_base + 8:idx_base + 16]
            # CCPが [0, 1] の範囲内に収まるようクリップ
            ccp1_b = np.clip(ccp1_b, 0.01, 0.99)
            ccp2_b = np.clip(ccp2_b, 0.01, 0.99)

    return (np.column_stack([EstCCP1, EstCCP2]), EstTransition, sol_b.x)


print("Estimation_AM_bootstrap 関数定義完了")
"""))

# ============================================================
# Cell 17: Bootstrap SE
# ============================================================
cells.append(nbf.v4.new_markdown_cell(
    "## Bootstrapによる標準誤差の計算\n"
    "\n"
    "マーケット単位でリサンプリングを行い(500市場)、\n"
    "100回のBootstrapで標準誤差を計算する。"
))

cells.append(nbf.v4.new_code_cell("""\
np.random.seed(2023)

numBootSample = 100
NumSimPeriodsB = int(FakeData[:, 1].max())
NumSimMarketsB = int(FakeData[:, 0].max())

# 各Bootstrapサンプルで用いるマーケットのインデックスを生成
bootindex = np.random.randint(1, NumSimMarketsB + 1,
                              size=(NumSimMarketsB, numBootSample))

# 結果保存用
bootresult_transition = np.full((2, numBootSample), np.nan)
bootresult_CCP1 = np.full((8, numBootSample), np.nan)
bootresult_CCP2 = np.full((8, numBootSample), np.nan)
bootresult_payoff = np.full((5, numBootSample), np.nan)

print(f"Bootstrap開始: {numBootSample} サンプル")
start_time = time.time()

for b in range(numBootSample):
    if (b + 1) % 10 == 0 or b == 0:
        elapsed_b = time.time() - start_time
        print(f"  Bootstrap {b+1}/{numBootSample} (経過: {elapsed_b:.1f}秒)")

    # Bootstrapサンプルの構築 (マーケット単位でリサンプリング)
    boot_sample_list = []
    for m in range(NumSimMarketsB):
        mk = bootindex[m, b]
        boot_sample_list.append(FakeData[FakeData[:, 0] == mk])
    bootsample = np.vstack(boot_sample_list)

    # AM推定を実行
    output = Estimation_AM_bootstrap(bootsample, beta)

    if output is not None:
        bootresult_CCP1[:, b] = output[0][:, 0]
        bootresult_CCP2[:, b] = output[0][:, 1]
        bootresult_transition[:, b] = np.diag(output[1])
        bootresult_payoff[:, b] = output[2]
    # output が None の場合は NaN のまま

elapsed = time.time() - start_time
print(f"\\nBootstrap完了: {elapsed:.1f} 秒")

# NAの処理
na_columns = np.any(np.isnan(bootresult_CCP1), axis=0)
n_na = np.sum(na_columns)
print(f"収束しなかったサンプル数: {n_na}")

if n_na > 30:
    print("WARNING: NAが多すぎます (> 30)")
else:
    # NAを除外
    valid_cols = ~na_columns
    bootresult_CCP1_valid = bootresult_CCP1[:, valid_cols]
    bootresult_CCP2_valid = bootresult_CCP2[:, valid_cols]
    bootresult_transition_valid = bootresult_transition[:, valid_cols]
    bootresult_payoff_valid = bootresult_payoff[:, valid_cols]

    # 最大100列に制限
    n_valid = min(bootresult_CCP1_valid.shape[1], 100)
    bootresult_CCP1_valid = bootresult_CCP1_valid[:, :n_valid]
    bootresult_CCP2_valid = bootresult_CCP2_valid[:, :n_valid]
    bootresult_transition_valid = bootresult_transition_valid[:, :n_valid]
    bootresult_payoff_valid = bootresult_payoff_valid[:, :n_valid]

    print(f"有効なBootstrapサンプル数: {n_valid}")
"""))

# ============================================================
# Cell 18: Summary Tables
# ============================================================
cells.append(nbf.v4.new_markdown_cell(
    "## 結果の表示と保存\n"
    "\n"
    "推定結果をTab11_2, Tab11_3, Tab11_4 として表示・保存する。"
))

cells.append(nbf.v4.new_code_cell("""\
# --- Tab 11.2: 遷移確率の推定結果 ---
true_trans = np.array([0.7, 0.6])
se_trans = np.std(bootresult_transition_valid, axis=1, ddof=0)

tab11_2 = pd.DataFrame({
    'Transition Probability (GG and BB): True': true_trans,
    'Estimated': np.diag(EstimatedTransition),
    'SE': se_trans
})
tab11_2.index = range(1, 3)

print("=" * 70)
print("表 11.2: 遷移確率の推定結果")
print("=" * 70)
print(tab11_2.to_string())
tab11_2.to_csv(output_dir / 'Tab11_2_Transition.csv')
print("\\n保存: output/Tab11_2_Transition.csv")
"""))

cells.append(nbf.v4.new_code_cell("""\
# --- Tab 11.3: CCP (企業1) の推定結果 ---
true_ccp1 = CCP1UpdatedMat[:, 1]
se_ccp1 = np.std(bootresult_CCP1_valid, axis=1, ddof=0)

tab11_3_firm1 = pd.DataFrame({
    'CCP for firm 1: True': true_ccp1,
    'Estimated': EstimatedCCP1,
    'SE': se_ccp1
})
tab11_3_firm1.index = range(1, 9)

print("=" * 70)
print("表 11.3: CCP の推定結果 (企業1)")
print("=" * 70)
print(tab11_3_firm1.to_string())
tab11_3_firm1.to_csv(output_dir / 'Tab11_3_CCP_firm1.csv')
print("\\n保存: output/Tab11_3_CCP_firm1.csv")

# --- Tab 11.3: CCP (企業2) の推定結果 ---
true_ccp2 = CCP2UpdatedMat[:, 1]
se_ccp2 = np.std(bootresult_CCP2_valid, axis=1, ddof=0)

tab11_3_firm2 = pd.DataFrame({
    'CCP for firm 2: True': true_ccp2,
    'Estimated': EstimatedCCP2,
    'SE': se_ccp2
})
tab11_3_firm2.index = range(1, 9)

print("\\n" + "=" * 70)
print("表 11.3: CCP の推定結果 (企業2)")
print("=" * 70)
print(tab11_3_firm2.to_string())
tab11_3_firm2.to_csv(output_dir / 'Tab11_3_CCP_firm2.csv')
print("\\n保存: output/Tab11_3_CCP_firm2.csv")
"""))

cells.append(nbf.v4.new_code_cell("""\
# --- Tab 11.4: AM2007 利潤パラメータの推定結果 ---
# 真のパラメータ (5-dim): [theta1, theta2, rival, z_good, entry_cost]
true_params = np.array([Parameters[0], Parameters[1], Parameters[2],
                        Parameters[3], Parameters[5]])

# 正規化 (Aguirregabiria-Suzuki normalization)
# theta1_normalized = theta1 - (1-beta)/beta * theta5 (exit cost)
# theta2_normalized = theta2 - (1-beta)/beta * theta5
# theta5_normalized = 0 (exit cost normalized to 0)
# theta6_normalized = theta6 + theta5 (entry + exit)
# rival, z_good は変わらない
normalized_params = np.array([
    Parameters[0] - (1 - beta) / beta * Parameters[4],   # theta1 - (1-b)/b * exit
    Parameters[1] - (1 - beta) / beta * Parameters[4],   # theta2 - (1-b)/b * exit
    Parameters[2],                                         # rival (unchanged)
    Parameters[3],                                         # z_good (unchanged)
    Parameters[5] + Parameters[4]                          # entry + exit
])

# Bootstrap SE
se_payoff = np.std(bootresult_payoff_valid, axis=1, ddof=0)

tab11_4 = pd.DataFrame({
    'Payoff parameter: True': true_params,
    'Normalized true': normalized_params,
    'Estimated': sol.x,
    'SE': se_payoff
})
tab11_4.index = range(1, 6)

print("=" * 70)
print("表 11.4: AM2007 利潤パラメータの推定結果")
print("=" * 70)
print(tab11_4.to_string())

print("\\nパラメータの解釈:")
print("  Row 1: 企業1のベース利潤")
print("  Row 2: 企業2のベース利潤")
print("  Row 3: 顧客収奪効果 (ライバルの影響)")
print("  Row 4: 景気が良い時の追加的利潤")
print("  Row 5: 参入のためのコスト (正規化後: entry + exit)")
print(f"\\n注: 退出コストは0に正規化されている (theta_exit = {Parameters[4]})")
print(f"    正規化: theta1_norm = theta1 - (1-beta)/beta * theta_exit")
print(f"           = {Parameters[0]} - {(1-beta)/beta} * ({Parameters[4]}) = {normalized_params[0]}")

tab11_4.to_csv(output_dir / 'Tab11_4_AM2007.csv')
print("\\n保存: output/Tab11_4_AM2007.csv")
"""))

# ============================================================
# Cell 19: Summary
# ============================================================
cells.append(nbf.v4.new_markdown_cell("## まとめ"))

cells.append(nbf.v4.new_code_cell("""\
print("=" * 70)
print("第11章: AM (Aguirregabiria-Mira 2007) 推定 まとめ")
print("=" * 70)
print(f"状態空間: 8 状態 (2x2x2: 景気 x 企業1 x 企業2)")
print(f"各企業の行動: 3種類 (退出=-1, 現状維持=0, 参入=1)")
print(f"割引因子 beta = {beta}")
print(f"遷移行列: P(G|G)={TransitionMat[0,0]}, P(B|B)={TransitionMat[1,1]}")
print(f"\\n推定方法: Aguirregabiria-Mira (2007) CCP反復推定量")
print(f"  - 5パラメータを推定 (退出コスト=0に正規化)")
print(f"  - Nelder-Mead法による擬似尤度最大化")
print(f"  - CCP反復で収束まで繰り返し")
print(f"\\nBootstrap: {numBootSample} サンプル (マーケット単位リサンプリング)")
print(f"\\n出力ファイル:")
print(f"  - Tab11_2_Transition.csv: 遷移確率の推定結果")
print(f"  - Tab11_3_CCP_firm1.csv: 企業1のCCP推定結果")
print(f"  - Tab11_3_CCP_firm2.csv: 企業2のCCP推定結果")
print(f"  - Tab11_4_AM2007.csv:    利潤パラメータの推定結果")
print("\\n完了")
"""))

nb.cells = cells
nbf.write(nb, 'main_ch11_AM.ipynb')
print("Generated: main_ch11_AM.ipynb")
