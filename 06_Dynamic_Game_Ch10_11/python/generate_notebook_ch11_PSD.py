"""Generate main_ch11_PSD.ipynb for Chapter 11: P-SD (Pesendorfer-Schmidt-Dengler) Estimation."""
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
    "# 第11章 Pesendorfer-Schmidt-Dengler (P-SD) によるパラメータ推定\n"
    "\n"
    "2企業の参入退出の動的ゲームにおいて、\n"
    "Pesendorfer and Schmidt-Dengler (2008) の方法でペイオフパラメータを推定する。\n"
    "\n"
    "- Step 1: データからCCPと遷移確率を推定\n"
    "- Step 2: P-SD推定量（最小二乗法）によるパラメータ推定\n"
    "- Bootstrapによる標準誤差の計算"
))

# ============================================================
# Cell 2: Setup
# ============================================================
cells.append(nbf.v4.new_code_cell("""\
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.optimize import minimize
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
# Cell 3: Parameters
# ============================================================
cells.append(nbf.v4.new_markdown_cell("## パラメータの設定"))

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

# シミュレーション設定
NumSimMarkets = 500   # 市場の数
NumSimPeriods = 50    # 期間の数

print("パラメータ設定完了")
print(f"beta = {beta}")
print(f"遷移行列:\\n{TransitionMat}")
print(f"Parameters = {Parameters}")
"""))

# ============================================================
# Cell 4: Helper Functions - Profit
# ============================================================
cells.append(nbf.v4.new_markdown_cell("## ヘルパー関数の定義"))

cells.append(nbf.v4.new_code_cell("""\
def pi1gen(theta):
    \"\"\"企業1の利潤行列を生成する (8x3)
    行: 状態 (G00, G01, G10, G11, B00, B01, B10, B11)
    列: 行動 (a1=-1, a1=0, a1=1)
    \"\"\"
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
    invdiv = np.array([theta[3], 0, theta[4]])
    output = np.tile(base.reshape(-1, 1), (1, 3)) + np.tile(invdiv.reshape(1, -1), (8, 1))
    return output


def pi2gen(theta):
    \"\"\"企業2の利潤行列を生成する (8x3)\"\"\"
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
    invdiv = np.array([theta[8], 0, theta[9]])
    output = np.tile(base.reshape(-1, 1), (1, 3)) + np.tile(invdiv.reshape(1, -1), (8, 1))
    return output


def CCP1Transform(x):
    \"\"\"CCP1ベクトル (8,) を CCP1行列 (8x3) に変換する\"\"\"
    output = np.array([
        [0,        x[0], 1-x[0]],
        [0,        x[1], 1-x[1]],
        [1-x[2],   x[2], 0      ],
        [1-x[3],   x[3], 0      ],
        [0,        x[4], 1-x[4]],
        [0,        x[5], 1-x[5]],
        [1-x[6],   x[6], 0      ],
        [1-x[7],   x[7], 0      ],
    ])
    return output


def CCP2Transform(x):
    \"\"\"CCP2ベクトル (8,) を CCP2行列 (8x3) に変換する\"\"\"
    output = np.array([
        [0,        x[0], 1-x[0]],
        [1-x[1],   x[1], 0      ],
        [0,        x[2], 1-x[2]],
        [1-x[3],   x[3], 0      ],
        [0,        x[4], 1-x[4]],
        [1-x[5],   x[5], 0      ],
        [0,        x[6], 1-x[6]],
        [1-x[7],   x[7], 0      ],
    ])
    return output


def CCP1LogTransform(x):
    \"\"\"CCP1ベクトル (8,) を log CCP1行列 (8x3) に変換する\"\"\"
    output = np.array([
        [0,           np.log(x[0]), np.log(1-x[0])],
        [0,           np.log(x[1]), np.log(1-x[1])],
        [np.log(1-x[2]), np.log(x[2]), 0           ],
        [np.log(1-x[3]), np.log(x[3]), 0           ],
        [0,           np.log(x[4]), np.log(1-x[4])],
        [0,           np.log(x[5]), np.log(1-x[5])],
        [np.log(1-x[6]), np.log(x[6]), 0           ],
        [np.log(1-x[7]), np.log(x[7]), 0           ],
    ])
    return output


def CCP2LogTransform(x):
    \"\"\"CCP2ベクトル (8,) を log CCP2行列 (8x3) に変換する\"\"\"
    output = np.array([
        [0,           np.log(x[0]), np.log(1-x[0])],
        [np.log(1-x[1]), np.log(x[1]), 0           ],
        [0,           np.log(x[2]), np.log(1-x[2])],
        [np.log(1-x[3]), np.log(x[3]), 0           ],
        [0,           np.log(x[4]), np.log(1-x[4])],
        [np.log(1-x[5]), np.log(x[5]), 0           ],
        [0,           np.log(x[6]), np.log(1-x[6])],
        [np.log(1-x[7]), np.log(x[7]), 0           ],
    ])
    return output


print("利潤関数・CCP変換関数 定義完了")
"""))

# ============================================================
# Cell 5: Helper Functions - Transition Matrices
# ============================================================
cells.append(nbf.v4.new_code_cell("""\
def fP(Matrix1, Vec1, Vec2):
    \"\"\"CCPの下での状態遷移行列 F^{P,sigma} を構築する (8x8)\"\"\"
    TempMat0 = np.kron(Matrix1, np.ones((4, 4)))

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


print("遷移行列・期待利潤関数 定義完了")
"""))

# ============================================================
# Cell 6: MPE Solver
# ============================================================
cells.append(nbf.v4.new_markdown_cell("## MPEソルバーの定義"))

cells.append(nbf.v4.new_code_cell("""\
def f_MPE(TransitionMat, pi1, pi2, beta, eulergamma, CCP1Adjuster, CCP2Adjuster, tol=1e-12):
    \"\"\"マルコフ完全均衡（MPE）を固定点反復で計算する\"\"\"
    # Step 1: CCPの初期値
    CCP1 = np.full(8, 0.5)
    CCP2 = np.full(8, 0.5)

    CCP1Mat = CCP1Transform(CCP1)
    CCP2Mat = CCP2Transform(CCP2)

    # Step 2: 事前価値関数の計算
    fPsigma = fP(TransitionMat, CCP1, CCP2)
    pi1Psigma = pi1PsigmaGen(pi1, CCP2Mat)
    pi2Psigma = pi2PsigmaGen(pi2, CCP1Mat)
    eP1 = eulergamma - CCP1LogTransform(CCP1)
    eP2 = eulergamma - CCP2LogTransform(CCP2)

    inv_mat = np.linalg.inv(np.eye(8) - beta * fPsigma)
    ExanteV1 = inv_mat @ np.sum(CCP1Mat * (pi1Psigma + eP1), axis=1)
    ExanteV2 = inv_mat @ np.sum(CCP2Mat * (pi2Psigma + eP2), axis=1)

    # Step 3: CCPの更新
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

    # Step 4: 事前価値関数の再計算
    fPsigma = fP(TransitionMat, CCP1Updated, CCP2Updated)
    pi1Psigma = pi1PsigmaGen(pi1, CCP2UpdatedMat)
    pi2Psigma = pi2PsigmaGen(pi2, CCP1UpdatedMat)
    eP1 = eulergamma - CCP1LogTransform(CCP1Updated)
    eP2 = eulergamma - CCP2LogTransform(CCP2Updated)

    inv_mat = np.linalg.inv(np.eye(8) - beta * fPsigma)
    ExanteV1Updated = inv_mat @ np.sum(CCP1UpdatedMat * (pi1Psigma + eP1), axis=1)
    ExanteV2Updated = inv_mat @ np.sum(CCP2UpdatedMat * (pi2Psigma + eP2), axis=1)

    # Step 5: 収束するまで繰り返す
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
# Cell 7: CCP_to_Value_to_prediction
# ============================================================
cells.append(nbf.v4.new_markdown_cell("## CCP_to_Value_to_prediction 関数の定義\n"
    "\n"
    "CCPから事前の価値関数をインバージョンし、その価値関数とパラメータからCCPを予測する。"))

cells.append(nbf.v4.new_code_cell("""\
def CCP_to_Value_to_prediction(theta, CCP1, CCP2, TransitionMat, beta):
    \"\"\"CCPから価値関数をインバージョンし、CCPを予測する

    Parameters
    ----------
    theta : (10,) パラメータベクトル
    CCP1 : (8,) 企業1のCCPベクトル (P(a1=0))
    CCP2 : (8,) 企業2のCCPベクトル (P(a2=0))
    TransitionMat : (2,2) 景気の遷移行列
    beta : float 割引因子

    Returns
    -------
    output : (8,2) 予測されたCCP [CCP1Updated, CCP2Updated]
    \"\"\"
    # Step 1: Setup
    CCP1Adj = np.array([
        [0, 1, 1], [0, 1, 1], [1, 1, 0], [1, 1, 0],
        [0, 1, 1], [0, 1, 1], [1, 1, 0], [1, 1, 0],
    ], dtype=float)
    CCP2Adj = np.array([
        [0, 1, 1], [1, 1, 0], [0, 1, 1], [1, 1, 0],
        [0, 1, 1], [1, 1, 0], [0, 1, 1], [1, 1, 0],
    ], dtype=float)

    pi1_local = pi1gen(theta) * CCP1Adj
    pi2_local = pi2gen(theta) * CCP2Adj

    CCP1Mat = CCP1Transform(CCP1)
    CCP2Mat = CCP2Transform(CCP2)

    # Step 2: 価値関数のインバージョン
    fPsigma = fP(TransitionMat, CCP1, CCP2)

    pi1Psigma = pi1PsigmaGen(pi1_local, CCP2Mat)
    pi2Psigma = pi2PsigmaGen(pi2_local, CCP1Mat)

    eP1 = eulergamma - CCP1LogTransform(CCP1)
    eP2 = eulergamma - CCP2LogTransform(CCP2)

    inv_mat = np.linalg.inv(np.eye(8) - beta * fPsigma)
    ExanteV1 = inv_mat @ np.sum(CCP1Mat * (pi1Psigma + eP1), axis=1)
    ExanteV2 = inv_mat @ np.sum(CCP2Mat * (pi2Psigma + eP2), axis=1)

    # Step 3: CCPの予測
    fP_a1_list = fP_a1given(TransitionMat, CCP2)
    fP_a2_list = fP_a2given(TransitionMat, CCP1)

    future1 = np.column_stack([fP_a1_list[k] @ ExanteV1 for k in range(3)])
    NewSigmaSeed1 = (pi1Psigma + beta * future1) * CCP1Adj
    NewSigmaDeno1 = np.sum(np.exp(NewSigmaSeed1), axis=1) - 1.0
    NewSigma1 = np.exp(NewSigmaSeed1) / NewSigmaDeno1.reshape(-1, 1)
    CCP1UpdatedMat = NewSigma1 * CCP1Adj
    CCP1Updated = CCP1UpdatedMat[:, 1]

    future2 = np.column_stack([fP_a2_list[k] @ ExanteV2 for k in range(3)])
    NewSigmaSeed2 = (pi2Psigma + beta * future2) * CCP2Adj
    NewSigmaDeno2 = np.sum(np.exp(NewSigmaSeed2), axis=1) - 1.0
    NewSigma2 = np.exp(NewSigmaSeed2) / NewSigmaDeno2.reshape(-1, 1)
    CCP2UpdatedMat = NewSigma2 * CCP2Adj
    CCP2Updated = CCP2UpdatedMat[:, 1]

    output = np.column_stack([CCP1Updated, CCP2Updated])
    return output


print("CCP_to_Value_to_prediction 定義完了")
"""))

# ============================================================
# Cell 8: obj_fun and Estimation_PS_bootstrap
# ============================================================
cells.append(nbf.v4.new_markdown_cell("## 目的関数とBootstrap推定関数の定義"))

cells.append(nbf.v4.new_code_cell("""\
def obj_fun(param, EstimatedCCP1, EstimatedCCP2, TransitionProb, beta):
    \"\"\"P-SD推定の目的関数: 推定CCPと予測CCPの二乗和 (SSE)

    Parameters
    ----------
    param : (10,) パラメータベクトル
    EstimatedCCP1 : (8,) 推定された企業1のCCP
    EstimatedCCP2 : (8,) 推定された企業2のCCP
    TransitionProb : (2,2) 推定された遷移確率行列
    beta : float 割引因子

    Returns
    -------
    obj : float 目的関数の値 (SSE)
    \"\"\"
    output = CCP_to_Value_to_prediction(param, EstimatedCCP1, EstimatedCCP2, TransitionProb, beta)
    obj = np.sum((output[:, 0] - EstimatedCCP1)**2) + np.sum((output[:, 1] - EstimatedCCP2)**2)
    return obj


def Estimation_PS_bootstrap(FakeData_boot, beta):
    \"\"\"Bootstrapサンプルに対するP-SD推定

    Parameters
    ----------
    FakeData_boot : (N, 8) Bootstrapサンプル
    beta : float 割引因子

    Returns
    -------
    list: [CCP行列 (8,2), 遷移確率行列 (2,2), 推定パラメータ (5,), 目的関数値]
    \"\"\"
    # Step 1: CCPと遷移確率の推定
    EstCCP1 = np.zeros(8)
    EstCCP2 = np.zeros(8)

    for s in range(1, 9):  # 1-indexed states
        sub = FakeData_boot[FakeData_boot[:, 2] == s]
        n_total = len(sub)
        if n_total > 0:
            EstCCP1[s - 1] = np.sum(sub[:, 6] == 0) / n_total
            EstCCP2[s - 1] = np.sum(sub[:, 7] == 0) / n_total

    EstTrans = np.zeros((2, 2))
    n_rows = len(FakeData_boot)
    lag_demand = np.zeros(n_rows)
    lag_demand[1:] = FakeData_boot[:-1, 3]
    mask = FakeData_boot[:, 1] != 1
    data_with_lag = np.column_stack([FakeData_boot[mask], lag_demand[mask]])

    for z in range(1, 3):
        sub_z = data_with_lag[data_with_lag[:, 3] == z]
        sub_zz = sub_z[sub_z[:, 8] == z]
        if len(sub_z) > 0:
            EstTrans[z - 1, z - 1] = len(sub_zz) / len(sub_z)
    EstTrans[0, 1] = 1 - EstTrans[0, 0]
    EstTrans[1, 0] = 1 - EstTrans[1, 1]

    # Step 2: P-SD推定量 (Nelder-Mead)
    def obj_wrapper(x):
        param10 = np.array([x[0], x[2], x[3], 0.0, x[4],
                            x[1], x[2], x[3], 0.0, x[4]])
        return obj_fun(param10, EstCCP1, EstCCP2, EstTrans, beta)

    initial = np.array([0.3, 0.2, -0.27, 0.45, -2.1])
    sol = minimize(obj_wrapper, initial, method='Nelder-Mead')

    output = [
        np.column_stack([EstCCP1, EstCCP2]),
        EstTrans,
        sol.x,
        sol.fun
    ]
    return output


print("obj_fun, Estimation_PS_bootstrap 定義完了")
"""))

# ============================================================
# Cell 9: Compute Equilibrium + Fake Data
# ============================================================
cells.append(nbf.v4.new_markdown_cell(
    "## 1. 下準備: 均衡CCPの計算と疑似データの生成\n"
    "\n"
    "Ch10の均衡計算と疑似データ生成をインラインで実行する。"
))

cells.append(nbf.v4.new_code_cell("""\
# 真のパラメータで利潤行列を計算
pi1 = pi1gen(TrueParameterValues) * CCP1Adjuster
pi2 = pi2gen(TrueParameterValues) * CCP2Adjuster

# 均衡CCPの計算
start_time = time.time()
CCP1UpdatedMat, CCP2UpdatedMat, ExanteV1, ExanteV2 = f_MPE(
    TransitionMat, pi1, pi2, beta, eulergamma, CCP1Adjuster, CCP2Adjuster
)
elapsed = time.time() - start_time
print(f"均衡CCP計算時間: {elapsed:.3f} 秒")

# 均衡CCPの表示
state_labels = ['G00', 'G01', 'G10', 'G11', 'B00', 'B01', 'B10', 'B11']
print("\\n均衡CCP1 (P(a1=0 | s)):", np.round(CCP1UpdatedMat[:, 1], 6))
print("均衡CCP2 (P(a2=0 | s)):", np.round(CCP2UpdatedMat[:, 1], 6))
"""))

# ============================================================
# Cell 10: Simulate FakeData
# ============================================================
cells.append(nbf.v4.new_code_cell("""\
# 疑似データの生成
np.random.seed(2023)

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
        FakeData[row_idx, 0] = m + 1  # Market ID (1-indexed)
        FakeData[row_idx, 1] = t + 1  # Time (1-indexed)

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
# Cell 11: Load Matlab FakeData
# ============================================================
cells.append(nbf.v4.new_markdown_cell(
    "## 2. データの読み込み\n"
    "\n"
    "紙面の点推定値を再現するために、Matlabで生成したFakeDataを読み込む。"
))

cells.append(nbf.v4.new_code_cell("""\
# Matlabで生成されたFakeDataを読み込む
matlab_data_path = data_dir / 'FakeData_Matlab.csv'
try:
    FakeDataMatlab = pd.read_csv(matlab_data_path, header=None).values
    print(f"FakeData_Matlab.csv 読み込み成功: shape = {FakeDataMatlab.shape}")
    FakeData = FakeDataMatlab
except FileNotFoundError:
    print("FakeData_Matlab.csv が見つからないため、Python生成データを使用")

print(f"使用データ: shape = {FakeData.shape}")
print(f"先頭5行:\\n{FakeData[:5].astype(int)}")
"""))

# ============================================================
# Cell 12: Step 1 - Estimate CCP and Transition
# ============================================================
cells.append(nbf.v4.new_markdown_cell(
    "## 3. 推定 Step 1: CCPと遷移確率の推定\n"
    "\n"
    "疑似データからCCPと景気の遷移確率を推定する。"
))

cells.append(nbf.v4.new_code_cell("""\
# CCPの推定
EstimatedCCP1 = np.zeros(8)
EstimatedCCP2 = np.zeros(8)

for s in range(1, 9):  # 1-indexed states
    sub = FakeData[FakeData[:, 2] == s]
    n_total = len(sub)
    n_a1_zero = np.sum(sub[:, 6] == 0)
    n_a2_zero = np.sum(sub[:, 7] == 0)
    EstimatedCCP1[s - 1] = n_a1_zero / n_total
    EstimatedCCP2[s - 1] = n_a2_zero / n_total

print("推定されたCCP:")
print(np.column_stack([EstimatedCCP1, EstimatedCCP2]))

# 遷移確率の推定
EstimatedTransition = np.zeros((2, 2))

n_rows = len(FakeData)
lag_demand = np.zeros(n_rows)
lag_demand[1:] = FakeData[:-1, 3]
# t=1 を除外 (市場間のまたぎを避ける)
mask = FakeData[:, 1] != 1
data_with_lag = np.column_stack([FakeData[mask], lag_demand[mask]])

for z in range(1, 3):
    sub_z = data_with_lag[data_with_lag[:, 3] == z]
    sub_zz = sub_z[sub_z[:, 8] == z]
    EstimatedTransition[z - 1, z - 1] = len(sub_zz) / len(sub_z)
    EstimatedTransition[z - 1, 2 - z] = 1 - EstimatedTransition[z - 1, z - 1]

print("\\n推定されたzの遷移確率行列:")
print(EstimatedTransition)
"""))

# ============================================================
# Cell 13: Step 2-1 - Sanity Check 1
# ============================================================
cells.append(nbf.v4.new_markdown_cell(
    "## 4. 推定 Step 2: P-SD推定\n"
    "\n"
    "### Step 2-1: Sanity Check 1\n"
    "\n"
    "DGPにおける真のCCP、遷移確率、真のパラメータを入れて出てくる予測が、\n"
    "DGPのCCPと一致しているか確認する。"
))

cells.append(nbf.v4.new_code_cell("""\
# 真のCCP (均衡CCP)
CCP1true = CCP1UpdatedMat[:, 1]
CCP2true = CCP2UpdatedMat[:, 1]

# 真のパラメータと真のCCPで予測
output = CCP_to_Value_to_prediction(TrueParameterValues, CCP1true, CCP2true, TransitionMat, beta)

# 真のCCPとの差分 (ゼロに近いはず)
diff = output - np.column_stack([CCP1true, CCP2true])
print("Sanity Check 1: 真のCCPと予測CCPの差分")
print("(ゼロに近ければ正しく実装されている)")
print(diff)
print(f"\\n最大絶対差: {np.abs(diff).max():.2e}")
"""))

# ============================================================
# Cell 14: Step 2-2 - Sanity Check 2
# ============================================================
cells.append(nbf.v4.new_markdown_cell(
    "### Step 2-2: Sanity Check 2\n"
    "\n"
    "Aguirregabiria-Suzukiの正規化を施したパラメータでも同じ予測が得られることを確認する。"
))

cells.append(nbf.v4.new_code_cell("""\
# 正規化されたパラメータ
# theta4 (退出コスト) を0に正規化し、他のパラメータを調整
Normalized_TrueParam = np.array([
    Parameters[0] - ((1 - beta) / beta) * Parameters[4],  # 企業1のベース利潤
    Parameters[2],                                          # 顧客収奪効果
    Parameters[3],                                          # 好景気効果
    0.0,                                                    # 退出コスト (0に正規化)
    Parameters[5] + Parameters[4],                          # 参入コスト
    Parameters[1] - ((1 - beta) / beta) * Parameters[4],  # 企業2のベース利潤
    Parameters[2],                                          # 顧客収奪効果
    Parameters[3],                                          # 好景気効果
    0.0,                                                    # 退出コスト (0に正規化)
    Parameters[5] + Parameters[4],                          # 参入コスト
])

output_normalized = CCP_to_Value_to_prediction(Normalized_TrueParam, CCP1true, CCP2true, TransitionMat, beta)

diff2 = output_normalized - output
print("Sanity Check 2: 真のパラメータと正規化パラメータの予測CCPの差分")
print("(ゼロに近ければ正規化が正しい)")
print(diff2)
print(f"\\n最大絶対差: {np.abs(diff2).max():.2e}")

print("\\n正規化された真のパラメータ (10x1):")
print(Normalized_TrueParam)
"""))

# ============================================================
# Cell 15: Step 2-3 - P-SD Estimation
# ============================================================
cells.append(nbf.v4.new_markdown_cell(
    "### Step 2-3: P-SD最小二乗法による推定\n"
    "\n"
    "推定されたCCPと遷移確率を用いて、ペイオフパラメータを最小二乗法で推定する。\n"
    "10回のランダム初期値で推定を行い、目的関数が最小のものを選択する。"
))

cells.append(nbf.v4.new_code_cell("""\
# 目的関数のラッパー: 5パラメータ → 10パラメータに変換
def obj(x):
    param10 = np.array([x[0], x[2], x[3], 0.0, x[4],
                        x[1], x[2], x[3], 0.0, x[4]])
    return obj_fun(param10, EstimatedCCP1, EstimatedCCP2, EstimatedTransition, beta)

# 初期値の設定 (真のパラメータに基づく)
# [theta1, theta2, rival, z, entry_cost]
initial = np.array([Parameters[0], Parameters[1], Parameters[2], Parameters[3], Parameters[5]])
print(f"基本初期値: {initial}")

# 10個のランダム初期値 (0.6倍~1.2倍の揺らぎ)
np.random.seed(42)
n_restarts = 10
mat_initial = np.tile(initial.reshape(-1, 1), (1, n_restarts)) * \\
              np.random.uniform(0.6, 1.2, size=(5, n_restarts))

# 各初期値で最適化 (Nelder-Mead)
result = np.zeros((6, n_restarts))  # 5 params + 1 obj value

print(f"\\n{n_restarts}回のランダム初期値で最適化中...")
start_time = time.time()

for i in range(n_restarts):
    sol = minimize(obj, mat_initial[:, i], method='Nelder-Mead')
    result[:5, i] = sol.x
    result[5, i] = sol.fun
    print(f"  初期値 {i+1}: obj = {sol.fun:.6e}, params = {np.round(sol.x, 4)}")

elapsed = time.time() - start_time
print(f"\\n最適化完了: {elapsed:.1f} 秒")

# 最良の結果を選択
best_idx = np.argmin(result[5, :])
result_pick = result[:5, best_idx]

print(f"\\n最良の結果 (初期値 {best_idx+1}):")
print(f"  目的関数値: {result[5, best_idx]:.6e}")
print(f"  推定パラメータ (theta1, theta2, rival, z, entry_cost):")
print(f"  {result_pick}")
"""))

# ============================================================
# Cell 16: Compare with Normalized True
# ============================================================
cells.append(nbf.v4.new_code_cell("""\
# 正規化された真のパラメータ (5次元版)
normalized_trueparam = np.array([
    Normalized_TrueParam[0],  # theta1 (firm1 base profit, normalized)
    Normalized_TrueParam[5],  # theta2 (firm2 base profit, normalized)
    Normalized_TrueParam[1],  # rival effect
    Normalized_TrueParam[2],  # z effect
    Normalized_TrueParam[4],  # entry cost
])

print("Aguirregabiria-Suzuki正規化された真のパラメータ:")
print(normalized_trueparam)

print("\\n推定パラメータ:")
print(result_pick)

print("\\n差分 (推定 - 正規化真):")
print(result_pick - normalized_trueparam)
"""))

# ============================================================
# Cell 17: Bootstrap SE
# ============================================================
cells.append(nbf.v4.new_markdown_cell(
    "## 5. Bootstrapによる標準誤差の計算\n"
    "\n"
    "市場単位でリサンプリングを行い、100回のBootstrapで標準誤差を計算する。"
))

cells.append(nbf.v4.new_code_cell("""\
# Bootstrapの設定
np.random.seed(2023)

numBootSample = 100

# 各Bootstrapサンプルで用いる市場のインデックス (NumSimMarkets x numBootSample)
bootindex = np.random.randint(1, NumSimMarkets + 1, size=(NumSimMarkets, numBootSample))

# 結果格納用
bootresult_transition = np.zeros((2, numBootSample))
bootresult_CCP1 = np.zeros((8, numBootSample))
bootresult_CCP2 = np.zeros((8, numBootSample))
bootresult_payoff = np.zeros((5, numBootSample))

print(f"Bootstrap: {numBootSample} 回のリサンプリング")
start_time = time.time()

for b in range(numBootSample):
    if (b + 1) % 10 == 0 or b == 0:
        print(f"  Bootstrap {b + 1}/{numBootSample}...")

    # Bootstrapサンプルの構築 (市場単位でリサンプリング)
    bootsample_list = []
    for m_idx in range(NumSimMarkets):
        mk = bootindex[m_idx, b]
        temp = FakeData[FakeData[:, 0] == mk]
        bootsample_list.append(temp)
    bootsample = np.vstack(bootsample_list)

    # P-SD推定
    output = Estimation_PS_bootstrap(bootsample, beta)

    # 結果を保存
    bootresult_CCP1[:, b] = output[0][:, 0]
    bootresult_CCP2[:, b] = output[0][:, 1]
    bootresult_transition[:, b] = np.diag(output[1])
    bootresult_payoff[:, b] = output[2]

elapsed = time.time() - start_time
print(f"\\nBootstrap完了: {elapsed:.1f} 秒")
"""))

# ============================================================
# Cell 18: Summary Tables - CCP
# ============================================================
cells.append(nbf.v4.new_markdown_cell("## 結果の表示"))

cells.append(nbf.v4.new_code_cell("""\
# CCP player 1
Summary_CCP1 = pd.DataFrame({
    'CCP for firm 1: True': CCP1UpdatedMat[:, 1],
    'Estimated': EstimatedCCP1,
    'SE': np.std(bootresult_CCP1, axis=1, ddof=0)
}, index=range(1, 9))

print("CCP 企業1:")
print(Summary_CCP1.to_string())

# CCP player 2
Summary_CCP2 = pd.DataFrame({
    'CCP for firm 2: True': CCP2UpdatedMat[:, 1],
    'Estimated': EstimatedCCP2,
    'SE': np.std(bootresult_CCP2, axis=1, ddof=0)
}, index=range(1, 9))

print("\\nCCP 企業2:")
print(Summary_CCP2.to_string())

# Transition Probability
Summary_Transition = pd.DataFrame({
    'Transition Probability (GG and BB): True': [0.7, 0.6],
    'Estimated': np.diag(EstimatedTransition),
    'SE': np.std(bootresult_transition, axis=1, ddof=0)
}, index=range(1, 3))

print("\\n遷移確率:")
print(Summary_Transition.to_string())
"""))

# ============================================================
# Cell 19: Summary Tables - Payoff Parameters
# ============================================================
cells.append(nbf.v4.new_markdown_cell("## 表11.4: P-SD推定の結果"))

cells.append(nbf.v4.new_code_cell("""\
# ペイオフパラメータの比較表
true_param = np.array([Parameters[0], Parameters[1], Parameters[2], Parameters[3], Parameters[5]])
# = [0.3, 0.2, -0.27, 0.45, -2.1]

normalized_param = np.array([
    0.3 - (1 - beta) / beta * (-0.15),    # theta1
    0.2 - (1 - beta) / beta * (-0.15),    # theta2
    -0.27,                                  # rival
    0.45,                                   # z
    -2.1 + (-0.15)                          # entry cost
])

se_payoff = np.std(bootresult_payoff, axis=1, ddof=0)

Summary_Payoff = pd.DataFrame({
    'Payoff parameter: True': true_param,
    'Normalized true': normalized_param,
    'Estimated': result_pick,
    'SE': se_payoff
}, index=range(1, 6))

print("Tab 11.4: P-SDによるペイオフパラメータの推定結果")
print(Summary_Payoff.to_string())

# CSV保存
Summary_Payoff.to_csv(output_dir / 'Tab11_4_PSD.csv')
print("\\n保存しました: output/Tab11_4_PSD.csv")
"""))

# ============================================================
# Cell 20: Summary
# ============================================================
cells.append(nbf.v4.new_markdown_cell("## まとめ"))

cells.append(nbf.v4.new_code_cell("""\
print("=" * 70)
print("第11章 P-SD推定 まとめ")
print("=" * 70)
print()
print("パラメータ名       真の値    正規化真   推定値     SE")
print("-" * 70)
param_names = ['theta1 (firm1)', 'theta2 (firm2)', 'rival effect',
               'z effect', 'entry cost']
for i in range(5):
    print(f"{param_names[i]:20s} {true_param[i]:8.4f}  {normalized_param[i]:8.4f}  "
          f"{result_pick[i]:8.4f}  {se_payoff[i]:8.4f}")
print("-" * 70)
print()
print("注:")
print("- P-SDでは退出コストtheta4が0に正規化されるため、")
print("  推定値は正規化された真のパラメータと比較すべきである。")
print("- 標準誤差は100回のBootstrap (市場単位リサンプリング) により計算。")
print()
print("出力ファイル:")
print("  - output/Tab11_4_PSD.csv")
print()
print("完了")
"""))

nb.cells = cells
nbf.write(nb, 'main_ch11_PSD.ipynb')
print("Generated: main_ch11_PSD.ipynb")
