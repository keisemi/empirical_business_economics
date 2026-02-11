"""Generate main_ch11_policy_sim.ipynb for Chapter 11: Counterfactual Policy Simulation."""
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
    "# 第11章 動学ゲーム: 反実仮想シミュレーション\n"
    "\n"
    "2企業の参入退出の動的ゲームにおいて、パラメータを変化させた\n"
    "反実仮想シミュレーション（Counterfactual Policy Simulation）を行う。\n"
    "\n"
    "- シナリオ1: 企業1のベース利潤が上昇し、顧客奪取効果が両企業とも0になった場合\n"
    "- シナリオ2: 企業1のベース利潤が上昇し、顧客奪取効果に企業間で異質性がある場合"
))

# ============================================================
# Cell 2: Setup
# ============================================================
cells.append(nbf.v4.new_code_cell("""\
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# 日本語フォント設定
import matplotlib
font_candidates = [
    'Hiragino Sans',
    'Hiragino Kaku Gothic ProN',
    'IPAexGothic',
    'Noto Sans CJK JP',
    'Yu Gothic',
    'MS Gothic',
]
available_fonts = [f.name for f in matplotlib.font_manager.fontManager.ttflist]
jp_font = None
for fc in font_candidates:
    if fc in available_fonts:
        jp_font = fc
        break
if jp_font:
    plt.rcParams['font.family'] = jp_font
    print(f"日本語フォント: {jp_font}")
else:
    print("警告: 日本語フォントが見つかりません。文字化けの可能性があります。")
plt.rcParams['axes.unicode_minus'] = False

# パス設定
base_dir = Path('..')
output_dir = base_dir / 'output'
output_dir.mkdir(exist_ok=True)

print("Setup complete.")
"""))

# ============================================================
# Cell 3: Constants
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

# ベースラインパラメータ (10x1)
# theta[0]: 企業1のベース利潤
# theta[1]: ライバルの店舗数が企業1の利潤に与える影響
# theta[2]: 景気が良い時の企業1への追加的利潤
# theta[3]: 企業1の退出のためのコスト
# theta[4]: 企業1の出店のためのコスト
# theta[5]-theta[9]: 企業2について同様
BaselineParameterValues = np.array([
    Parameters[0],   # 企業1のベース利潤
    Parameters[2],   # ライバルの店舗数が企業1の利潤に与える影響
    Parameters[3],   # 景気が良い時の企業1への追加的利潤
    Parameters[4],   # 企業1の退出のためのコスト
    Parameters[5],   # 企業1の出店のためのコスト
    Parameters[1],   # 企業2のベース利潤
    Parameters[2],   # ライバルの店舗数が企業2の利潤に与える影響
    Parameters[3],   # 景気が良い時の企業2への追加的利潤
    Parameters[4],   # 企業2の退出のためのコスト
    Parameters[5],   # 企業2の出店のためのコスト
])

print("パラメータ設定完了")
print(f"beta = {beta}")
print(f"遷移行列:\\n{TransitionMat}")
print(f"BaselineParameterValues = {BaselineParameterValues}")
"""))

# ============================================================
# Cell 4: CCP Adjusters
# ============================================================
cells.append(nbf.v4.new_markdown_cell("## CCP Adjuster の定義"))

cells.append(nbf.v4.new_code_cell("""\
# CCP Adjuster: 各状態で選択可能な行動に1、不可能な行動に0を割り当てる行列 (8x3)
# 列は a_i = -1, 0, 1 に対応
# 状態: G00, G01, G10, G11, B00, B01, B10, B11

# 企業1について:
#   n1=0 のとき: a1=-1 は不可, a1=0 は可, a1=1 は可 → [0,1,1]
#   n1=1 のとき: a1=-1 は可, a1=0 は可, a1=1 は不可 → [1,1,0]
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

print("CCP1Adjuster (8x3):")
print(CCP1Adjuster)
print("\\nCCP2Adjuster (8x3):")
print(CCP2Adjuster)
"""))

# ============================================================
# Cell 5: Helper Functions - Profit and CCP Transform
# ============================================================
cells.append(nbf.v4.new_markdown_cell("## 利潤関数・CCP変換関数の定義"))

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


print("利潤関数・CCP変換関数定義完了")
"""))

# ============================================================
# Cell 6: Helper Functions - Transition Matrices
# ============================================================
cells.append(nbf.v4.new_markdown_cell("## 遷移行列の構築関数"))

cells.append(nbf.v4.new_code_cell("""\
def fP(Matrix1, Vec1, Vec2):
    \"\"\"CCPの下での状態遷移行列 F^{P,sigma} を構築する (8x8)

    Parameters
    ----------
    Matrix1 : (2,2) 景気の遷移行列
    Vec1 : (8,) 企業1のCCPベクトル (P(stay))
    Vec2 : (8,) 企業2のCCPベクトル (P(stay))
    \"\"\"
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

    # a1 = -1: 退出 → n1'=0
    vec_m1 = np.concatenate([np.zeros(16), np.ones(16), np.zeros(16), np.ones(16)])
    MatAdjustMinus1 = vec_m1.reshape(8, 8, order='C')
    vec_m2 = np.concatenate([np.ones(16), np.zeros(16), np.ones(16), np.zeros(16)])
    MatAdjustMinus2 = vec_m2.reshape(8, 8, order='F')
    output1 = TempMat0 * TempMat2 * MatAdjustMinus1 * MatAdjustMinus2

    # a1 = 0: 現状維持
    ForZero = np.array([[1, 0], [0, 1]])
    block = np.kron(ForZero, np.ones((2, 2)))
    MatAdjustZero = np.hstack([block, block])
    MatAdjustZero = np.vstack([MatAdjustZero, MatAdjustZero])
    output2 = TempMat0 * TempMat2 * MatAdjustZero

    # a1 = 1: 参入 → n1'=1
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


print("遷移行列・期待利潤関数定義完了")
"""))

# ============================================================
# Cell 7: f_MPE function
# ============================================================
cells.append(nbf.v4.new_markdown_cell("## MPEソルバーの定義"))

cells.append(nbf.v4.new_code_cell("""\
def f_MPE(TransitionMat, pi1, pi2, beta, tol=1e-12):
    \"\"\"マルコフ完全均衡（MPE）を固定点反復で計算する

    Parameters
    ----------
    TransitionMat : (2,2) 景気の遷移行列
    pi1 : (8,3) 企業1の利潤行列 (CCP Adjuster 適用済み)
    pi2 : (8,3) 企業2の利潤行列 (CCP Adjuster 適用済み)
    beta : float 割引因子
    tol : float 収束判定の閾値

    Returns
    -------
    CCP1UpdatedMat : (8,3) 企業1の均衡CCP行列
    CCP2UpdatedMat : (8,3) 企業2の均衡CCP行列
    ExanteV1 : (8,) 企業1の事前価値関数
    ExanteV2 : (8,) 企業2の事前価値関数
    \"\"\"
    # Step 1: CCPの初期値
    CCP1 = np.full(8, 0.5)
    CCP2 = np.full(8, 0.5)

    CCP1Mat = CCP1Transform(CCP1)
    CCP2Mat = CCP2Transform(CCP2)

    # Step 2: 事前の価値関数を計算
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
# Cell 8: Baseline Equilibrium
# ============================================================
cells.append(nbf.v4.new_markdown_cell("## ベースラインの均衡"))

cells.append(nbf.v4.new_code_cell("""\
# 利潤行列の計算
pi1 = pi1gen(BaselineParameterValues) * CCP1Adjuster
pi2 = pi2gen(BaselineParameterValues) * CCP2Adjuster

print("企業1の利潤行列 pi1 (8x3):")
print(pi1)
print("\\n企業2の利潤行列 pi2 (8x3):")
print(pi2)

# MPEを解く
CCP1base, CCP2base, V1base, V2base = f_MPE(TransitionMat, pi1, pi2, beta)

state_labels = ['G00', 'G01', 'G10', 'G11', 'B00', 'B01', 'B10', 'B11']

print("\\nベースライン均衡CCP:")
eq_df = pd.DataFrame(
    np.hstack([CCP1base, CCP2base]),
    columns=['a1=-1', 'a1=0', 'a1=1', 'a2=-1', 'a2=0', 'a2=1'],
    index=state_labels
)
print(eq_df.to_string())

print("\\nベースライン事前価値関数:")
v_df = pd.DataFrame({
    'State': state_labels,
    'V1base': V1base,
    'V2base': V2base
})
print(v_df.to_string(index=False))
"""))

# ============================================================
# Cell 9: Forward state distribution simulation (baseline)
# ============================================================
cells.append(nbf.v4.new_markdown_cell(
    "## ベースラインの状態分布シミュレーション\n"
    "\n"
    "初期状態を (G, n1=0, n2=0) として、15期間の状態分布の推移を\n"
    "遷移行列を用いて計算する。"
))

cells.append(nbf.v4.new_code_cell("""\
NumSimPeriods = 15

# ベースラインの遷移行列 (CCP1base[:,1] = P(a1=0), CCP2base[:,1] = P(a2=0))
fPsigma_base = fP(TransitionMat, CCP1base[:, 1], CCP2base[:, 1])

# 初期状態: state 0 (G,n1=0,n2=0) に確率1
initial_state = np.zeros(8)
initial_state[0] = 1.0

# 状態分布の推移
transitionpath = np.zeros((NumSimPeriods, 8))
for t in range(NumSimPeriods):
    if t == 0:
        transitionpath[t] = initial_state
    else:
        transitionpath[t] = fPsigma_base.T @ transitionpath[t - 1]

# 企業の店舗存在確率 (n1=1 の状態に対応するマスク)
n1_mask = np.array([0, 0, 1, 1, 0, 0, 1, 1])
n2_mask = np.array([0, 1, 0, 1, 0, 1, 0, 1])

n1 = (transitionpath * n1_mask).sum(axis=1)
n2 = (transitionpath * n2_mask).sum(axis=1)

print("ベースラインの店舗存在確率:")
baseline_df = pd.DataFrame({
    '期': range(1, NumSimPeriods + 1),
    '企業1 (n1)': n1,
    '企業2 (n2)': n2
})
print(baseline_df.to_string(index=False))
"""))

# ============================================================
# Cell 10: Counterfactual Scenario 1
# ============================================================
cells.append(nbf.v4.new_markdown_cell(
    "## 反実仮想シナリオ1\n"
    "\n"
    "企業1のベース利潤が0.5に上昇し、顧客奪取効果がどちらの企業も0になった場合。\n"
    "差別化戦略が成功し、ライバルの存在が利潤に影響しなくなるシナリオ。"
))

cells.append(nbf.v4.new_code_cell("""\
# シナリオ1のパラメータ
CounterfactualParameterValues1 = np.array([
    0.5,             # 企業1のベース利潤（上昇）
    0,               # ライバルの店舗数が企業1の利潤に与える影響（0に）
    Parameters[3],   # 景気が良い時の企業1への追加的利潤
    Parameters[4],   # 企業1の退出のためのコスト
    Parameters[5],   # 企業1の出店のためのコスト
    Parameters[1],   # 企業2のベース利潤
    0,               # ライバルの店舗数が企業2の利潤に与える影響（0に）
    Parameters[3],   # 景気が良い時の企業2への追加的利潤
    Parameters[4],   # 企業2の退出のためのコスト
    Parameters[5],   # 企業2の出店のためのコスト
])

print("シナリオ1パラメータ:", CounterfactualParameterValues1)

# 利潤行列の計算
pi1_cf1 = pi1gen(CounterfactualParameterValues1) * CCP1Adjuster
pi2_cf1 = pi2gen(CounterfactualParameterValues1) * CCP2Adjuster

# MPEを解く
CCP1cf1, CCP2cf1, V1cf1, V2cf1 = f_MPE(TransitionMat, pi1_cf1, pi2_cf1, beta)

# シナリオ1の遷移行列
fPsigma_cf1 = fP(TransitionMat, CCP1cf1[:, 1], CCP2cf1[:, 1])

# 状態分布のシミュレーション
# 注意: R code ではカウンターファクチュアルの遷移に、前期の状態分布として
# ベースラインの transitionpath を使用している
transitionpath_cf1 = np.zeros((NumSimPeriods, 8))
for t in range(NumSimPeriods):
    if t == 0:
        transitionpath_cf1[t] = initial_state
    else:
        transitionpath_cf1[t] = fPsigma_cf1.T @ transitionpath[t - 1]

# 店舗存在確率
n1_cf_s1 = (transitionpath_cf1 * n1_mask).sum(axis=1)
n2_cf_s1 = (transitionpath_cf1 * n2_mask).sum(axis=1)

print("\\nシナリオ1の店舗存在確率:")
cf1_df = pd.DataFrame({
    '期': range(1, NumSimPeriods + 1),
    '企業1 (n1)': n1_cf_s1,
    '企業2 (n2)': n2_cf_s1
})
print(cf1_df.to_string(index=False))
"""))

# ============================================================
# Cell 11: Scenario 1 plot
# ============================================================
cells.append(nbf.v4.new_markdown_cell("### シナリオ1: 店舗存在確率のプロット"))

cells.append(nbf.v4.new_code_cell("""\
periods = range(1, NumSimPeriods + 1)

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# 企業1
axes[0].plot(periods, n1, 'k-', linewidth=1.5, label='ベースライン')
axes[0].plot(periods, n1_cf_s1, 'k--D', markersize=4, linewidth=1, label='シナリオ1')
axes[0].set_xlim(1, 15)
axes[0].set_ylim(0, 0.8)
axes[0].set_xticks(range(1, 16))
axes[0].set_title('企業1の店舗存在確率')
axes[0].set_xlabel('')
axes[0].set_ylabel('')
axes[0].legend()

# 企業2
axes[1].plot(periods, n2, 'k-', linewidth=1.5, label='ベースライン')
axes[1].plot(periods, n2_cf_s1, 'k--D', markersize=4, linewidth=1, label='シナリオ1')
axes[1].set_xlim(1, 15)
axes[1].set_ylim(0, 0.8)
axes[1].set_xticks(range(1, 16))
axes[1].set_title('企業2の店舗存在確率')
axes[1].set_xlabel('')
axes[1].set_ylabel('')
axes[1].legend()

plt.tight_layout()
plt.savefig(output_dir / 'ProbEntry1Original.png', dpi=150, bbox_inches='tight')
plt.show()
print("保存しました: output/ProbEntry1Original.png")
"""))

# ============================================================
# Cell 12: Scenario 1 value comparison
# ============================================================
cells.append(nbf.v4.new_markdown_cell("### シナリオ1: 価値関数の比較"))

cells.append(nbf.v4.new_code_cell("""\
diff_value_1_s1 = V1cf1 - V1base
diff_value_2_s1 = V2cf1 - V2base

print("シナリオ1: 価値関数の比較")
mat_scenario1 = np.column_stack([V1base, V2base, V1cf1, diff_value_1_s1, V2cf1])

scenario1_df = pd.DataFrame(
    mat_scenario1,
    columns=['V1(ベースライン)', 'V2(ベースライン)', 'V1(シナリオ1)', 'V1差分', 'V2(シナリオ1)'],
    index=state_labels
)
print(scenario1_df.to_string())
"""))

# ============================================================
# Cell 13: Counterfactual Scenario 2
# ============================================================
cells.append(nbf.v4.new_markdown_cell(
    "## 反実仮想シナリオ2\n"
    "\n"
    "企業1のベース利潤が0.5に上昇し、顧客奪取効果に企業間で異質性がある場合。\n"
    "企業1の顧客奪取効果は-0.1、企業2は-0.2。"
))

cells.append(nbf.v4.new_code_cell("""\
# シナリオ2のパラメータ
CounterfactualParameterValues2 = np.array([
    0.5,             # 企業1のベース利潤（上昇）
    -0.1,            # ライバルの店舗数が企業1の利潤に与える影響
    Parameters[3],   # 景気が良い時の企業1への追加的利潤
    Parameters[4],   # 企業1の退出のためのコスト
    Parameters[5],   # 企業1の出店のためのコスト
    Parameters[1],   # 企業2のベース利潤
    -0.2,            # ライバルの店舗数が企業2の利潤に与える影響
    Parameters[3],   # 景気が良い時の企業2への追加的利潤
    Parameters[4],   # 企業2の退出のためのコスト
    Parameters[5],   # 企業2の出店のためのコスト
])

print("シナリオ2パラメータ:", CounterfactualParameterValues2)

# 利潤行列の計算
pi1_cf2 = pi1gen(CounterfactualParameterValues2) * CCP1Adjuster
pi2_cf2 = pi2gen(CounterfactualParameterValues2) * CCP2Adjuster

# MPEを解く
CCP1cf2, CCP2cf2, V1cf2, V2cf2 = f_MPE(TransitionMat, pi1_cf2, pi2_cf2, beta)

# シナリオ2の遷移行列
fPsigma_cf2 = fP(TransitionMat, CCP1cf2[:, 1], CCP2cf2[:, 1])

# 状態分布のシミュレーション
# R code と同様にベースラインの transitionpath を使用
transitionpath_cf2 = np.zeros((NumSimPeriods, 8))
for t in range(NumSimPeriods):
    if t == 0:
        transitionpath_cf2[t] = initial_state
    else:
        transitionpath_cf2[t] = fPsigma_cf2.T @ transitionpath[t - 1]

# 店舗存在確率
n1_cf_s2 = (transitionpath_cf2 * n1_mask).sum(axis=1)
n2_cf_s2 = (transitionpath_cf2 * n2_mask).sum(axis=1)

print("\\nシナリオ2の店舗存在確率:")
cf2_df = pd.DataFrame({
    '期': range(1, NumSimPeriods + 1),
    '企業1 (n1)': n1_cf_s2,
    '企業2 (n2)': n2_cf_s2
})
print(cf2_df.to_string(index=False))
"""))

# ============================================================
# Cell 14: Scenario 2 plot
# ============================================================
cells.append(nbf.v4.new_markdown_cell("### シナリオ2: 店舗存在確率のプロット"))

cells.append(nbf.v4.new_code_cell("""\
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# 企業1
axes[0].plot(periods, n1, 'k-', linewidth=1.5, label='ベースライン')
axes[0].plot(periods, n1_cf_s2, 'k--D', markersize=4, linewidth=1, label='シナリオ2')
axes[0].set_xlim(1, 15)
axes[0].set_ylim(0, 0.8)
axes[0].set_xticks(range(1, 16))
axes[0].set_title('企業1の店舗存在確率')
axes[0].set_xlabel('')
axes[0].set_ylabel('')
axes[0].legend()

# 企業2
axes[1].plot(periods, n2, 'k-', linewidth=1.5, label='ベースライン')
axes[1].plot(periods, n2_cf_s2, 'k--D', markersize=4, linewidth=1, label='シナリオ2')
axes[1].set_xlim(1, 15)
axes[1].set_ylim(0, 0.8)
axes[1].set_xticks(range(1, 16))
axes[1].set_title('企業2の店舗存在確率')
axes[1].set_xlabel('')
axes[1].set_ylabel('')
axes[1].legend()

plt.tight_layout()
plt.savefig(output_dir / 'ProbEntry2Original.png', dpi=150, bbox_inches='tight')
plt.show()
print("保存しました: output/ProbEntry2Original.png")
"""))

# ============================================================
# Cell 15: Scenario 2 value comparison
# ============================================================
cells.append(nbf.v4.new_markdown_cell("### シナリオ2: 価値関数の比較"))

cells.append(nbf.v4.new_code_cell("""\
diff_value_1_s2 = V1cf2 - V1base
diff_value_2_s2 = V2cf2 - V2base

print("シナリオ2: 価値関数の比較")
mat_scenario2 = np.column_stack([V1base, V2base, V1cf2, diff_value_1_s2, V2cf2])

scenario2_df = pd.DataFrame(
    mat_scenario2,
    columns=['V1(ベースライン)', 'V2(ベースライン)', 'V1(シナリオ2)', 'V1差分', 'V2(シナリオ2)'],
    index=state_labels
)
print(scenario2_df.to_string())
"""))

# ============================================================
# Cell 16: Combined results table
# ============================================================
cells.append(nbf.v4.new_markdown_cell(
    "## 表11.6: 反実仮想シミュレーション結果のまとめ\n"
    "\n"
    "シナリオ1・2の結果をまとめてCSVに保存する。"
))

cells.append(nbf.v4.new_code_cell("""\
# mat_all = cbind(mat_scenario1, mat_scenario2[:, 2:4])
# mat_scenario1: [V1base, V2base, V1cf1, diff1_s1, V2cf1]  (5列)
# mat_scenario2[:, 2:4]: [V1cf2, diff1_s2, V2cf2]           (3列)
mat_all = np.hstack([mat_scenario1, mat_scenario2[:, 2:5]])

# CSV保存 (R版と同じフォーマット)
mat_all_df = pd.DataFrame(
    mat_all,
    columns=['V1', 'V2', 'V3', 'V4', 'V5', 'V6', 'V7', 'V8']
)
mat_all_df.index = range(1, 9)
mat_all_df.to_csv(output_dir / 'Tab11_6_CF_simulation.csv')
print("保存しました: output/Tab11_6_CF_simulation.csv")

# 見やすい形式で表示
print("\\n表11.6: 反実仮想シミュレーション結果")
display_df = pd.DataFrame(
    mat_all,
    columns=[
        'V1(基準)', 'V2(基準)',
        'V1(S1)', 'V1差分(S1)', 'V2(S1)',
        'V1(S2)', 'V1差分(S2)', 'V2(S2)'
    ],
    index=state_labels
)
print(display_df.to_string())

# R出力との比較
try:
    r_tab = pd.read_csv(output_dir / 'Tab11_6_CF_simulation.csv', index_col=0)
    r_values = r_tab.values
    diff = np.abs(mat_all - r_values)
    print(f"\\nR出力との最大差: {diff.max():.2e}")
except Exception as e:
    print(f"\\nR出力との比較をスキップ: {e}")
"""))

# ============================================================
# Cell 17: Combined plot
# ============================================================
cells.append(nbf.v4.new_markdown_cell("## 結果のまとめプロット"))

cells.append(nbf.v4.new_code_cell("""\
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# 企業1
axes[0].plot(periods, n1, 'k-', linewidth=1.5, label='ベースライン')
axes[0].plot(periods, n1_cf_s1, 'k--o', markersize=5, markerfacecolor='none',
             linewidth=1, label='シナリオ1')
axes[0].plot(periods, n1_cf_s2, 'k--o', markersize=5, markerfacecolor='black',
             linewidth=1, label='シナリオ2')
axes[0].set_xlim(1, 15)
axes[0].set_ylim(0, 0.8)
axes[0].set_xticks(range(1, 16))
axes[0].set_title('企業1の店舗存在確率')
axes[0].set_xlabel('')
axes[0].set_ylabel('')
axes[0].legend()

# 企業2
axes[1].plot(periods, n2, 'k-', linewidth=1.5, label='ベースライン')
axes[1].plot(periods, n2_cf_s1, 'k--o', markersize=5, markerfacecolor='none',
             linewidth=1, label='シナリオ1')
axes[1].plot(periods, n2_cf_s2, 'k--o', markersize=5, markerfacecolor='black',
             linewidth=1, label='シナリオ2')
axes[1].set_xlim(1, 15)
axes[1].set_ylim(0, 0.8)
axes[1].set_xticks(range(1, 16))
axes[1].set_title('企業2の店舗存在確率')
axes[1].set_xlabel('')
axes[1].set_ylabel('')
axes[1].legend()

plt.tight_layout()
plt.savefig(output_dir / 'ProbEntry3Plots.png', dpi=150, bbox_inches='tight')
plt.show()
print("保存しました: output/ProbEntry3Plots.png")
"""))

# ============================================================
# Cell 18: Summary
# ============================================================
cells.append(nbf.v4.new_markdown_cell("## まとめ"))

cells.append(nbf.v4.new_code_cell("""\
print("=" * 60)
print("第11章: 動学ゲーム 反実仮想シミュレーション まとめ")
print("=" * 60)
print(f"状態空間: 8 状態 (景気{{G,B}} x n1{{0,1}} x n2{{0,1}})")
print(f"各企業の行動: 3種類 (退出=-1, 現状維持=0, 参入=1)")
print(f"割引因子 beta = {beta}")
print(f"遷移行列: P(G|G)={TransitionMat[0,0]}, P(B|B)={TransitionMat[1,1]}")
print()
print("シナリオ1: 企業1のベース利潤=0.5, 顧客奪取効果=0 (両企業)")
print(f"  ベースラインパラメータ: {BaselineParameterValues}")
print(f"  シナリオ1パラメータ:    {CounterfactualParameterValues1}")
print()
print("シナリオ2: 企業1のベース利潤=0.5, 顧客奪取効果: 企業1=-0.1, 企業2=-0.2")
print(f"  シナリオ2パラメータ:    {CounterfactualParameterValues2}")
print()
print("出力ファイル:")
print("  - ProbEntry1Original.png (シナリオ1 店舗存在確率)")
print("  - ProbEntry2Original.png (シナリオ2 店舗存在確率)")
print("  - ProbEntry3Plots.png    (まとめプロット)")
print("  - Tab11_6_CF_simulation.csv (表11.6)")
print()
print("完了")
"""))

nb.cells = cells
nbf.write(nb, 'main_ch11_policy_sim.ipynb')
print("Generated: main_ch11_policy_sim.ipynb")
