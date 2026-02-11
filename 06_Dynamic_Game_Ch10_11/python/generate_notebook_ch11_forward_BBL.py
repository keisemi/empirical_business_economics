"""Generate main_ch11_forward_BBL.ipynb for Chapter 11:
Forward Simulation-based P-SD estimation and BBL inequality estimator."""
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
    "# 第11章 動的ゲーム：Forward Simulation を用いた P-SD 推定および BBL 不等式推定量\n"
    "\n"
    "本ノートブックでは以下を実行する。\n"
    "\n"
    "1. 均衡 CCP の計算（MPE）\n"
    "2. 疑似データの読み込みと CCP/遷移確率の推定\n"
    "3. Forward Simulation を用いた P-SD 推定 + Bootstrap\n"
    "4. BBL 不等式推定量 + Bootstrap\n"
    "\n"
    "R ソース: `main4_Estimation_Forward_PSD_BBL.R`, `sub_5_Bootstrap_PSD_forward.R`"
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

# Matlab データを使用するか
isUseMatlabData = True

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
# Parameters[0]: 企業1のベース利潤, [1]: 企業2のベース利潤
# [2]: 顧客収奪効果, [3]: 景気が良い時の追加的利潤
# [4]: 退出のためのコスト, [5]: 参入のためのコスト
Parameters = np.array([0.3, 0.2, -0.27, 0.45, -0.15, -2.10])

# 10次元パラメータベクトルに並べ替え
TrueParameterValues = np.array([
    Parameters[0],   # theta[0]: 企業1のベース利潤
    Parameters[2],   # theta[1]: ライバル効果
    Parameters[3],   # theta[2]: 好景気効果
    Parameters[4],   # theta[3]: 退出コスト
    Parameters[5],   # theta[4]: 参入コスト
    Parameters[1],   # theta[5]: 企業2のベース利潤
    Parameters[2],   # theta[6]: ライバル効果
    Parameters[3],   # theta[7]: 好景気効果
    Parameters[4],   # theta[8]: 退出コスト
    Parameters[5],   # theta[9]: 参入コスト
])

# CCP Adjuster: 各状態で選択可能な行動に1、不可能な行動に0
# 列: a_i = -1, 0, 1
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

print("パラメータ設定完了")
print(f"beta = {beta}")
print(f"遷移行列:\\n{TransitionMat}")
print(f"TrueParameterValues = {TrueParameterValues}")
"""))

# ============================================================
# Cell 4: Helper Functions - Profit
# ============================================================
cells.append(nbf.v4.new_markdown_cell("## 利潤関数と CCP 変換関数"))

cells.append(nbf.v4.new_code_cell("""\
def pi1gen(theta):
    \"\"\"企業1の利潤行列 (8x3). 行:状態, 列:行動(a1=-1,0,1)\"\"\"
    base = np.array([
        0,                              # G00
        0,                              # G01
        theta[0] + theta[2],            # G10
        theta[0] + theta[1] + theta[2], # G11
        0,                              # B00
        0,                              # B01
        theta[0],                       # B10
        theta[0] + theta[1],            # B11
    ])
    invdiv = np.array([theta[3], 0, theta[4]])
    return np.tile(base.reshape(-1, 1), (1, 3)) + np.tile(invdiv.reshape(1, -1), (8, 1))


def pi2gen(theta):
    \"\"\"企業2の利潤行列 (8x3)\"\"\"
    base = np.array([
        0,                              # G00
        theta[5] + theta[7],            # G01
        0,                              # G10
        theta[5] + theta[6] + theta[7], # G11
        0,                              # B00
        theta[5],                       # B01
        0,                              # B10
        theta[5] + theta[6],            # B11
    ])
    invdiv = np.array([theta[8], 0, theta[9]])
    return np.tile(base.reshape(-1, 1), (1, 3)) + np.tile(invdiv.reshape(1, -1), (8, 1))


def CCP1Transform(x):
    \"\"\"CCP1ベクトル(8) -> CCP1行列(8x3). x[s]=P(a1=0|s).\"\"\"
    return np.array([
        [0,      x[0], 1-x[0]],
        [0,      x[1], 1-x[1]],
        [1-x[2], x[2], 0     ],
        [1-x[3], x[3], 0     ],
        [0,      x[4], 1-x[4]],
        [0,      x[5], 1-x[5]],
        [1-x[6], x[6], 0     ],
        [1-x[7], x[7], 0     ],
    ])


def CCP2Transform(x):
    \"\"\"CCP2ベクトル(8) -> CCP2行列(8x3). x[s]=P(a2=0|s).\"\"\"
    return np.array([
        [0,      x[0], 1-x[0]],
        [1-x[1], x[1], 0     ],
        [0,      x[2], 1-x[2]],
        [1-x[3], x[3], 0     ],
        [0,      x[4], 1-x[4]],
        [1-x[5], x[5], 0     ],
        [0,      x[6], 1-x[6]],
        [1-x[7], x[7], 0     ],
    ])


def CCP1LogTransform(x):
    \"\"\"CCP1ベクトル(8) -> log CCP1行列(8x3)\"\"\"
    return np.array([
        [0,           np.log(x[0]), np.log(1-x[0])],
        [0,           np.log(x[1]), np.log(1-x[1])],
        [np.log(1-x[2]), np.log(x[2]), 0           ],
        [np.log(1-x[3]), np.log(x[3]), 0           ],
        [0,           np.log(x[4]), np.log(1-x[4])],
        [0,           np.log(x[5]), np.log(1-x[5])],
        [np.log(1-x[6]), np.log(x[6]), 0           ],
        [np.log(1-x[7]), np.log(x[7]), 0           ],
    ])


def CCP2LogTransform(x):
    \"\"\"CCP2ベクトル(8) -> log CCP2行列(8x3)\"\"\"
    return np.array([
        [0,           np.log(x[0]), np.log(1-x[0])],
        [np.log(1-x[1]), np.log(x[1]), 0           ],
        [0,           np.log(x[2]), np.log(1-x[2])],
        [np.log(1-x[3]), np.log(x[3]), 0           ],
        [0,           np.log(x[4]), np.log(1-x[4])],
        [np.log(1-x[5]), np.log(x[5]), 0           ],
        [0,           np.log(x[6]), np.log(1-x[6])],
        [np.log(1-x[7]), np.log(x[7]), 0           ],
    ])


# 真のパラメータで利潤行列を計算
pi1 = pi1gen(TrueParameterValues) * CCP1Adjuster
pi2 = pi2gen(TrueParameterValues) * CCP2Adjuster

print("利潤関数・CCP変換関数 定義完了")
print("企業1の利潤行列 pi1:")
print(pi1)
print("\\n企業2の利潤行列 pi2:")
print(pi2)
"""))

# ============================================================
# Cell 5: Transition functions
# ============================================================
cells.append(nbf.v4.new_markdown_cell("## 遷移行列関数"))

cells.append(nbf.v4.new_code_cell("""\
def fP(Matrix1, Vec1, Vec2):
    \"\"\"CCPの下での状態遷移行列 (8x8)\"\"\"
    TempMat0 = np.kron(Matrix1, np.ones((4, 4)))

    rows1 = np.array([
        [Vec1[0], 1-Vec1[0]], [Vec1[1], 1-Vec1[1]],
        [1-Vec1[2], Vec1[2]], [1-Vec1[3], Vec1[3]],
        [Vec1[4], 1-Vec1[4]], [Vec1[5], 1-Vec1[5]],
        [1-Vec1[6], Vec1[6]], [1-Vec1[7], Vec1[7]],
    ])
    TempMat1 = np.kron(rows1, np.ones((1, 2)))
    TempMat1 = np.hstack([TempMat1, TempMat1])

    rows2 = np.array([
        [Vec2[0], 1-Vec2[0]], [1-Vec2[1], Vec2[1]],
        [Vec2[2], 1-Vec2[2]], [1-Vec2[3], Vec2[3]],
        [Vec2[4], 1-Vec2[4]], [1-Vec2[5], Vec2[5]],
        [Vec2[6], 1-Vec2[6]], [1-Vec2[7], Vec2[7]],
    ])
    TempMat2 = np.kron(np.ones((1, 4)), rows2)

    return TempMat0 * TempMat1 * TempMat2


def fP_a1given(Matrix1, Vec2):
    \"\"\"企業1の行動を条件付けた遷移行列リスト [a1=-1, a1=0, a1=1] (各8x8)\"\"\"
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
    \"\"\"企業2の行動を条件付けた遷移行列リスト [a2=-1, a2=0, a2=1] (各8x8)\"\"\"
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


print("遷移行列関数 定義完了")
"""))

# ============================================================
# Cell 6: piPsigmaGen and f_MPE
# ============================================================
cells.append(nbf.v4.new_markdown_cell("## 期待利潤関数と MPE ソルバー"))

cells.append(nbf.v4.new_code_cell("""\
def pi1PsigmaGen(pi1, Mat2):
    \"\"\"企業2のCCP行列の下での企業1の期待利潤 (8x3)\"\"\"
    ones3 = np.ones(3)
    pi1_dec = (pi1[:, 0].reshape(-1, 1) * Mat2) @ ones3
    pi1_0   = (pi1[:, 1].reshape(-1, 1) * Mat2) @ ones3
    pi1_inc = (pi1[:, 2].reshape(-1, 1) * Mat2) @ ones3
    return np.column_stack([pi1_dec, pi1_0, pi1_inc])


def pi2PsigmaGen(pi2, Mat1):
    \"\"\"企業1のCCP行列の下での企業2の期待利潤 (8x3)\"\"\"
    ones3 = np.ones(3)
    pi2_dec = (pi2[:, 0].reshape(-1, 1) * Mat1) @ ones3
    pi2_0   = (pi2[:, 1].reshape(-1, 1) * Mat1) @ ones3
    pi2_inc = (pi2[:, 2].reshape(-1, 1) * Mat1) @ ones3
    return np.column_stack([pi2_dec, pi2_0, pi2_inc])


def f_MPE(TransitionMat, pi1, pi2, beta, tol=1e-12):
    \"\"\"MPEを固定点反復で計算する.
    Returns: (CCP1UpdatedMat, CCP2UpdatedMat, ExanteV1, ExanteV2)
    \"\"\"
    CCP1 = np.full(8, 0.5)
    CCP2 = np.full(8, 0.5)

    CCP1Mat = CCP1Transform(CCP1)
    CCP2Mat = CCP2Transform(CCP2)

    # Step 2: 事前価値関数
    fPsigma = fP(TransitionMat, CCP1, CCP2)
    pi1Psigma = pi1PsigmaGen(pi1, CCP2Mat)
    pi2Psigma = pi2PsigmaGen(pi2, CCP1Mat)
    eP1 = eulergamma - CCP1LogTransform(CCP1)
    eP2 = eulergamma - CCP2LogTransform(CCP2)

    inv_mat = np.linalg.inv(np.eye(8) - beta * fPsigma)
    ExanteV1 = inv_mat @ np.sum(CCP1Mat * (pi1Psigma + eP1), axis=1)
    ExanteV2 = inv_mat @ np.sum(CCP2Mat * (pi2Psigma + eP2), axis=1)

    # Step 3: CCP更新
    fP_a1 = fP_a1given(TransitionMat, CCP2)
    fP_a2 = fP_a2given(TransitionMat, CCP1)

    future1 = np.column_stack([fP_a1[k] @ ExanteV1 for k in range(3)])
    NewSigmaSeed1 = (pi1Psigma + beta * future1) * CCP1Adjuster
    NewSigmaDeno1 = np.sum(np.exp(NewSigmaSeed1), axis=1) - 1.0
    CCP1UpdatedMat = np.exp(NewSigmaSeed1) / NewSigmaDeno1.reshape(-1, 1) * CCP1Adjuster
    CCP1Updated = CCP1UpdatedMat[:, 1]

    future2 = np.column_stack([fP_a2[k] @ ExanteV2 for k in range(3)])
    NewSigmaSeed2 = (pi2Psigma + beta * future2) * CCP2Adjuster
    NewSigmaDeno2 = np.sum(np.exp(NewSigmaSeed2), axis=1) - 1.0
    CCP2UpdatedMat = np.exp(NewSigmaSeed2) / NewSigmaDeno2.reshape(-1, 1) * CCP2Adjuster
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

    # Step 5: 収束まで反復
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
        CCP1UpdatedMat = np.exp(NewSigmaSeed1) / NewSigmaDeno1.reshape(-1, 1) * CCP1Adjuster
        CCP1Updated = CCP1UpdatedMat[:, 1]

        future2 = np.column_stack([fP_a2[k] @ ExanteV2 for k in range(3)])
        NewSigmaSeed2 = (pi2Psigma + beta * future2) * CCP2Adjuster
        NewSigmaDeno2 = np.sum(np.exp(NewSigmaSeed2), axis=1) - 1.0
        CCP2UpdatedMat = np.exp(NewSigmaSeed2) / NewSigmaDeno2.reshape(-1, 1) * CCP2Adjuster
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

    print(f"MPE 収束: {iteration} 回の反復, 差分 = {DiffExanteV:.2e}")
    return CCP1UpdatedMat, CCP2UpdatedMat, ExanteV1Updated, ExanteV2Updated


print("期待利潤関数・MPEソルバー 定義完了")
"""))

# ============================================================
# Cell 7: Compute Equilibrium
# ============================================================
cells.append(nbf.v4.new_markdown_cell("## 均衡 CCP の計算"))

cells.append(nbf.v4.new_code_cell("""\
start_time = time.time()
CCP1UpdatedMat, CCP2UpdatedMat, ExanteV1, ExanteV2 = f_MPE(TransitionMat, pi1, pi2, beta)
elapsed = time.time() - start_time
print(f"計算時間: {elapsed:.3f} 秒")

print("\\n均衡 CCP (企業1):")
print(CCP1UpdatedMat)
print("\\n均衡 CCP (企業2):")
print(CCP2UpdatedMat)
print("\\nExanteV1:", ExanteV1)
print("ExanteV2:", ExanteV2)
"""))

# ============================================================
# Cell 8: Data Loading & CCP Estimation
# ============================================================
cells.append(nbf.v4.new_markdown_cell(
    "## データの読み込みと CCP/遷移確率の推定\n"
    "\n"
    "Matlab で生成された FakeData (500 市場 x 50 期間) を読み込み、\n"
    "CCP と遷移確率を推定する。"
))

cells.append(nbf.v4.new_code_cell("""\
# FakeData の読み込み
if isUseMatlabData:
    FakeData = pd.read_csv(data_dir / 'FakeData_Matlab.csv', header=None).values
    print(f"FakeData_Matlab.csv 読み込み: shape = {FakeData.shape}")
else:
    raise NotImplementedError("Python生成データは別途 ch10 ノートブックで生成してください")

# 列の意味: [Market_ID, Time, State(1-8), Demand, n1, n2, a1, a2]

# --- CCP の推定 ---
EstimatedCCP1 = np.zeros(8)
EstimatedCCP2 = np.zeros(8)

for s in range(1, 9):  # 1-indexed
    sub = FakeData[FakeData[:, 2] == s]
    EstimatedCCP1[s-1] = np.sum(sub[:, 6] == 0) / len(sub)
    EstimatedCCP2[s-1] = np.sum(sub[:, 7] == 0) / len(sub)

print("推定 CCP1 (P(a1=0)):", np.round(EstimatedCCP1, 6))
print("推定 CCP2 (P(a2=0)):", np.round(EstimatedCCP2, 6))

# --- 遷移確率の推定 ---
EstimatedTransition = np.zeros((2, 2))

# 1期ラグの景気状態
n_rows = len(FakeData)
lag_demand = np.zeros(n_rows)
lag_demand[1:] = FakeData[:-1, 3]

# t != 1 のデータのみ使用
mask = FakeData[:, 1] != 1
data_with_lag = np.column_stack([FakeData[mask], lag_demand[mask]])

for z in range(1, 3):
    sub_z = data_with_lag[data_with_lag[:, 3] == z]
    sub_zz = sub_z[sub_z[:, 8] == z]
    EstimatedTransition[z-1, z-1] = len(sub_zz) / len(sub_z)
    EstimatedTransition[z-1, 2-z] = 1 - EstimatedTransition[z-1, z-1]

print("\\n推定遷移確率行列:")
print(EstimatedTransition)
"""))

# ============================================================
# Cell 9: VSigmaGeneration
# ============================================================
cells.append(nbf.v4.new_markdown_cell(
    "## Forward Simulation 関数 (VSigmaGeneration)\n"
    "\n"
    "CCP と遷移確率の下で Forward Simulation を行い、\n"
    "Value function の基底 W (6 x 8) を生成する。\n"
    "6つの基底: [base_profit, rival_effect, boom_effect, exit/entry_action, 0(unused), random_profit]"
))

cells.append(nbf.v4.new_code_cell("""\
def VSigmaGeneration(CCP1, CCP2, EstimatedTransition, EVrandom, UNIrandom,
                     InitialState, NumSimMarkets, NumSimulations, NumSimPeriods, beta):
    \"\"\"Forward Simulation による Value function の基底生成 (ベクトル化版)

    Parameters
    ----------
    CCP1, CCP2 : (8,) CCP ベクトル (P(stay))
    EVrandom : (NumSimMarkets, NumSimPeriods, 2, NumSimulations, 8, 3) Gumbel ショック
    UNIrandom : (NumSimMarkets, NumSimPeriods, NumSimulations) 景気遷移用一様乱数
    InitialState : (NumSimMarkets, 2) 初期状態 (2列目を使用, 1-indexed)

    Returns
    -------
    W1out : (6, NumSimMarkets) 企業1の基底 (シミュレーション平均)
    W2out : (6, NumSimMarkets) 企業2の基底
    \"\"\"
    W1 = np.zeros((6, NumSimMarkets, NumSimulations))
    W2 = np.zeros((6, NumSimMarkets, NumSimulations))

    # 閾値: log(CCP) - log(1-CCP)
    ThresholdValue1 = np.log(CCP1) - np.log(1 - CCP1)
    ThresholdValue2 = np.log(CCP2) - np.log(1 - CCP2)

    # State encoding: state = 4*boom + 2*n1 + n2
    # n1 presence: states 2,3,6,7 -> n1=1; states 0,1,4,5 -> n1=0
    n1_present = np.array([0, 0, 1, 1, 0, 0, 1, 1])
    n2_present = np.array([0, 1, 0, 1, 0, 1, 0, 1])
    is_boom = np.array([0, 0, 0, 0, 1, 1, 1, 1])

    # W1Seed base_profit lookup tables (for states 0-7)
    # base_profit_indicator for firm1: 1 if n1=1
    w1_base = n1_present.astype(float)
    # rival_present for firm1: 1 if n2=1 AND n1=1 (states 3,7)
    w1_rival = np.array([0, 0, 0, 1, 0, 0, 0, 1], dtype=float)
    # boom_indicator for firm1: (1-is_boom) if n1=1 else 0
    w1_boom = np.array([0, 0, 1, 1, 0, 0, 0, 0], dtype=float)
    # W2Seed lookup tables
    w2_base = n2_present.astype(float)
    w2_rival = np.array([0, 0, 0, 1, 0, 0, 0, 1], dtype=float)
    w2_boom = np.array([0, 1, 0, 1, 0, 0, 0, 0], dtype=float)

    for mrkt in range(NumSimMarkets):
        # All simulations start from same state (vectorized over sim dim)
        States = np.full(NumSimulations, int(InitialState[mrkt, 1]) - 1, dtype=int)

        for t in range(NumSimPeriods):
            bt = beta ** t
            sim_idx = np.arange(NumSimulations)

            # Get EV shocks for current states
            # EVrandom[mrkt, t, firm, sim, state, action]
            ev1 = EVrandom[mrkt, t, 0]  # (NumSim, 8, 3)
            ev2 = EVrandom[mrkt, t, 1]  # (NumSim, 8, 3)

            # Get shocks for current states
            ev1_s = ev1[sim_idx, States, :]  # (NumSim, 3)
            ev2_s = ev2[sim_idx, States, :]  # (NumSim, 3)

            # Thresholds for current states
            tv1 = ThresholdValue1[States]  # (NumSim,)
            tv2 = ThresholdValue2[States]  # (NumSim,)

            # Firm 1 actions
            is_n1_zero = (n1_present[States] == 0)  # (NumSim,)
            # If n1=0: entry decision (action 2 vs 1)
            DiffFirm1_entry = tv1 - (-ev1_s[:, 1] + ev1_s[:, 2])
            # If n1=1: exit decision (action 0 vs 1)
            DiffFirm1_exit = tv1 - (-ev1_s[:, 1] + ev1_s[:, 0])

            a1 = np.where(is_n1_zero,
                         (DiffFirm1_entry < 0).astype(int),
                         (DiffFirm1_exit < 0).astype(int))

            e1 = np.where(is_n1_zero,
                         np.where(a1 == 1, ev1_s[:, 2], ev1_s[:, 1]),
                         np.where(a1 == 1, ev1_s[:, 0], ev1_s[:, 1]))

            # Firm 2 actions
            is_n2_zero = (n2_present[States] == 0)
            DiffFirm2_entry = tv2 - (-ev2_s[:, 1] + ev2_s[:, 2])
            DiffFirm2_exit = tv2 - (-ev2_s[:, 1] + ev2_s[:, 0])

            a2 = np.where(is_n2_zero,
                         (DiffFirm2_entry < 0).astype(int),
                         (DiffFirm2_exit < 0).astype(int))

            e2 = np.where(is_n2_zero,
                         np.where(a2 == 1, ev2_s[:, 2], ev2_s[:, 1]),
                         np.where(a2 == 1, ev2_s[:, 0], ev2_s[:, 1]))

            # W seeds from lookup tables
            W1[:, mrkt, :] += bt * np.array([
                w1_base[States],
                w1_rival[States],
                w1_boom[States],
                np.where(n1_present[States] == 1, a1.astype(float), 0.0),
                np.where(n1_present[States] == 0, a1.astype(float), 0.0),
                e1
            ])

            W2[:, mrkt, :] += bt * np.array([
                w2_base[States],
                w2_rival[States],
                w2_boom[States],
                np.where(n2_present[States] == 1, a2.astype(float), 0.0),
                np.where(n2_present[States] == 0, a2.astype(float), 0.0),
                e2
            ])

            # State transition
            # Next n1: if currently n1=0, new_n1=a1; if n1=1, new_n1=(1-a1)=stay
            new_n1 = np.where(is_n1_zero, a1, 1 - a1)
            new_n2 = np.where(is_n2_zero, a2, 1 - a2)

            # Exogenous state transition
            cur_boom = is_boom[States]
            uni = UNIrandom[mrkt, t, :]  # (NumSim,)
            # If good (boom=0): transition to bad with prob 1-trans[0,0]
            # If bad (boom=1): stay bad with prob trans[1,1]
            new_boom = np.where(cur_boom == 0,
                               (uni > EstimatedTransition[0, 0]).astype(int),
                               (uni < EstimatedTransition[1, 1]).astype(int))

            States = 4 * new_boom + 2 * new_n1 + new_n2

    W1out = np.mean(W1, axis=2)  # (6, NumSimMarkets)
    W2out = np.mean(W2, axis=2)

    return W1out, W2out


print("VSigmaGeneration 関数 定義完了")
"""))

# ============================================================
# Cell 10: Forward simulation setup and random draws
# ============================================================
cells.append(nbf.v4.new_markdown_cell(
    "## Forward Simulation のセットアップ\n"
    "\n"
    "P-SD 推定用の Forward Simulation パラメータと乱数を準備する。"
))

cells.append(nbf.v4.new_code_cell("""\
# Forward Simulation パラメータ
NumSimPeriods = 100
NumSimFirms = 2
NumSimulations = 1000

# Initial State: 8通りの初期状態
InitialState = np.column_stack([np.arange(1, 9), np.arange(1, 9)])  # (8, 2)
NumSimMarkets = 8

# 乱数の生成
np.random.seed(2023)

# Gumbel ショック: F(x) = exp(-exp(-x)), 逆関数法
total_ev = NumSimMarkets * NumSimPeriods * NumSimFirms * NumSimulations * 8 * 3
EVrandom = -np.log(-np.log(np.random.uniform(size=total_ev)))
EVrandom = EVrandom.reshape((NumSimMarkets, NumSimPeriods, NumSimFirms, NumSimulations, 8, 3), order='F')

# 景気遷移用一様乱数
total_uni = NumSimMarkets * NumSimPeriods * NumSimulations
UNIrandom = np.random.uniform(size=total_uni)
UNIrandom = UNIrandom.reshape((NumSimMarkets, NumSimPeriods, NumSimulations), order='F')

print(f"EVrandom shape: {EVrandom.shape}")
print(f"UNIrandom shape: {UNIrandom.shape}")
print(f"NumSimPeriods={NumSimPeriods}, NumSimulations={NumSimulations}, NumSimMarkets={NumSimMarkets}")
"""))

# ============================================================
# Cell 11: Sanity check with true parameters
# ============================================================
cells.append(nbf.v4.new_markdown_cell(
    "## Sanity Check: 真のパラメータでの Forward Simulation\n"
    "\n"
    "DGP における ExanteV と Forward Simulation した Value が一致するか確認する。"
))

cells.append(nbf.v4.new_code_cell("""\
print("Sanity check: 真のCCPでForward Simulation...")
start_time = time.time()
W1star_true, W2star_true = VSigmaGeneration(
    CCP1UpdatedMat[:, 1], CCP2UpdatedMat[:, 1], TransitionMat, EVrandom,
    UNIrandom, InitialState, NumSimMarkets, NumSimulations, NumSimPeriods, beta
)
elapsed = time.time() - start_time
print(f"計算時間: {elapsed:.1f} 秒")

# W1star_true は (6, 8), W2star_true は (6, 8)
param1 = np.append(TrueParameterValues[:5], 1.0)
param2 = np.append(TrueParameterValues[5:10], 1.0)

V1_sim = W1star_true.T @ param1
V2_sim = W2star_true.T @ param2

print("\\n企業1: ExanteV (true), Forward Sim, 差分")
for i in range(8):
    print(f"  State {i+1}: {ExanteV1[i]:10.6f}  {V1_sim[i]:10.6f}  {ExanteV1[i] - V1_sim[i]:10.6f}")

print("\\n企業2: ExanteV (true), Forward Sim, 差分")
for i in range(8):
    print(f"  State {i+1}: {ExanteV2[i]:10.6f}  {V2_sim[i]:10.6f}  {ExanteV2[i] - V2_sim[i]:10.6f}")

# 正規化後の確認
Normalized_TrueParam = np.array([
    Parameters[0] - (1-beta)/beta * Parameters[4],
    Parameters[2],
    Parameters[3],
    0,
    Parameters[5] + Parameters[4],
    Parameters[1] - (1-beta)/beta * Parameters[4],
    Parameters[2],
    Parameters[3],
    0,
    Parameters[5] + Parameters[4],
])
normparam1 = np.append(Normalized_TrueParam[:5], 1.0)
normparam2 = np.append(Normalized_TrueParam[5:10], 1.0)

V1_norm = W1star_true.T @ normparam1
V2_norm = W2star_true.T @ normparam2

print("\\n--- 正規化パラメータでの比較 ---")
print("企業1: ExanteV, Simulated(original), Simulated(normalized)")
for i in range(8):
    print(f"  State {i+1}: {ExanteV1[i]:10.6f}  {V1_sim[i]:10.6f}  {V1_norm[i]:10.6f}")
print("\\n差分 (original - normalized):")
print(V1_sim - V1_norm)
"""))

# ============================================================
# Cell 12: Forward Simulation with Estimated CCPs
# ============================================================
cells.append(nbf.v4.new_markdown_cell("## 推定 CCP の下での Forward Simulation"))

cells.append(nbf.v4.new_code_cell("""\
print("推定CCPでForward Simulation...")
start_time = time.time()
W1star, W2star = VSigmaGeneration(
    EstimatedCCP1, EstimatedCCP2, EstimatedTransition, EVrandom,
    UNIrandom, InitialState, NumSimMarkets, NumSimulations, NumSimPeriods, beta
)
elapsed = time.time() - start_time
print(f"計算時間: {elapsed:.1f} 秒")
print(f"W1star shape: {W1star.shape}")
print(f"W2star shape: {W2star.shape}")
"""))

# ============================================================
# Cell 13: Estimation_forward_PSD function
# ============================================================
cells.append(nbf.v4.new_markdown_cell(
    "## P-SD 推定関数の定義\n"
    "\n"
    "Forward Simulation した Value を用いた P-SD の目的関数を定義する。"
))

cells.append(nbf.v4.new_code_cell("""\
def Estimation_forward_PSD(param10, W1star, W2star, TransitionMat, CCP1, CCP2, beta):
    \"\"\"Forward Simulation を用いた P-SD の目的関数 (SSE)

    Parameters
    ----------
    param10 : (10,) パラメータベクトル
    W1star, W2star : (6, 8) Forward Simulation の基底
    CCP1, CCP2 : (8,) 推定された CCP
    Returns : float (目的関数値)
    \"\"\"
    param1 = np.append(param10[:5], 1.0)
    param2 = np.append(param10[5:10], 1.0)

    ExanteV1 = W1star.T @ param1  # (8,)
    ExanteV2 = W2star.T @ param2  # (8,)

    # 利潤行列
    pi1_local = pi1gen(param10) * CCP1Adjuster
    pi2_local = pi2gen(param10) * CCP2Adjuster

    CCP1Mat = CCP1Transform(CCP1)
    CCP2Mat = CCP2Transform(CCP2)

    pi1Psigma = pi1PsigmaGen(pi1_local, CCP2Mat)
    pi2Psigma = pi2PsigmaGen(pi2_local, CCP1Mat)

    fP_a1 = fP_a1given(TransitionMat, CCP2)
    fP_a2 = fP_a2given(TransitionMat, CCP1)

    # CCP の更新
    future1 = np.column_stack([fP_a1[k] @ ExanteV1 for k in range(3)])
    NewSigmaSeed1 = (pi1Psigma + beta * future1) * CCP1Adjuster
    NewSigmaDeno1 = np.sum(np.exp(NewSigmaSeed1), axis=1) - 1.0
    CCP1UpdatedMat_local = np.exp(NewSigmaSeed1) / NewSigmaDeno1.reshape(-1, 1) * CCP1Adjuster
    CCP1Updated = CCP1UpdatedMat_local[:, 1]

    future2 = np.column_stack([fP_a2[k] @ ExanteV2 for k in range(3)])
    NewSigmaSeed2 = (pi2Psigma + beta * future2) * CCP2Adjuster
    NewSigmaDeno2 = np.sum(np.exp(NewSigmaSeed2), axis=1) - 1.0
    CCP2UpdatedMat_local = np.exp(NewSigmaSeed2) / NewSigmaDeno2.reshape(-1, 1) * CCP2Adjuster
    CCP2Updated = CCP2UpdatedMat_local[:, 1]

    obj = np.sum((CCP1Updated - CCP1)**2) + np.sum((CCP2Updated - CCP2)**2)
    return obj


def obj_forward_PSD(x):
    \"\"\"5パラメータ -> 10パラメータへの変換ラッパー\"\"\"
    theta10 = np.array([x[0], x[2], x[3], 0, x[4], x[1], x[2], x[3], 0, x[4]])
    return Estimation_forward_PSD(theta10, W1star, W2star, EstimatedTransition,
                                   EstimatedCCP1, EstimatedCCP2, beta)


print("P-SD 推定関数 定義完了")
"""))

# ============================================================
# Cell 14: P-SD Optimization
# ============================================================
cells.append(nbf.v4.new_markdown_cell("## P-SD 推定の実行"))

cells.append(nbf.v4.new_code_cell("""\
# 初期値
initial_PSD = np.array([0.3375, 0.2375, -0.27, 0.45, -2.25])

print("P-SD 推定開始...")
start_time = time.time()
res_PSD = minimize(obj_forward_PSD, initial_PSD, method='Nelder-Mead',
                   options={'maxiter': 10000, 'xatol': 1e-10, 'fatol': 1e-10})
elapsed = time.time() - start_time

opt_forwardPSD = res_PSD.x
print(f"計算時間: {elapsed:.1f} 秒")
print(f"収束: {res_PSD.success}")
print(f"目的関数値: {res_PSD.fun:.10f}")
print(f"推定パラメータ: {opt_forwardPSD}")
print(f"真の値:        {np.array([0.3, 0.2, -0.27, 0.45, -2.1])}")
"""))

# ============================================================
# Cell 15: Bootstrap for Forward PSD
# ============================================================
cells.append(nbf.v4.new_markdown_cell(
    "## P-SD with Forward Simulation の Bootstrap\n"
    "\n"
    "市場単位でリサンプリングし、100回の Bootstrap を実行する。"
))

cells.append(nbf.v4.new_code_cell("""\
def Bootstrap_PS_forward(bootsample, beta, EVrandom, UNIrandom, InitialState,
                         NumSimMarkets, NumSimulations, NumSimPeriods):
    \"\"\"Bootstrap で P-SD (Forward Simulation) の推定を行う

    Returns: (optim_result, EstCCP1, EstCCP2, EstTransition)
    \"\"\"
    # Step 1: CCP と遷移確率の推定
    EstCCP1 = np.zeros(8)
    EstCCP2 = np.zeros(8)
    for s in range(1, 9):
        sub = bootsample[bootsample[:, 2] == s]
        EstCCP1[s-1] = np.sum(sub[:, 6] == 0) / len(sub)
        EstCCP2[s-1] = np.sum(sub[:, 7] == 0) / len(sub)

    EstTransition = np.zeros((2, 2))
    n_rows = len(bootsample)
    lag_demand = np.zeros(n_rows)
    lag_demand[1:] = bootsample[:-1, 3]
    mask = bootsample[:, 1] != 1
    data_lag = np.column_stack([bootsample[mask], lag_demand[mask]])
    for z in range(1, 3):
        sub_z = data_lag[data_lag[:, 3] == z]
        sub_zz = sub_z[sub_z[:, 8] == z]
        EstTransition[z-1, z-1] = len(sub_zz) / len(sub_z)
        EstTransition[z-1, 2-z] = 1 - EstTransition[z-1, z-1]

    # Step 2: Forward Simulation
    W1s, W2s = VSigmaGeneration(
        EstCCP1, EstCCP2, EstTransition, EVrandom,
        UNIrandom, InitialState, NumSimMarkets, NumSimulations, NumSimPeriods, beta
    )

    # Step 3: 最適化
    def obj_boot(x):
        theta10 = np.array([x[0], x[2], x[3], 0, x[4], x[1], x[2], x[3], 0, x[4]])
        return Estimation_forward_PSD(theta10, W1s, W2s, EstTransition, EstCCP1, EstCCP2, beta)

    initial = np.array([0.3, 0.2, -0.27, 0.45, -2.1])
    result = minimize(obj_boot, initial, method='Nelder-Mead',
                      options={'maxiter': 10000, 'xatol': 1e-10, 'fatol': 1e-10})
    return result, EstCCP1, EstCCP2, EstTransition


# Bootstrap の実行
np.random.seed(2023)
# 注意: R版では numBootSample=100 だが、Python版では計算時間の制約のため10に削減
# R版と同じ結果を得るには numBootSample=100 に変更してください
numBootSample = 10
NumMarketsData = 500
NumPeriodsData = 50

# Bootstrap インデックス (500 x numBootSample)
bootindex = np.random.randint(1, NumMarketsData + 1, size=(NumMarketsData, numBootSample))

bootresult_payoff = np.zeros((5, numBootSample))

print(f"P-SD Forward Bootstrap 開始 ({numBootSample} 回)...")
start_time = time.time()

for b in range(numBootSample):
    # Bootstrap サンプルの構築
    bootsample = np.zeros((NumMarketsData * NumPeriodsData, 8))
    for m in range(NumMarketsData):
        mk = bootindex[m, b]
        temp = FakeData[FakeData[:, 0] == mk]
        bootsample[m*NumPeriodsData:(m+1)*NumPeriodsData, :] = temp

    result, _, _, _ = Bootstrap_PS_forward(
        bootsample, beta, EVrandom, UNIrandom, InitialState,
        NumSimMarkets, NumSimulations, NumSimPeriods
    )
    bootresult_payoff[:, b] = result.x

    if (b + 1) % 10 == 0:
        elapsed = time.time() - start_time
        print(f"  Bootstrap {b+1}/{numBootSample} 完了 ({elapsed:.1f} 秒)")

total_elapsed = time.time() - start_time
print(f"Bootstrap 完了: {total_elapsed:.1f} 秒")
"""))

# ============================================================
# Cell 16: Forward PSD Results Table
# ============================================================
cells.append(nbf.v4.new_markdown_cell("## Forward P-SD 推定結果 (Tab 11.5)"))

cells.append(nbf.v4.new_code_cell("""\
true_param = np.array([0.3, 0.2, -0.27, 0.45, -2.1])
normalized_param = np.array([
    0.3 - (1-beta)/beta * (-0.15),
    0.2 - (1-beta)/beta * (-0.15),
    -0.27,
    0.45,
    -2.1 + (-0.15)
])
se_PSD = np.std(bootresult_payoff, axis=1, ddof=0)

param_names = ['theta_base1', 'theta_base2', 'theta_rival', 'theta_boom', 'theta_entry']

print("=" * 70)
print("Forward P-SD 推定結果 (Tab 11.5)")
print("=" * 70)
print(f"{'Parameter':<14} {'True':>10} {'Normalized':>12} {'Estimated':>12} {'SE':>10}")
print("-" * 70)
for i in range(5):
    print(f"{param_names[i]:<14} {true_param[i]:10.4f} {normalized_param[i]:12.4f} "
          f"{opt_forwardPSD[i]:12.4f} {se_PSD[i]:10.4f}")

# CSV 保存
estmat_PSD_forward = np.column_stack([true_param, normalized_param, opt_forwardPSD, se_PSD])
df_PSD = pd.DataFrame(estmat_PSD_forward,
                       columns=['True', 'Normalized', 'Estimated', 'SE'],
                       index=range(1, 6))
df_PSD.to_csv(output_dir / 'Tab11_5_forward_PSD.csv')
print(f"\\n保存しました: output/Tab11_5_forward_PSD.csv")
"""))

# ============================================================
# Cell 17: BBL Setup
# ============================================================
cells.append(nbf.v4.new_markdown_cell(
    "## BBL 不等式推定量のセットアップ\n"
    "\n"
    "計算速度の観点から、BBL ではシミュレーション期間を 30 期間に短縮する。\n"
    "また、CCP を Perturbation して 200 個の代替政策を生成する。"
))

cells.append(nbf.v4.new_code_cell("""\
# BBL 用パラメータ
NumSimPeriods_BBL = 30

# 乱数の再生成
np.random.seed(2023)
total_ev_BBL = NumSimMarkets * NumSimPeriods_BBL * NumSimFirms * NumSimulations * 8 * 3
EVrandom_BBL = -np.log(-np.log(np.random.uniform(size=total_ev_BBL)))
EVrandom_BBL = EVrandom_BBL.reshape(
    (NumSimMarkets, NumSimPeriods_BBL, NumSimFirms, NumSimulations, 8, 3), order='F'
)

total_uni_BBL = NumSimMarkets * NumSimPeriods_BBL * NumSimulations
UNIrandom_BBL = np.random.uniform(size=total_uni_BBL)
UNIrandom_BBL = UNIrandom_BBL.reshape(
    (NumSimMarkets, NumSimPeriods_BBL, NumSimulations), order='F'
)

print(f"BBL用乱数生成完了")
print(f"EVrandom_BBL shape: {EVrandom_BBL.shape}")
print(f"UNIrandom_BBL shape: {UNIrandom_BBL.shape}")
"""))

# ============================================================
# Cell 18: Forward Sim with estimated CCPs (BBL version)
# ============================================================
cells.append(nbf.v4.new_markdown_cell("## 推定 CCP の下での Forward Simulation (BBL 用)"))

cells.append(nbf.v4.new_code_cell("""\
print("BBL 用 Forward Simulation (推定 CCP)...")
start_time = time.time()
W1star_BBL, W2star_BBL = VSigmaGeneration(
    EstimatedCCP1, EstimatedCCP2, EstimatedTransition, EVrandom_BBL,
    UNIrandom_BBL, InitialState, NumSimMarkets, NumSimulations, NumSimPeriods_BBL, beta
)
elapsed = time.time() - start_time
print(f"計算時間: {elapsed:.1f} 秒")
"""))

# ============================================================
# Cell 19: CCP Perturbation
# ============================================================
cells.append(nbf.v4.new_markdown_cell("## CCP の Perturbation (200 回)"))

cells.append(nbf.v4.new_code_cell("""\
NumPerturbations = 200

# Perturbation 用の乱数
# R の set.seed(2023) に対応する乱数はここでは独自に生成
# (Matlab の .mat ファイルがないため)
np.random.seed(2023)
noise_CCP1_raw = np.random.normal(0, 0.1, size=(8, NumPerturbations))
noise_CCP2_raw = np.random.normal(0, 0.1, size=(8, NumPerturbations))

PerturbedCCP1 = np.tile(EstimatedCCP1.reshape(-1, 1), (1, NumPerturbations)) + noise_CCP1_raw
PerturbedCCP2 = np.tile(EstimatedCCP2.reshape(-1, 1), (1, NumPerturbations)) + noise_CCP2_raw

# [0.001, 0.999] にクリップ
PerturbedCCP1 = np.clip(PerturbedCCP1, 0.001, 0.999)
PerturbedCCP2 = np.clip(PerturbedCCP2, 0.001, 0.999)

print(f"Perturbed CCP 生成完了: shape = {PerturbedCCP1.shape}")
print(f"PerturbedCCP1 range: [{PerturbedCCP1.min():.4f}, {PerturbedCCP1.max():.4f}]")
print(f"PerturbedCCP2 range: [{PerturbedCCP2.min():.4f}, {PerturbedCCP2.max():.4f}]")
"""))

# ============================================================
# Cell 20: Forward Simulation for all perturbations
# ============================================================
cells.append(nbf.v4.new_markdown_cell(
    "## Perturbation した CCP での Forward Simulation\n"
    "\n"
    "200 個の Perturbed CCP それぞれに対して Forward Simulation を実行する。\n"
    "この処理は時間がかかる。"
))

cells.append(nbf.v4.new_code_cell("""\
W1_all = np.zeros((6, NumSimMarkets, NumPerturbations))
W2_all = np.zeros((6, NumSimMarkets, NumPerturbations))

print(f"BBL Perturbation Forward Simulation 開始 ({NumPerturbations} 回)...")
start_time = time.time()

for per in range(NumPerturbations):
    # 企業1の CCP を Perturb して Forward Simulation
    W1_p, _ = VSigmaGeneration(
        PerturbedCCP1[:, per], EstimatedCCP2, EstimatedTransition,
        EVrandom_BBL, UNIrandom_BBL, InitialState,
        NumSimMarkets, NumSimulations, NumSimPeriods_BBL, beta
    )
    W1_all[:, :, per] = W1_p

    # 企業2の CCP を Perturb して Forward Simulation
    _, W2_p = VSigmaGeneration(
        EstimatedCCP1, PerturbedCCP2[:, per], EstimatedTransition,
        EVrandom_BBL, UNIrandom_BBL, InitialState,
        NumSimMarkets, NumSimulations, NumSimPeriods_BBL, beta
    )
    W2_all[:, :, per] = W2_p

    if (per + 1) % 10 == 0:
        elapsed = time.time() - start_time
        print(f"  Perturbation {per+1}/{NumPerturbations} 完了 ({elapsed:.1f} 秒)")

total_elapsed = time.time() - start_time
print(f"Perturbation Forward Simulation 完了: {total_elapsed:.1f} 秒")
"""))

# ============================================================
# Cell 21: BBL Objective Function
# ============================================================
cells.append(nbf.v4.new_markdown_cell("## BBL 目的関数の定義"))

cells.append(nbf.v4.new_code_cell("""\
def BBLobjective_NLS(theta10, NumPerturbations, W1star, W2star, W1_all, W2_all):
    \"\"\"BBL 不等式推定量の目的関数 (NLS 形式)

    均衡 CCP での Value >= Perturbed CCP での Value
    を不等式制約として、違反分の二乗和を最小化する。
    \"\"\"
    Param1 = np.append(theta10[:5], 1.0)
    Param2 = np.append(theta10[5:10], 1.0)

    # W1star を NumPerturbations 回複製
    # fa1: (6, 8, NumPerturbations)
    fa1 = np.tile(W1star[:, :, np.newaxis], (1, 1, NumPerturbations))
    fa2 = np.tile(W2star[:, :, np.newaxis], (1, 1, NumPerturbations))

    diff1 = fa1 - W1_all  # (6, 8, NumPerturbations)
    diff2 = fa2 - W2_all

    # (6, 8*NumPerturbations) に reshape
    temp1 = diff1.reshape(6, -1, order='F')
    temp2 = diff2.reshape(6, -1, order='F')

    val1 = temp1.T @ Param1  # (8*NumPerturbations,)
    val2 = temp2.T @ Param2

    # max(0, -val) → pmin(val, 0) の二乗和
    val1 = np.minimum(val1, 0)
    val2 = np.minimum(val2, 0)

    obj = np.sum(val1**2) + np.sum(val2**2)
    return obj


def obj_forward_BBL_NLS(x):
    \"\"\"5パラメータ -> 10パラメータへの変換ラッパー (BBL)\"\"\"
    theta10 = np.array([x[0], x[2], x[3], 0, x[4], x[1], x[2], x[3], 0, x[4]])
    return BBLobjective_NLS(theta10, NumPerturbations, W1star_BBL, W2star_BBL, W1_all, W2_all)


print("BBL 目的関数 定義完了")
"""))

# ============================================================
# Cell 22: BBL Optimization
# ============================================================
cells.append(nbf.v4.new_markdown_cell("## BBL 推定の実行"))

cells.append(nbf.v4.new_code_cell("""\
initial_BBL = np.array([0.3, 0.2, -0.27, 0.45, -2.1])

print("BBL 推定開始...")
start_time = time.time()
opt_BBL = minimize(obj_forward_BBL_NLS, initial_BBL, method='Nelder-Mead',
                   options={'maxiter': 10000, 'xatol': 1e-10, 'fatol': 1e-10})
elapsed = time.time() - start_time

print(f"計算時間: {elapsed:.1f} 秒")
print(f"収束: {opt_BBL.success}")
print(f"目的関数値: {opt_BBL.fun:.10f}")
print(f"推定パラメータ: {opt_BBL.x}")
print(f"真の値:        {np.array([0.3, 0.2, -0.27, 0.45, -2.1])}")
"""))

# ============================================================
# Cell 23: Bootstrap for BBL
# ============================================================
cells.append(nbf.v4.new_markdown_cell(
    "## BBL の Bootstrap\n"
    "\n"
    "BBL 推定を 100 回の Bootstrap で標準誤差を計算する。\n"
    "注意: この処理は非常に長い時間がかかる。"
))

cells.append(nbf.v4.new_code_cell("""\
def Bootstrap_BBL(bootsample, beta, EVrandom, UNIrandom, InitialState,
                  NumSimMarkets, NumSimulations, NumSimPeriods,
                  NumPerturbations, noise_CCP1, noise_CCP2):
    \"\"\"BBL の Bootstrap 1 回分\"\"\"
    # Step 1: CCP と遷移確率の推定
    EstCCP1 = np.zeros(8)
    EstCCP2 = np.zeros(8)
    for s in range(1, 9):
        sub = bootsample[bootsample[:, 2] == s]
        EstCCP1[s-1] = np.sum(sub[:, 6] == 0) / len(sub)
        EstCCP2[s-1] = np.sum(sub[:, 7] == 0) / len(sub)

    EstTransition = np.zeros((2, 2))
    n_rows = len(bootsample)
    lag_demand = np.zeros(n_rows)
    lag_demand[1:] = bootsample[:-1, 3]
    mask = bootsample[:, 1] != 1
    data_lag = np.column_stack([bootsample[mask], lag_demand[mask]])
    for z in range(1, 3):
        sub_z = data_lag[data_lag[:, 3] == z]
        sub_zz = sub_z[sub_z[:, 8] == z]
        EstTransition[z-1, z-1] = len(sub_zz) / len(sub_z)
        EstTransition[z-1, 2-z] = 1 - EstTransition[z-1, z-1]

    # Step 2-1: Forward Simulation with estimated CCP
    W1s, W2s = VSigmaGeneration(
        EstCCP1, EstCCP2, EstTransition, EVrandom,
        UNIrandom, InitialState, NumSimMarkets, NumSimulations, NumSimPeriods, beta
    )

    # Step 2-2: Perturbed CCP での Forward Simulation
    PertCCP1 = np.tile(EstCCP1.reshape(-1, 1), (1, NumPerturbations)) + noise_CCP1
    PertCCP2 = np.tile(EstCCP2.reshape(-1, 1), (1, NumPerturbations)) + noise_CCP2
    PertCCP1 = np.clip(PertCCP1, 0.001, 0.999)
    PertCCP2 = np.clip(PertCCP2, 0.001, 0.999)

    W1_all_b = np.zeros((6, NumSimMarkets, NumPerturbations))
    W2_all_b = np.zeros((6, NumSimMarkets, NumPerturbations))

    for per in range(NumPerturbations):
        W1_p, _ = VSigmaGeneration(
            PertCCP1[:, per], EstCCP2, EstTransition,
            EVrandom, UNIrandom, InitialState,
            NumSimMarkets, NumSimulations, NumSimPeriods, beta
        )
        W1_all_b[:, :, per] = W1_p

        _, W2_p = VSigmaGeneration(
            EstCCP1, PertCCP2[:, per], EstTransition,
            EVrandom, UNIrandom, InitialState,
            NumSimMarkets, NumSimulations, NumSimPeriods, beta
        )
        W2_all_b[:, :, per] = W2_p

    # Step 3: 最適化
    def obj_boot_BBL(x):
        theta10 = np.array([x[0], x[2], x[3], 0, x[4], x[1], x[2], x[3], 0, x[4]])
        return BBLobjective_NLS(theta10, NumPerturbations, W1s, W2s, W1_all_b, W2_all_b)

    initial = np.array([0.3, 0.2, -0.27, 0.45, -2.1])
    result = minimize(obj_boot_BBL, initial, method='L-BFGS-B',
                      bounds=[(0, None), (0, None), (None, 0), (0, None), (None, 0)],
                      options={'maxiter': 100000, 'ftol': 1e-10})
    return result, EstCCP1, EstCCP2, EstTransition


# Bootstrap の実行
np.random.seed(2023)
# 注意: R版では numBootSample_BBL=100 だが、Python版では計算時間の制約のため10に削減
# R版と同じ結果を得るには numBootSample_BBL=100 に変更してください
numBootSample_BBL = 10

# Bootstrap インデックスの再生成
bootindex_BBL = np.random.randint(1, NumMarketsData + 1, size=(NumMarketsData, numBootSample_BBL))

# Perturbation 用のノイズ (Bootstrap 全体で共通)
np.random.seed(2023)
noise_CCP1_boot = np.random.normal(0, 0.1, size=(8, NumPerturbations))
noise_CCP2_boot = np.random.normal(0, 0.1, size=(8, NumPerturbations))

bootresult_BBL = np.zeros((5, numBootSample_BBL))

print(f"BBL Bootstrap 開始 ({numBootSample_BBL} 回)...")
print("注意: 各 Bootstrap サンプルで 200 回の Forward Simulation を行うため非常に時間がかかります")
start_time = time.time()

for b in range(numBootSample_BBL):
    # Bootstrap サンプルの構築
    bootsample = np.zeros((NumMarketsData * NumPeriodsData, 8))
    for m in range(NumMarketsData):
        mk = bootindex_BBL[m, b]
        temp = FakeData[FakeData[:, 0] == mk]
        bootsample[m*NumPeriodsData:(m+1)*NumPeriodsData, :] = temp

    result_b, _, _, _ = Bootstrap_BBL(
        bootsample, beta,
        EVrandom_BBL, UNIrandom_BBL, InitialState,
        NumSimMarkets, NumSimulations, NumSimPeriods_BBL,
        NumPerturbations, noise_CCP1_boot, noise_CCP2_boot
    )
    bootresult_BBL[:, b] = result_b.x

    if (b + 1) % 10 == 0:
        elapsed = time.time() - start_time
        print(f"  Bootstrap {b+1}/{numBootSample_BBL} 完了 ({elapsed:.1f} 秒)")

total_elapsed = time.time() - start_time
print(f"BBL Bootstrap 完了: {total_elapsed:.1f} 秒")
"""))

# ============================================================
# Cell 24: BBL Results Table
# ============================================================
cells.append(nbf.v4.new_markdown_cell("## BBL 推定結果 (Tab 11.5)"))

cells.append(nbf.v4.new_code_cell("""\
true_param = np.array([0.3, 0.2, -0.27, 0.45, -2.1])
normalized_param = np.array([
    0.3 - (1-beta)/beta * (-0.15),
    0.2 - (1-beta)/beta * (-0.15),
    -0.27,
    0.45,
    -2.1 + (-0.15)
])
se_BBL = np.std(bootresult_BBL, axis=1, ddof=0)

param_names = ['theta_base1', 'theta_base2', 'theta_rival', 'theta_boom', 'theta_entry']

print("=" * 70)
print("BBL 不等式推定量の結果 (Tab 11.5)")
print("=" * 70)
print(f"{'Parameter':<14} {'True':>10} {'Normalized':>12} {'Estimated':>12} {'SE':>10}")
print("-" * 70)
for i in range(5):
    print(f"{param_names[i]:<14} {true_param[i]:10.4f} {normalized_param[i]:12.4f} "
          f"{opt_BBL.x[i]:12.4f} {se_BBL[i]:10.4f}")

# CSV 保存
estmat_BBL = np.column_stack([true_param, normalized_param, opt_BBL.x, se_BBL])
df_BBL = pd.DataFrame(estmat_BBL,
                       columns=['True', 'Normalized', 'Estimated', 'SE'],
                       index=range(1, 6))
df_BBL.to_csv(output_dir / 'Tab11_5_forward_BBL.csv')
print(f"\\n保存しました: output/Tab11_5_forward_BBL.csv")
"""))

# ============================================================
# Cell 25: Summary
# ============================================================
cells.append(nbf.v4.new_markdown_cell("## まとめ"))

cells.append(nbf.v4.new_code_cell("""\
print("=" * 70)
print("第11章: Forward Simulation を用いた P-SD 推定および BBL 不等式推定量")
print("=" * 70)

print("\\n--- Forward P-SD 推定結果 ---")
print(f"{'Parameter':<14} {'True':>10} {'Estimated':>12} {'SE':>10}")
print("-" * 50)
for i in range(5):
    print(f"{param_names[i]:<14} {true_param[i]:10.4f} {opt_forwardPSD[i]:12.4f} {se_PSD[i]:10.4f}")

print("\\n--- BBL 不等式推定量の結果 ---")
print(f"{'Parameter':<14} {'True':>10} {'Estimated':>12} {'SE':>10}")
print("-" * 50)
for i in range(5):
    print(f"{param_names[i]:<14} {true_param[i]:10.4f} {opt_BBL.x[i]:12.4f} {se_BBL[i]:10.4f}")

print("\\n出力ファイル:")
print("  - Tab11_5_forward_PSD.csv")
print("  - Tab11_5_forward_BBL.csv")
print("\\n完了")
"""))

nb.cells = cells
nbf.write(nb, 'main_ch11_forward_BBL.ipynb')
print("Generated: main_ch11_forward_BBL.ipynb")
