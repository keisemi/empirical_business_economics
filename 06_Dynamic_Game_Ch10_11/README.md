# 動学ゲームの推定コードについて

## はじめに

- 第10章・11章の「動学ゲームの推定」については、Matlab、R、及びPython(Jupyter Notebook)でコードを提供している。経済セミナー連載時(第9回から第12回)にはMatlabを用いて数値計算を行っていたが、書籍化にあたってRでコードを作成した。さらにPython版をJupyter Notebookとして作成した。
- しかしながら、連載時の結果と同じ結果を得るために、Matlabで生成したデータ(乱数を含む)をRで読み込んで用いる形にしている。
- それでもなお、一部の結果については、連載時の結果と異なる場合があることに注意されたい。  
  
## ディレクトリの構成

- `R/`: Rで作成したコード群
  - `main1_Computation_Equilibrium.R` 〜 `main5_Policy_Simulation.R`: メインコード
  - `sub_X_YYYY.R`: サブコード群。メインコードから呼び出しされる
- `python/`: Python(Jupyter Notebook)版のコード群
  - `generate_notebook_ch10_ch11.py`: 全ノートブックを一括生成・実行するマスタースクリプト
  - `generate_notebook_*.py`: 各ノートブックを生成するスクリプト
  - `main_ch10_equilibrium.ipynb` 〜 `main_ch11_policy_sim.ipynb`: 生成されたJupyter Notebook
- `code_from_matlab/`: 連載時に使用したMatlabのコード。説明は後述
- `data_from_matlab/`: 連載時にMatlabで生成したデータ。R・Pythonコードで読み込んで使用する
- `functions_R/`: Rで作成した関数群
- `output/`: R・Pythonコードのアウトプットを保存するディレクトリ

## 実行する順番 (R版)

1. `data_from_matlab`に以下のファイルがあること確認する。
  - `FakeData_Matlab.csv`: 疑似データ. これはGithub上にも保存されている。
  - `random_number_matlab_BBL.mat`: BBLの推定に用いる乱数。ファイルサイズから別場所にアップロードされている。
  - `random_number_matlab_PSD.mat`: P-SDの推定に用いる乱数。ファイルサイズから別場所にアップロードされている。
2. `R/main1_Computation_Equilibrium.R`: 均衡計算を行うコード
3. `R/main2_Estimation_AM.R`: Aguirregabiria and Mira (2007)の方法による推定を行うコード。
4. `R/main3_Estimation_PSD.R`: P-SDによる推定を行うコード
5. `R/main4_Estimation_Forward_PSD_BBL.R`: Forward simulation を用いたP-SDとBBLによる推定を行うコード
6. `R/main5_Policy_Simulation.R`: 反実仮想分析を行うコード

## 実行する順番 (Python版)

Python版はJupyter Notebookとして提供されている。`python/`ディレクトリで以下を実行する。

### ノートブックの一括生成

```bash
python generate_notebook_ch10_ch11.py
```

### ノートブックの一括生成・実行

```bash
python generate_notebook_ch10_ch11.py --execute        # 全ノートブックを生成・実行
python generate_notebook_ch10_ch11.py --execute 1 5    # 全ノートブックを生成し、#1と#5のみ実行
```

### 各ノートブックの内容

1. `main_ch10_equilibrium.ipynb`: 均衡計算と疑似データの生成
2. `main_ch11_AM.ipynb`: Aguirregabiria and Mira (2007)の方法による推定
3. `main_ch11_PSD.ipynb`: P-SDによる推定
4. `main_ch11_forward_BBL.ipynb`: Forward simulationを用いたP-SDとBBLによる推定
5. `main_ch11_policy_sim.ipynb`: 反実仮想分析

### Python版の注意事項

- Python版ではR版と同じMatlabデータ(`data_from_matlab/`)を使用する。
- BBLのブートストラップ標準誤差の計算は、R版が並列計算(10コア)を用いるのに対し、Python版では逐次計算のため、ブートストラップ回数をR版の100回から10回に削減している。
- 一部の結果(ブートストラップ標準誤差など)はR版と異なる場合がある。

## 計算時間

- デフォルト環境: 
  - MacBook Air M2 2022 (Apple M2, 8コアCPU, メモリ16GB) 
  - Sequioia 15.6
  - R 4.4.3

- ワークステーション環境: 
  - Mac Studio 2025 Apple M4 Max CPU 16 コア CPU メモリ 128GB 
  - Sequioia 15.6
  - R 4.4.3

- main1,2,3,5はデフォルト環境で実行した。(並列計算なし)
  - main1,3,5は数分程度(10分未満)
  - main2 は15-20分程度

- main4はワークステーション環境で実行した。
  - 並列計算に用いたコア数は10.
  - Forward simulationを用いたP-SD: 5分未満
  - BBLの点推定値: 5分未満
  - BBLのブートストラップ標準誤差: 約1時間半


### Matlabコードについて

- ディレクトリ`code_from_matlab`に連載時に使用したMatlabコードを保存している。第12回 (https://sites.google.com/view/keisemi-ebiz/%E7%AC%AC12%E5%9B%9E?authuser=0) のコードを若干修正したもの。なお、第12回のコードには、通常のP-SD (Forward simulationをしないP-SD)のコードが含まれていない。こちらは第11回 (https://sites.google.com/view/keisemi-ebiz/%E7%AC%AC11%E5%9B%9E?authuser=0) を参照されたい。
- MatlabコードはRコードのためのデータ生成に用いるのみで、Rコードの実行には不要。
- 実行する際には以下の順番で行う。


1. `main1_Estimation.m`
  - 疑似データの作成、Forward simulationに基づくP-SDによる推定、BBLによる推定を行う。
  - アウトプット: `result_BBL_bootstrap.mat`
  - アウトプット: `FakeData_Matlab.csv` `random_number_matlab_BBL.mat` `random_number_matlab_PSD.mat`
2. `main2_Policy_Simulation.m`
  - アウトプット: `ProbEntry.jpg`

もし連載時におけるBBLの標準誤差を計算したい場合には、以下を実行すること。

```
library(R.matlab)
data <- readMat("06_Dynamic_Game_Ch10_11/code_from_matlab/result_BBL_bootstrap.mat")
apply(data$bootresult.payoff,2, sd)
```
