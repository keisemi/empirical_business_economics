# 実証ビジネス・エコノミクス — Python版

> 本リポジトリは [keisemi/empirical_business_economics](https://github.com/keisemi/empirical_business_economics) の fork です。原著のRコードを **Python (Jupyter Notebook)** に変換したものを提供しています。原著のRコード・データ・オンライン補足もそのまま含まれています。

<p align="center">
  <img src=".\src\09633.jpg" width="250">
</p>

# 書籍情報

- 上武康亮・遠山祐太・若森直樹・渡辺安虎/著『[実証ビジネス・エコノミクス](https://nippyo.co.jp/shop/book/9633.html)』（日本評論社、2025年12月刊）
- **経済理論とデータの力で、ビジネスはもっと強くなる！**
  - 理論とデータを融合した「**構造推定**」。 ビジネスの現場でも活躍するこの実証手法を、 製品やサービスのプライシング（価格戦略）、ライバル会社と合併、新規事業への参入や既存事業からの退出の意思決定、ブランド戦略などなど、経営戦略・マーケティング戦略を考える豊富な具体例とともに実践的に解説！
- 本書の「**はしがき**」や**内容紹介**などを公開しています【 **[リンク](https://note.com/keisemi/n/ne05895107194)** 】
- 本書の一部をPDFで「**立ち読み**」いただけます【 **[リンク](https://nippyo.co.jp/shop/img/content_pdf/09633.pdf)** 】(PDFが開きます)
- サポート情報のご案内【 **[リンク](https://keisemi.github.io/empirical_business_economics/)** 】

# Python Notebook 一覧

| 章 | トピック | ディレクトリ | ノートブック |
|---|---|---|---|
| 第2章 | 離散選択モデル入門 | `01_Discrete_Choice_Ch02/python/` | `main_ch02.ipynb` |
| 第3章 | ロジットモデルによる需要推定 | `02_BLP_Ch03_04_05/python/` | `main_ch03.ipynb` |
| 第4章 | BLPモデルによる需要推定 | `02_BLP_Ch03_04_05/python/` | `main_ch04.ipynb` |
| 第5章 | 需要推定の応用：価格戦略・合併シミュレーション | `02_BLP_Ch03_04_05/python/` | `main_ch05.ipynb` |
| 第6章 | 動的単一エージェントモデルの推定 | `03_Dynamic_Single_Agent_Ch06_07/python/` | `main_ch06.ipynb` |
| 第7章 | 動的単一エージェントモデルの反実仮想シミュレーション | `03_Dynamic_Single_Agent_Ch06_07/python/` | `main_ch07.ipynb` |
| 第8章(1) | 静的参入ゲーム：Bresnahan & Reiss (1991) | `04_Static_Game_MRI_Ch08/python/` | `main_ch08_01_BR1991.ipynb` |
| 第8章(2) | 静的参入ゲーム：Berry (1992) | `04_Static_Game_MRI_Ch08/python/` | `main_ch08_02_Berry1992.ipynb` |
| 第9章 | 静的ゲームの応用：航空路線参入 | `05_Static_Game_Airline_Ch09/python/` | `main_ch09.ipynb` |
| 第10章 | 動的ゲームの均衡計算 | `06_Dynamic_Game_Ch10_11/python/` | `main_ch10_equilibrium.ipynb` |
| 第11章(1) | 動的ゲームの推定：Aguirregabiria & Mira (2007) | `06_Dynamic_Game_Ch10_11/python/` | `main_ch11_AM.ipynb` |
| 第11章(2) | 動的ゲームの推定：Pesendorfer & Schmidt-Dengler (2008) | `06_Dynamic_Game_Ch10_11/python/` | `main_ch11_PSD.ipynb` |
| 第11章(3) | 動的ゲームの推定：Forward-looking BBL | `06_Dynamic_Game_Ch10_11/python/` | `main_ch11_forward_BBL.ipynb` |
| 第11章(4) | 動的ゲームの政策シミュレーション | `06_Dynamic_Game_Ch10_11/python/` | `main_ch11_policy_sim.ipynb` |

# 実行方法

各章の `python/` ディレクトリ内で以下の手順で実行できます。

```bash
# 依存ライブラリのインストール
pip install -r requirements.txt

# ノートブックの生成
python generate_notebook_chXX.py

# ノートブックの実行（実行済みノートブックを出力）
jupyter nbconvert --to notebook --execute --ExecutePreprocessor.timeout=7200 \
  --output main_chXX_executed.ipynb main_chXX.ipynb
```

生成済みの `*_executed.ipynb` も同梱しているため、コードを実行せずに結果を確認することもできます。

# 実行環境

- Python 3.11
- 主な依存ライブラリ: NumPy, SciPy, pandas, statsmodels, matplotlib, linearmodels

各章の `python/requirements.txt` に必要なパッケージが記載されています。

# 注意事項

- 一部の章ではブートストラップ回数を削減しています（Ch08_02: 100→20, Ch11 forward BBL: 100→10）。計算結果はR版と若干異なる場合があります。
- 原著のRコード・データ・オンライン補足はそのまま含まれています。Rコードの実行環境については[原著リポジトリ](https://github.com/keisemi/empirical_business_economics)を参照してください。

# ご利用に際してのお断り

本GitHubリポジトリや本書（『実証ビジネス・エコノミクス』）の内容、およびサンプルコード等の資料は、情報提供のみを目的としています。運用に際しては十分にご確認をいただき、お客様ご自身の責任とご判断に基づいて行ってください。これらの情報を運用した結果により損害等が生じた場合でも、著者・日本評論社はいかなる責任も負うことはできませんので、ご留意ください。
