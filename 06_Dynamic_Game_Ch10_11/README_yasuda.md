# Dynamic Game

## 注意

-   本`README.md`はあくまで一般公開を想定したものではなく、編集者向けの内容であるため、編集者向けの留意点も含むことに注意すること。

-   以下、連載当初用いられていたパラメターを旧パラメター、新しく用いようとしたパラメターを新パラメターと呼ぶ。

## ディレクトリ構成

-   `AM2007_R`: Aguirregabiria and Mira (2007) に基づく推定
    -   概ね経済セミナー第11回の[サポートコード](https://sites.google.com/view/keisemi-ebiz/%E7%AC%AC11%E5%9B%9E?authuser=0)と同内容
    -   `main.R`の`4. Aguirregabiria and Mira (2007)の方法によるパラメターの推定`以降の内容
    -   追加した自作関数群
        -   `check_convergence.R`: CCPと更新されたCCPが収束しているか確認する関数。

        -   `Estimation_AM_bootstrap.R`: Aguirregabiria and Mira (2007)の推定をブートストラップで繰り返すために必要な処理をまとめた関数。従来の推定より結果の標準出力を極力省略している。

        -   `obj_lik.R`: 目的関数として疑似対数尤度関数を定義した関数。
-   `KS12_R`: 経済セミナー第12回に基づく反実仮想分析
    -   概ね経済セミナー第12回の[サポートコード](https://sites.google.com/view/keisemi-ebiz/%E7%AC%AC12%E5%9B%9E?authuser=0)と同内容
    - `main1_Estimation.R`の内容はMatLab版データで実行されている（本誌サポートコードではRで生成したデータを用いていた）。
    - `main2_Policy_Simulation.R`は本誌とは異なるパラメターのもとで反実仮想分析のプロットを生成している。
    -   `result`内にはCounter Factualで出力された店舗存在確率のグラフが格納されている。
        -   `ProbEntry.png`: 何の意味もなさないファイル
        -   `ProbEntry1.png`: 新パラメターに基づく企業1のCounter Factual
        -   `ProbEntry1Original.png`: 旧パラメターに基づく企業1のCounter Factual
        -   `ProbEntry2.png`: 新パラメターに基づく企業2のCounter Factual
        -   `ProbEntry2Original.png`: 旧パラメターに基づく企業2のCounter Factual
        -   `ProbEntry3Plots`: 企業ごとにベースライン、シナリオ1、シナリオ2をまとめたグラフ。実線がベースライン、点線に白丸がシナリオ1、点線に黒丸がシナリオ2を指している

## Aguirregabiria and Mira (2007)

-   アルゴリズムを実行するにあたり、パラメターとCCPについては初期値が必要となる。ここでは、以下のように初期値を設定した。

    -   パラメターについて、真のパラメターに近い`InitialParameters1`と真のパラメターから遠い`InitialParameters2`が存在する。`InitialParameters2`では推定が収束しなかったため、`InitialParameters1`を用いた。

    -   CCPについて、実装当初はデータから推定されたCCPから遠いCCPも初期値として設定していたが、収束しなかったため不採用とし、データから推定されたCCPをそのまま用いた。

-   なお、ブートストラップ法による標準誤差の出力においても同様の初期値を用いている。詳細は`EM_AM_bootstrap.R`を参照のこと。

## Policy Simulation

-   シナリオ1

    -   企業1のベース利潤を`0.3`から`0.5`に、顧客奪取効果を企業1, 2共に`-0.27`から`0`に変更したものである。

-   シナリオ2

    -   企業1のベース利潤を`0.3`から`0.5`に、

    -   ライバルの店舗数が企業1の利潤に与える影響を`-0.27`から`-0.1`に、

    -   ライバルの店舗数が企業2の利潤に与える影響を`-0.27`から`-0.2`に変更したものである。

## Pending

-   `AM2007` の`main.R` について、新パラメターを記載しつつも採用しなかったという表記にしてある。存在自体を抹消する場合、該当部分を削除する必要がある。

-   `AM2007` の`main.R` について、`set.seed(2023)` を複数回用いている。特に推定上の問題はないが、冗長である。

-   `KS12_R` のCounter Factualについて、採用したプロット以外のプロットが不要である場合、そのプロットのコードを削除することを推奨する。
