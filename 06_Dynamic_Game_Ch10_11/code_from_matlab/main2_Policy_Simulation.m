% 経済セミナー連載「実証ビジネス・エコノミクス」
% 第12回「高級化路線の長期的価値:動学ゲームの推定［入門編4］」
% 上武康亮・遠山祐太・若森直樹・渡辺安虎
% 最終更新：2023年5月25日

%% 0. はじめに %%
% このmファイルは経済セミナー連載「実証ビジネス・エコノミクス」における
% 第12回「高級化路線の長期的価値:動学ゲームの推定［入門編4］」
% に付随するmatlabのmファイルの紹介になります。

% 本連載の内容、およびサンプルコード等の資料は、情報の提供のみを目的として
% いますので、運用につきましては十分にご確認をいただき、お客様ご自身の責任
% とご判断によって行ってください。これらの情報の運用結果により損害が生じた
% 場合でも、日本評論社および著者はいかなる責任を負うことはできませんので、
% ご留意ください。

% 今回の連載におけるメインファイルは以下の２つになっています。
% main1_Estimation.m: パラメタの推定
% main2_Policy_Simulation.m: 反実仮想シミュレーション（本ファイル）

% 本ファイルの流れは以下のようになっています
% 1. 下準備
% 2. ベースラインの均衡計算
% 3. 反実仮想シミュレーション: 新ブランド展開

%% 1. 下準備 %%

% データをクリアする
clear
close all

% 割引因子
beta = 0.8;

% 景気の遷移行列
TransitionMat = [ 0.7 0.3 ; 0.4 0.6];

% 今回のモデルにおける、各企業の持っている全選択要素数
NumMaxChoices = 3;

% パラメターの設定 (ベースライン)
BaselineParameterValues = [ 0.3;  % 企業1のベース利潤
    -0.27;  % ライバルの店舗数が企業1の利潤に与える影響
    0.45;  % 景気が良い時の企業1への追加的利潤
    -0.15;  % 企業1の退出のためのコスト
    -2.1;  % 企業1の出店のためのコスト
    0.2;  % 企業2のベース利潤
    -0.27;  % ライバルの店舗数が企業2の利潤に与える影響
    0.45;  % 景気が良い時の企業2への追加的利潤
    -0.15;  % 企業2の退出のためのコスト
    -2.10 ]; % 企業2の出店のためのコスト

CCP1Adjuster = [ 0 1 1 ; 0 1 1 ; 1 1 0 ; 1 1 0 ; ...
    0 1 1 ; 0 1 1 ; 1 1 0 ; 1 1 0 ];
CCP2Adjuster = [ 0 1 1 ; 1 1 0 ; 0 1 1 ; 1 1 0 ; ...
    0 1 1 ; 1 1 0 ; 0 1 1 ; 1 1 0 ];

pi1 = pi1gen(BaselineParameterValues).*CCP1Adjuster;
pi2 = pi2gen(BaselineParameterValues).*CCP2Adjuster;


%% 2. ベースラインの均衡

% MPEを解く
[CCP1base, CCP2base, V1base, V2base] = f_MPE(TransitionMat, pi1, pi2, beta);

% dataをシミュレーションする
NumSimPeriods = 15;

% 均衡における状態の遷移行列を計算する。fPsigma.m を利用する。
fPsigma = fP(TransitionMat, CCP1base(:,2), CCP2base(:,2));

% 初期の状態をセットする。n1 = n2 = 0, z = G
initial_state = [1; zeros(7,1)];

transitionpath = zeros(NumSimPeriods, 8);

for tt = 1:NumSimPeriods
    
    if (tt == 1)
        transitionpath(tt, :) = initial_state'; 
    elseif tt >= 2 
        transitionpath(tt, :) = ( fPsigma'*transitionpath(tt-1,:)' )';
    end
end

% 参入状況をカウントする。
n1 = sum( transitionpath.*repmat( [0,0,1,1,0,0,1,1],NumSimPeriods,1 ), 2);
n2 = sum( transitionpath.*repmat( [0,1,0,1,0,1,0,1],NumSimPeriods,1 ), 2);


    
%% 3. 反実仮想シミュレーション: 新ブランド展開

% 新ブランド展開により企業１のベース利潤が0.6になると想定する。  

CounterfactualParameterValues = [ 0.6;  % 企業1のベース利潤
    -0.27;  % ライバルの店舗数が企業1の利潤に与える影響
    0.45;  % 景気が良い時の企業1への追加的利潤
    -0.15;  % 企業1の退出のためのコスト
    -2.1;  % 企業1の出店のためのコスト
    0.2;  % 企業2のベース利潤
    -0.27;  % ライバルの店舗数が企業2の利潤に与える影響
    0.45;  % 景気が良い時の企業2への追加的利潤
    -0.15;  % 企業2の退出のためのコスト
    -2.10 ]; % 企業2の出店のためのコスト

% Profitを作成
pi1_cf1 = pi1gen(CounterfactualParameterValues).*CCP1Adjuster;
pi2_cf1 = pi2gen(CounterfactualParameterValues).*CCP2Adjuster;

% MPEを解く
[CCP1cf, CCP2cf, V1cf, V2cf] = f_MPE(TransitionMat, pi1_cf1, pi2_cf1, beta);

% 均衡における状態の遷移行列を計算する。fPｂsigma.m を利用する。
fPsigma = fP(TransitionMat, CCP1cf(:,2), CCP2cf(:,2));

initialstate = [1; zeros(7,1)];

% 初期の状態をセットする。n1 = n2 = 0, z = G
transitionpathcf = zeros(NumSimPeriods, 8);

for tt = 1:NumSimPeriods
    
    if (tt == 1)
        transitionpathcf(tt, :) = initial_state'; 
    elseif tt >= 2 
        transitionpathcf(tt, :) = ( fPsigma'*transitionpathcf(tt-1,:)' )';
    end
end

% 参入状況をカウントする。
n1_cf = sum( transitionpathcf.*repmat( [0,0,1,1,0,0,1,1],NumSimPeriods,1 ), 2);
n2_cf = sum( transitionpathcf.*repmat( [0,1,0,1,0,1,0,1],NumSimPeriods,1 ), 2);


% プロットを並べるためのfigureを作成する
figure;

% subplot関数を使用して2つのプロットを並べる
subplot(1,2,1); % 1行2列のsubplotの1つ目
plot(1:NumSimPeriods, n1, '--', 1:NumSimPeriods, n1_cf, "-o");
legend('現状維持', '新ブランド展開', 'Location', 'south');
title('企業1の店舗存在確率');

subplot(1,2,2); % 1行2列のsubplotの2つ目
plot(1:NumSimPeriods, n2, '--', 1:NumSimPeriods, n2_cf, "-o");
legend('現状維持', '新ブランド展開', 'Location', 'south');
title('企業2の店舗存在確率');

saveas(gcf, 'ProbEntry.jpg', 'jpg');

% ベースラインと差別化戦略

% 均衡におけるValue functionを比較すれば十分でした。。。

diff_value_1 = V1cf - V1base;
display([V1base, V1cf, diff_value_1])

