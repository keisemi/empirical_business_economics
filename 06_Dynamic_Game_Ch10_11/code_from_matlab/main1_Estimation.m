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
% main1_Estimation.m: パラメタの推定（本ファイル）
% main2_Policy_Simulation.m: 反実仮想シミュレーション

% 本ファイルの流れは以下のようになっています
% 1. 下準備
% 2. 疑似データの生成 (第10回の内容に相当)
% 3. 推定 Step 1: CCPとTransition
% 4. 推定 Step 2-1: 推定したCCPの元でのForward Simulation
% 5. Forward Simulationを用いたPesendorer and Schmidt-Dengler (以下P-SD) 
% 6. BBLの不等式推定量
% 7. BBLの不等式推定量におけるBootstrap


%% 1. 下準備 %%

% データをクリアする
clear

% 下準備のコード。第11回と同じ。
run sub_1_prepare.m

%% 2. 疑似データの生成 %%

% 疑似データ作成のコード。第11回と同じ。
run sub_2_DGP.m

% FakeDataを保存する。
writematrix(FakeData, "FakeData_Matlab.csv")

%% 3. 推定 Step 1: CCPとTransition

% CCPと状態変数zの遷移確率を推定する。第11回と同じだが再掲。

EstimatedCCP1 = zeros(8,1);
EstimatedCCP2 = zeros(8,1);

for s=1:8
    SubData = FakeData(FakeData(:,3) == s,:);
    Subseta1iszero = SubData(SubData(:,7)==0,:);
    Subseta2iszero = SubData(SubData(:,8)==0,:);
    EstimatedCCP1(s,1) = size(Subseta1iszero,1)/size(SubData,1);
    EstimatedCCP2(s,1) = size(Subseta2iszero,1)/size(SubData,1);
end

EstimatedTransition = zeros(2,2);

DataL1 = [ 0 ; FakeData(1:end-1,4) ];
SubData = [FakeData DataL1];
SubData = SubData(SubData(:,2)~=1,:);

for z=1:2
    SubDataZ     = SubData(SubData(:,4) == z,:);
    SubDataZnext = SubDataZ(SubDataZ(:,9) == z,:);
    EstimatedTransition(z,z) = size(SubDataZnext,1)/size(SubDataZ,1);
end

EstimatedTransition(1,2) = 1-EstimatedTransition(1,1);
EstimatedTransition(2,1) = 1-EstimatedTransition(2,2);

% Check
disp("推定されたCCP")
disp([ EstimatedCCP1, EstimatedCCP2 ])

disp("推定されたzの遷移確率行列")
disp(EstimatedTransition)

%% 4. 推定 Step 2-1: 推定したCCPの元でのForward Simulation

% Note: State index and its corresponding values
% 1: G(1) 0 0 
% 2: G(1) 0 1  
% 3: G(1) 1 0
% 4: G(1) 1 1 
% 5: B(2) 0 0 
% 6: B(2) 0 1  
% 7: B(2) 1 0
% 8: B(2) 1 1 

% Forward simulationのパラメタ設定
NumSimPeriods  = 100;
NumSimFirms    = 2;
NumSimulations = 1000;

% Forward simulationにおけるInitialのStateを設定。
% 今回は状態変数が取りうる値が８通りなので、その８通りをInitialとする。
%InitialState = FakeData(FakeData(:,2) == 1,[1,3]);
InitialState = repmat( [1:8]', 1, 2);
NumSimMarkets = 8;

% Forward simulationで用いるランダムショックを準備する。

% Idiosyncratic shock (ロジットショック)のドローを準備する。
% F\left(\varepsilon_{n j}\right)=e^{-e^{-\varepsilon_{n j}}}
% F(x) = exp(-exp(-x)) を使う。
% -log(-log(F(x)) ) = x となるので、F(x)の部分を0-1で乱数ドローして、そのInversionを取る。
EVrandom     = -log(-log(rand([NumSimMarkets,NumSimPeriods,NumSimFirms,NumSimulations,8,3])));

% 景気のtransitionのシミュレーションに用いるドローを用意する。
UNIrandom = rand(NumSimMarkets,NumSimPeriods,NumSimulations);

% Sanity check として、DGPにおけるExante value functionと、Forward
% SimulationしたValueが一致するかをチェックする。

% DGPにおける CCP と Transition の元でForward Simulationする。
disp('Forward simulation1回にかかる時間')
tic
[W1star, W2star] = VsigmaGeneration(CCP1UpdatedMat(:,2), CCP2UpdatedMat(:,2), TransitionMat, ...
    EVrandom, UNIrandom, InitialState, NumSimMarkets, NumSimulations, NumSimPeriods, beta);
toc

param1 = [ TrueParameterValues(1:5) ; 1 ];
param2 = [ TrueParameterValues(6:10); 1 ];

display("True value, Simulated value, and their difference")
fa1 = [ExanteV1Updated, W1star'*param1, ExanteV1Updated - W1star'*param1]
fa2 = [ExanteV2Updated, W2star'*param2, ExanteV2Updated - W2star'*param2]

% Check value after considering normalization a.la. Agg-Suzuki
Normalized_TrueParam =  [ 0.3 - ((1-beta)/beta)*(-0.15) ;  % base for 1  (fixed profit)
    -0.27  ;  % number of rival's shop
    0.45  ;  % good economic condition
    -0.15 - (-0.15)  ;  % Closing costs
    -2.1 + (-0.15) ;  % Opening costs
    0.2 - ((1-beta)/beta)*(-0.15)  ;  % base for 2 (fixed profit)
    -0.27  ;  % number of rival's shop
    0.45  ;  % good economic condition
    -0.15 - (-0.15)  ;  % Closing costs
    -2.1 + (-0.15) ]; % Opening costs

normparam1 = [ Normalized_TrueParam(1:5) ; 1 ];
normparam2 = [ Normalized_TrueParam(6:10); 1 ];

disp("True value, simulated value (original parameter), simulated value (noramlized parameter)")
fafa1 = [ExanteV1Updated, W1star'*param1, W1star'*normparam1]
fafa2 = [ExanteV2Updated, W2star'*param2, W2star'*normparam2]

disp("Difference b.w. two simulated values (with and without normalization")
W1star'*param1 - W1star'*normparam1
W2star'*param2 - W2star'*normparam2

% 著者メモ：各企業が店舗を開いている(n=1)のStateにおいて、
% -0.1875 * 0.8 (discount factor) = - 0.15 
% という差が生じている。
% これは、Closing costをゼロに基準化していることに起因する。

% 推定したCCPのもとで、Forward Simulationを行い、Value functionの基底を出す。
[W1star, W2star] = VsigmaGeneration(EstimatedCCP1, EstimatedCCP2, EstimatedTransition,...
    EVrandom, UNIrandom, InitialState, NumSimMarkets, NumSimulations, NumSimPeriods, beta);

%% 5. P-SD with Forward Simulation

% Forward simulationしたValueを用いて、P-SD流に推定を行う。

% まず目的関数の定義
obj_forward_PSD = @(x) Estimation_forward_PSD( [x(1); x(3:4);0; x(5);x(2); x(3:4);0;x(5)], ...
    W1star, W2star, EstimatedTransition, EstimatedCCP1, EstimatedCCP2, beta);

% 初期値はTrueにしておく。
initial = [0.3375, 0.2375, -0.27, 0.45, -2.25]';

% 推定
[opt, val]  = fminsearch(obj_forward_PSD, initial);

% Forward simulationを用いたP-SDについてBootstrapを行う。
% 以下のスクリプトを実行する。繰り返し計算を行うため多少時間かかる。

run sub_5_Bootstrap_PSD_forward.m

% Rのために保存
save random_number_matlab_PSD.mat EVrandom UNIrandom bootindex


%% 6. Estimation by BBL Inequality

% 計算速度の観点から、以下ではSimulationの期間を30期間とする。
NumSimPeriods  = 30;

% 乱数の発生が必要なため、乱数発生のためのシードを設定する
rng(2021)

% Idiosyncratic shock (ロジットショック)のドローを準備する。
% F\left(\varepsilon_{n j}\right)=e^{-e^{-\varepsilon_{n j}}}
% F(x) = exp(-exp(-x)) を使う。
% -log(-log(F(x)) ) = x となるので、F(x)の部分を0-1で乱数ドローして、そのInversionを取る。
EVrandom     = -log(-log(rand([NumSimMarkets,NumSimPeriods,NumSimFirms,NumSimulations,8,3])));

% 景気のtransitionのシミュレーションに用いるドローを用意する。
UNIrandom = rand(NumSimMarkets,NumSimPeriods,NumSimulations);

% 推定したCCPのもとで、Forward Simulationを行い、Value functionの基底を出す。
[W1star, W2star] = VsigmaGeneration(EstimatedCCP1, EstimatedCCP2, EstimatedTransition,...
    EVrandom, UNIrandom, InitialState, NumSimMarkets, NumSimulations, NumSimPeriods, beta);

% 乱数の発生が必要なため、乱数発生のためのシードを設定する
rng(2022)

% CCPのPerturbationを行う。
NumPerturbations = 200;
PerturbedCCP1    = repmat(EstimatedCCP1, 1, NumPerturbations) + 0.1*randn(8,NumPerturbations);
PerturbedCCP2    = repmat(EstimatedCCP2, 1, NumPerturbations) + 0.1*randn(8,NumPerturbations);

% To make CCP be inside of [0,1] or [0.001, 0.999]
for i=1:8
    for j = 1:NumPerturbations
        PerturbedCCP1(i,j) = max(PerturbedCCP1(i,j), 0.001);
        PerturbedCCP1(i,j) = min(PerturbedCCP1(i,j), 0.999);

        PerturbedCCP2(i,j) = max(PerturbedCCP2(i,j), 0.001);
        PerturbedCCP2(i,j) = min(PerturbedCCP2(i,j), 0.999);
    end
end

W1_all = zeros(6,NumSimMarkets,NumPerturbations);
W2_all = zeros(6,NumSimMarkets,NumPerturbations);

% PerturbationしたCCPについてForward simulationを行い、Valueを計算する。
for per = 1:NumPerturbations
    tic
    per
    [ W1_p, ~ ]  = VsigmaGeneration(PerturbedCCP1(:,per), EstimatedCCP2, EstimatedTransition, ...
        EVrandom, UNIrandom, InitialState, NumSimMarkets, NumSimulations, NumSimPeriods, beta);
    W1_all(:,:,per) = W1_p;

    [ ~, W2_p ]  = VsigmaGeneration(EstimatedCCP1, PerturbedCCP2(:,per), EstimatedTransition,...
        EVrandom, UNIrandom, InitialState, NumSimMarkets, NumSimulations, NumSimPeriods, beta);
    W2_all(:,:,per) = W2_p;
    toc 
end

% SimulationしたValueに基づいて目的関数を定義し、パラメタの推定を行う。

% Set initial value
initial = [0.3, 0.2, -0.27, 0.45, -2.1]';

% 推定：Nonlinear least squaresとして推定する。
obj_forward_BBL_NLS = @(x) BBLobjective_NLS( [x(1); x(3:4);0; x(5);x(2); x(3:4);0;x(5)], ...
    NumPerturbations, W1star, W2star, W1_all, W2_all);

options=optimset('Display',     'off',...
    'TolFun',1e-10,...
     'TolX',       1e-10,...
    'MaxIter',     10000,...
    'MaxFunEvals', 10000);

initial  = [0.1,0.1,0.1,0.1,0.1]';
[opt0, val0]  = lsqnonlin(obj_forward_BBL_NLS, initial, [], [], options);
display(val0)
display(opt0')

% 参考：fminserachでやっても同じ結果が得られる。
obj_forward_BBL = @(x) BBLobjective( [x(1); x(3:4);0; x(5);x(2); x(3:4);0;x(5)], ...
    NumSimMarkets, NumPerturbations, W1star, W2star, W1_all, W2_all);

[opt1, val1]  = fminsearch(obj_forward_BBL, initial);
[val1, opt1']

% Rのために保存
save random_number_matlab_BBL.mat EVrandom UNIrandom PerturbedCCP1 PerturbedCCP2

%% 7. BBL (inequality estimator)によるBootstrap

% Pertubationの回数
NumPerturbations = 200;

% Bootstrapのための乱数を固定
rng(2023)

% Bootstrap のリサンプリング回数
numBootSample = 100;

% 各Bootstrap sampleで用いるマーケットのインデックスを乱数から発生させる。
% 市場が500個であるため、１から５００の整数について、重複を許して５００個ドローする。
bootindex = randi(500, 500, numBootSample);

% 結果を保存するための行列
bootresult_payoff = zeros(5,numBootSample);


% Bootstrap を実行するループ

% 注意：：Bootstrapには非常に長い時間がかかるため注意。
% 計算時間＝(BBLの推定一回にかかる時間)*(Bootstrap回数)
% 執筆者のPCではBBLの推定１回に６分ー７分程度かかったため、
% 合計で、おおよそ１０時間－１２時間程度かかると思われる。

% なお、この部分の結果(変数bootresult_payoffに保存)はresult_BBL_bootstrap.matとして利用可能である。

for b = 1:numBootSample

    tic
    disp(append("Boostrap: ", num2str(b)))
    % Bootstrap sampleを構築する。
    bootsample = zeros(500*50, 8);
    for m = 1:500
        temp = FakeData( FakeData(:,1) == bootindex(m,b) ,:);
        bootsample( 1+50*(m-1):50*m  ,:) = temp;
    end

    % 構築したBootstrap sampleを用いて推定を行う。
    output = Bootstrap_BBL(bootsample, beta, ...
        NumSimMarkets, NumSimulations, NumSimPeriods, NumPerturbations);

    % 推定結果を保存する。
    bootresult_payoff(:,b) = output.estimatedpayoffparam;

    toc

end

% 結果を保存。
save result_BBL_bootstrap.mat bootresult_payoff


% 結果のまとめ
true = [0.3, 0.2, -0.27, 0.45, -2.1]';
disp("Payoff parameter: True, Normalized true, Estimated, SE ")
normalized_param = [ 0.3 - (1-beta)/beta*(-0.15), 0.2 - (1-beta)/beta*(-0.15),-0.27, 0.45, -2.1 + (-0.15)];
disp([ true, normalized_param',  std(bootresult_payoff')' ])

