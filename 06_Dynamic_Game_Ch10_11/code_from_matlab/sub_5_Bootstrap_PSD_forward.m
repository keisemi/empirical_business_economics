%% 5. P-SD + Forward SimulationにおけるBootstrap

% マーケット単位でリサンプリングを行う。マーケットは500個

% Bootstrapのための乱数を固定
rng(2023)

% Bootstrap のリサンプリング回数
numBootSample = 100;

% 各Bootstrap sampleで用いるマーケットのインデックスを乱数から発生させる。
% 市場が500個であるため、１から５００の整数について、重複を許して５００個ドローする。
bootindex = randi(500, 500, numBootSample);

% 結果を保存するための行列
bootresult_transition = zeros(2, numBootSample);
bootresult_CCP1 = zeros(8, numBootSample );
bootresult_CCP2 = zeros(8, numBootSample );
bootresult_payoff = zeros(5,numBootSample);

% Bootstrap を実行するループ
for b = 1:numBootSample

    disp(append("Boostrap: ", num2str(b)))
    tic
    % Bootstrap sampleを構築する。
    bootsample = zeros(500*50, 8);
    for m = 1:500
        temp = FakeData( FakeData(:,1) == bootindex(m,b) ,:);
        bootsample( 1+50*(m-1):50*m  ,:) = temp;
    end

    % 構築したBootstrap sampleを用いて推定を行う。
    output = Bootstrap_PS_forward(bootsample, beta, ...
    EVrandom, UNIrandom, InitialState, NumSimMarkets, NumSimulations, NumSimPeriods);
    
    % 推定結果を保存する。
    bootresult_payoff(:,b) = output.estimatedpayoffparam;
    %bootresult_CCP1(:,b) = output.estimatedCCP(:,1);
    %bootresult_CCP2(:,b) = output.estimatedCCP(:,2);
    %bootresult_transition(:,b) = diag(output.estimatedtransition);

    toc

end

% Payoff parameter
true = [0.3, 0.2, -0.27, 0.45, -2.1]';
disp("Payoff parameter: True, Normalized true, Estimated, SE ")
normalized_param = [ 0.3 - (1-beta)/beta*(-0.15), 0.2 - (1-beta)/beta*(-0.15),-0.27, 0.45, -2.1 + (-0.15)];
disp([ true, normalized_param', opt, std(bootresult_payoff')' ])
