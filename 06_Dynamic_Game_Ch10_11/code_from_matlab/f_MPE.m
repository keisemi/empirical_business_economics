function [outCCP1, outCCP2, outV1, outV2] = f_MPE(TransitionMat, pi1, pi2, beta)


% Step 1: CCPの初期値の設定（および設定したCCPをベクトルから行列へ変換する）

% 以下ではCCPを全て0.5とする
CCP1 = 0.5*ones(8,1);
CCP2 = 0.5*ones(8,1);

% 上で与えた初期値のベクトルを行列へと変換する。例えば、CCP1Matは
% CCP1Mat = [ 0   0.5 0.5
%             0   0.5 0.5
%             :    :   :
%            0.5  0.5  0
%            0.5  0.5  0  ] のような行列になる
CCP1Mat = CCP1Transform(CCP1);
CCP2Mat = CCP2Transform(CCP2);

% NumMaxChoices
NumMaxChoices = 3;

% CCP adjuster
CCP1Adjuster = [ 0 1 1 ; 0 1 1 ; 1 1 0 ; 1 1 0 ; ...
    0 1 1 ; 0 1 1 ; 1 1 0 ; 1 1 0 ];
CCP2Adjuster = [ 0 1 1 ; 1 1 0 ; 0 1 1 ; 1 1 0 ; ...
    0 1 1 ; 1 1 0 ; 0 1 1 ; 1 1 0 ];



% Step 2: Step 1 で与えられたCCPの初期値を基に、事前の価値関数を計算する

% 基本的には（４）式の計算を行いたいため、以下で必要なパーツを計算する

% まず、F^P^sigmaを求める。これは109ページの囲みで与えられているように、
% 遷移行列とCCP1とCCP2から作られる３つの行列の、要素ごとの積（アダマール
% 積）として計算できる
fPsigma = fP(TransitionMat, CCP1, CCP2);

% パラメターの下で、pi1Psigma と pi2Psigma を計算する
pi1Psigma = pi1PsigmaGen(pi1, CCP2Mat, NumMaxChoices);
pi2Psigma = pi2PsigmaGen(pi2, CCP1Mat, NumMaxChoices);

% e^P_i(a_i,s)を計算する
eP1 = 0.5772*ones(8,3) - CCP1LogTransform(CCP1);
eP2 = 0.5772*ones(8,3) - CCP2LogTransform(CCP2);

% 上で求めたものを（４）式に代入して、事前の価値関数を求める
ExanteV1 = (eye(8)-beta*fPsigma)\sum(CCP1Mat.*(pi1Psigma+eP1),2);
ExanteV2 = (eye(8)-beta*fPsigma)\sum(CCP2Mat.*(pi2Psigma+eP2),2);
% (pij + ePj) が何らかの値をとっていたとしても、CCP1MatやCCP2Matを
% かけることで、取ることのできない選択肢にはゼロが振られることになり、
% 分母の和は正しく計算できていることがわかる。例えば、
% CCP1Mat.*(pi1Psigma+eP1) を表示させてみると、実際にいくつかの
% 要素がゼロになっていることが確認できる。


% Step 3: 事前の価値関数を基に、CCPを計算する

% 各企業の行動が与えられたときの遷移行列を求める。これは、本文では
% f^P*(s'|s,a_i)として表現されている行列である。また、fP_a1isM1 は
% a_1 = -1 (a_1 is equal to minus 1), fP_a1is0 は a_1 = 0 (a_1 is
% equal to zero), fP_a1isP1 は a_1 = 1 (a_1 is equal to 1) の時の
% 遷移行列を表している。
[fP_a1isM1, fP_a1is0, fP_a1isP1] = fP_a1given(TransitionMat, CCP2);
[fP_a2isM1, fP_a2is0, fP_a2isP1] = fP_a2given(TransitionMat, CCP1);

% 以下では、求められた事前の価値関数をもとにCCPを更新する。最終的には、
% 計算の簡便化のため更新されたCCPを行列の形（CCP1UpdatedMat）とベクトル
% の形（CCP1Updated）の両方で求めておく。この計算は一見複雑だが、需要
% 関数の推定の時と似た計算であるため、説明は省略する。最初のブロックは企
% 業1について、次のブロックは企業2について更新されたCCPを求めている。
NewSigmaSeed1  =  (pi1Psigma + beta*[ fP_a1isM1*ExanteV1 fP_a1is0*ExanteV1 fP_a1isP1*ExanteV1 ]).*CCP1Adjuster;
NewSigmaDeno1  = sum(exp(NewSigmaSeed1), 2) - ones(8,1);
NewSigmaDeno1M = repmat(NewSigmaDeno1,1,3);
NewSigma1      = exp(NewSigmaSeed1)./NewSigmaDeno1M;
CCP1UpdatedMat = NewSigma1.*CCP1Adjuster;
CCP1Updated    = CCP1UpdatedMat(:,2);

NewSigmaSeed2 = (pi2Psigma + beta*[ fP_a2isM1*ExanteV2 fP_a2is0*ExanteV2 fP_a2isP1*ExanteV2 ]).*CCP2Adjuster;
NewSigmaDeno2 = sum(exp(NewSigmaSeed2), 2) - ones(8,1);
NewSigmaDeno2M = repmat(NewSigmaDeno2,1,3);
NewSigma2     = exp(NewSigmaSeed2)./NewSigmaDeno2M;
CCP2UpdatedMat = NewSigma2.*CCP2Adjuster;
CCP2Updated   = CCP2UpdatedMat(:,2);


% Step 4: アップデートされたCCPを基に、再び事前の価値関数を計算する。ただし、
% 　　　　Step 2と同一の作業であるため、説明は省略する。

fPsigma = fP(TransitionMat, CCP1Updated, CCP2Updated);

pi1Psigma = pi1PsigmaGen(pi1, CCP2UpdatedMat, NumMaxChoices);
pi2Psigma = pi2PsigmaGen(pi2, CCP1UpdatedMat, NumMaxChoices);

eP1 = 0.5772*ones(8,3) - CCP1LogTransform(CCP1Updated);
eP2 = 0.5772*ones(8,3) - CCP2LogTransform(CCP2Updated);

ExanteV1Updated = (eye(8)-beta*fPsigma)\sum(CCP1UpdatedMat.*(pi1Psigma+eP1),2);
ExanteV2Updated = (eye(8)-beta*fPsigma)\sum(CCP2UpdatedMat.*(pi2Psigma+eP2),2);


% Step 5: Step 3 と Step 4 を事前の価値関数が一致するまで繰り返す

% 前の事前の価値関数と更新された事前の価値関数の差を定義する
DiffExanteV = sum((ExanteV1Updated - ExanteV1).^2 + ...
    (ExanteV2Updated - ExanteV2).^2 );

% 上で定義した差が1.0e-12よりも小さくなるまで、繰り返す
while ( DiffExanteV > 1.0e-12 )

    % 更新されたCCPと事前の価値関数を、CCPとExanteVとして置き換える
    CCP1 = CCP1Updated;
    CCP2 = CCP2Updated;

    ExanteV1 = ExanteV1Updated;
    ExanteV2 = ExanteV2Updated;

    % Step 3 を再度実行する（上と全く一緒のため、説明は省略する）
    [fP_a1isM1, fP_a1is0, fP_a1isP1] = fP_a1given(TransitionMat, CCP2);
    [fP_a2isM1, fP_a2is0, fP_a2isP1] = fP_a2given(TransitionMat, CCP1);

    NewSigmaSeed1 =  (pi1Psigma + beta*[ fP_a1isM1*ExanteV1 fP_a1is0*ExanteV1 fP_a1isP1*ExanteV1 ]).*CCP1Adjuster;
    NewSigmaDeno1 = sum(exp(NewSigmaSeed1), 2) - ones(8,1);
    NewSigmaDeno1M = repmat(NewSigmaDeno1,1,3);
    NewSigma1     = exp(NewSigmaSeed1)./NewSigmaDeno1M;
    CCP1UpdatedMat = NewSigma1.*CCP1Adjuster;
    CCP1Updated   = CCP1UpdatedMat(:,2);

    NewSigmaSeed2 = (pi2Psigma + beta*[ fP_a2isM1*ExanteV2 fP_a2is0*ExanteV2 fP_a2isP1*ExanteV2 ]).*CCP2Adjuster;
    NewSigmaDeno2 = sum(exp(NewSigmaSeed2), 2) - ones(8,1);
    NewSigmaDeno2M = repmat(NewSigmaDeno2,1,3);
    NewSigma2     = exp(NewSigmaSeed2)./NewSigmaDeno2M;
    CCP2UpdatedMat = NewSigma2.*CCP2Adjuster;
    CCP2Updated   = CCP2UpdatedMat(:,2);

    % Step 4 を再度実行する（上と全く一緒のため、説明は省略する）
    fPsigma = fP(TransitionMat, CCP1Updated, CCP2Updated);

    pi1Psigma = pi1PsigmaGen(pi1, CCP2UpdatedMat, NumMaxChoices);
    pi2Psigma = pi2PsigmaGen(pi2, CCP1UpdatedMat, NumMaxChoices);

    eP1 = double(eulergamma)*ones(8,3) - CCP1LogTransform(CCP1Updated);
    eP2 = double(eulergamma)*ones(8,3) - CCP2LogTransform(CCP2Updated);

    ExanteV1Updated = (eye(8)-beta*fPsigma)\sum(CCP1UpdatedMat.*(pi1Psigma+eP1),2);
    ExanteV2Updated = (eye(8)-beta*fPsigma)\sum(CCP2UpdatedMat.*(pi2Psigma+eP2),2);

    % Step 5 の冒頭部分の差分の計算を再度実行する
    DiffExanteV = sum((ExanteV1Updated - ExanteV1).^2 + ...
        (ExanteV2Updated - ExanteV2).^2 );
end

outCCP1 = CCP1UpdatedMat;
outCCP2 = CCP2UpdatedMat;
outV1 = ExanteV1Updated;
outV2 = ExanteV2Updated;

end





