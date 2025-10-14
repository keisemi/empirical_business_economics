function output = Bootstrap_BBL(FakeData, beta, ...
    NumSimMarkets, NumSimulations, NumSimPeriods, NumPerturbations)

% Step 1: Estimateing CCP and transition probabilities
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


% Step 2: P-SD estimator with forward simulation

% NumSimPeriods  = 100;
NumSimFirms    = 2;
% NumSimulations = 1000;

InitialState = repmat( [1:8]', 1, 2);
% NumSimMarkets = 8;

% Step 2-1: Generating epsilons from etreme value distributions for
% choices and epsilons from uniform distribution for state
% transition

% 乱数の発生が必要なため、乱数発生のためのシードを設定する
rng(2021)

% Idiosyncratic shock (ロジットショック)のドローを準備する。
% F\left(\varepsilon_{n j}\right)=e^{-e^{-\varepsilon_{n j}}}
% F(x) = exp(-exp(-x)) を使う。
% -log(-log(F(x)) ) = x となるので、F(x)の部分を0-1で乱数ドローして、そのInversionを取る。
EVrandom     = -log(-log(rand([NumSimMarkets,NumSimPeriods,NumSimFirms,NumSimulations,8,3])));

% 景気のtransitionのシミュレーションに用いるドローを用意する。
UNIrandom = rand(NumSimMarkets,NumSimPeriods,NumSimulations);

% Step 2-2: Generating V for any sigma
% Given EstimatedCCPs (and parameter values), calculate a
% discounted sum of profit for each firm and basis functions

% True CCP and Transition
[W1star, W2star] = VsigmaGeneration(EstimatedCCP1, EstimatedCCP2, EstimatedTransition,...
        EVrandom, UNIrandom, InitialState, NumSimMarkets, NumSimulations, NumSimPeriods, beta);

% Step 2-3:

% 乱数の発生が必要なため、乱数発生のためのシードを設定する
rng(2022)

% NumPerturbations = 200;
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

for per = 1:NumPerturbations
    %tic
    %per
    [ W1_p, ~ ]  = VsigmaGeneration(PerturbedCCP1(:,per), EstimatedCCP2, EstimatedTransition,...
        EVrandom, UNIrandom, InitialState, NumSimMarkets, NumSimulations, NumSimPeriods, beta);
    W1_all(:,:,per) = W1_p;

    [ ~, W2_p ]  = VsigmaGeneration(EstimatedCCP1, PerturbedCCP2(:,per), EstimatedTransition,...
        EVrandom, UNIrandom, InitialState, NumSimMarkets, NumSimulations, NumSimPeriods, beta);
    W2_all(:,:,per) = W2_p;
    %toc 
end

% Step 2-4: Mimimizing the objective function

% Set initial
initial  = [0.1,0.1,0.1,0.1,0.1]';
initial = [0.3, 0.2, -0.27, 0.45, -2.1]';

% Backupするべき値
%    0.3375
%    0.2375
%   -0.2700
%    0.4500
%   -2.2500

% Nonlinear least squaresでやる。

obj_forward_BBL_NLS = @(x) BBLobjective_NLS( [x(1); x(3:4);0; x(5);x(2); x(3:4);0;x(5)], ...
    NumPerturbations, W1star, W2star, W1_all, W2_all);

options=optimset('Display',     'off',...
    'TolFun',1e-10,...
     'TolX',       1e-10,...
    'MaxIter',     10000,...
    'MaxFunEvals', 10000);

initial  = [0.1,0.1,0.1,0.1,0.1]';
[opt0, val0]  = lsqnonlin(obj_forward_BBL_NLS, initial, [], [], options);

% Output
output = struct;
output.PSobj = val0;
%output.estimatedCCP = [EstimatedCCP1, EstimatedCCP2];
%output.estimatedtransition = EstimatedTransition;
output.estimatedpayoffparam = opt0;

end