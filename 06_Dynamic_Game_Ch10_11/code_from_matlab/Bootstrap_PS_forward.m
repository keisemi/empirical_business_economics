function output = Bootstrap_PS_forward(FakeData, beta, ...
    EVrandom, UNIrandom, InitialState, NumSimMarkets, NumSimulations, NumSimPeriods)

% Step 1: Estimating CCP and transition probabilities
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

%NumSimPeriods  = 100;
%NumSimulations = 1000;
%InitialState = repmat( [1:8]', 1, 2);
%NumSimMarkets = 8;

% Step 2-1: Generating epsilons from etreme value distributions for
% choices and epsilons from uniform distribution for state
% transition


% Step 2-2: Generating V for any sigma
% Given EstimatedCCPs (and parameter values), calculate a
% discounted sum of profit for each firm and basis functions

% True CCP and Transition
[W1star, W2star] = VsigmaGeneration(EstimatedCCP1, EstimatedCCP2, EstimatedTransition, ...
    EVrandom, UNIrandom, InitialState, NumSimMarkets, NumSimulations, NumSimPeriods, beta);

obj_forward_PSD = @(x) Estimation_forward_PSD( [x(1); x(3:4);0; x(5);x(2); x(3:4);0;x(5)], ...
    W1star, W2star, EstimatedTransition, EstimatedCCP1, EstimatedCCP2, beta);

initial = [0.3, 0.2, -0.27, 0.45, -2.1]';

[opt, val]  = fminsearch(obj_forward_PSD, initial);

% Output
output = struct;
output.PSobj = val;
%output.estimatedCCP = [EstimatedCCP1, EstimatedCCP2];
%output.estimatedtransition = EstimatedTransition;
output.estimatedpayoffparam = opt;

end