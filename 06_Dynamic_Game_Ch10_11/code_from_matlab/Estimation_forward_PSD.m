function obj = Estimation_forward_PSD(param, W1star, W2star, TransitionMat, CCP1, CCP2, beta)

NumMaxChoices = 3;
    
% Forward simuからValueを作る。
param1 = [ param(1:5) ; 1 ];
param2 = [ param(6:10); 1 ];

ExanteV1 = W1star'*param1;
ExanteV2 = W2star'*param2;

% CCP1Adjuster and CCP2Adjuster are matrices where available
% choices are denoted by 1 and 0 otherwise
CCP1Adjuster = [ 0 1 1 ; 0 1 1 ; 1 1 0 ; 1 1 0 ; ...
    0 1 1 ; 0 1 1 ; 1 1 0 ; 1 1 0 ];
CCP2Adjuster = [ 0 1 1 ; 1 1 0 ; 0 1 1 ; 1 1 0 ; ...
    0 1 1 ; 1 1 0 ; 0 1 1 ; 1 1 0 ];


% Profit pi1 and pi2
pi1 = pi1gen(param).*CCP1Adjuster;
pi2 = pi2gen(param).*CCP2Adjuster;


% Transforming CCP vectors into matrices
CCP1Mat = CCP1Transform(CCP1);
CCP2Mat = CCP2Transform(CCP2);


% Given parameter values, calculate pi1Psigma and pi2Psigma
pi1Psigma = pi1PsigmaGen(pi1, CCP2Mat, NumMaxChoices);
pi2Psigma = pi2PsigmaGen(pi2, CCP1Mat, NumMaxChoices);



[fP_a1isM1, fP_a1is0, fP_a1isP1] = fP_a1given(TransitionMat, CCP2);
[fP_a2isM1, fP_a2is0, fP_a2isP1] = fP_a2given(TransitionMat, CCP1);

% Given ExanteVF1 and ExanteVF2, calculate CCPs again
NewSigmaSeed1  =  (pi1Psigma + beta*[ fP_a1isM1*ExanteV1 fP_a1is0*ExanteV1 fP_a1isP1*ExanteV1 ]).*CCP1Adjuster;
NewSigmaDeno1  = sum(exp(NewSigmaSeed1), 2) - ones(8,1);
NewSigmaDeno1M  = repmat(NewSigmaDeno1,1,3);
NewSigma1      = exp(NewSigmaSeed1)./NewSigmaDeno1M;
CCP1UpdatedMat = NewSigma1.*CCP1Adjuster;
CCP1Updated    = CCP1UpdatedMat(:,2);

NewSigmaSeed2 = (pi2Psigma + beta*[ fP_a2isM1*ExanteV2 fP_a2is0*ExanteV2 fP_a2isP1*ExanteV2 ]).*CCP2Adjuster;
NewSigmaDeno2 = sum(exp(NewSigmaSeed2), 2) - ones(8,1);
NewSigmaDeno2M = repmat(NewSigmaDeno2,1,3);
NewSigma2     = exp(NewSigmaSeed2)./NewSigmaDeno2M;
CCP2UpdatedMat = NewSigma2.*CCP2Adjuster;
CCP2Updated   = CCP2UpdatedMat(:,2);

% Output
CCP_updated = [CCP1Updated, CCP2Updated];


% Value function

obj = sum ((CCP_updated(:,1)-CCP1 ).^2) + sum ((CCP_updated(:,2)-CCP2 ).^2);


end