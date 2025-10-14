function output = CCP_to_Value_to_prediction(TrueParameterValues, CCP1, CCP2, TransitionProb,beta)

% CCP から ValueをInversionし、その上でCCPを予測するコード。
% Step 1: Setup
% Step 2: CCPを所与としたときの、ValueのInversion
% Step 3: Value及びパラメタを元に、CCPを予測する。

% Step 1: Setup
% Transition matrix (Good and Bad)
TransitionMat = TransitionProb;

NumMaxChoices = 3;

% CCP1Adjuster and CCP2Adjuster are matrices where available
% choices are denoted by 1 and 0 otherwise
CCP1Adjuster = [ 0 1 1 ; 0 1 1 ; 1 1 0 ; 1 1 0 ; ...
    0 1 1 ; 0 1 1 ; 1 1 0 ; 1 1 0 ];
CCP2Adjuster = [ 0 1 1 ; 1 1 0 ; 0 1 1 ; 1 1 0 ; ...
    0 1 1 ; 1 1 0 ; 0 1 1 ; 1 1 0 ];


% Profit pi1 and pi2
pi1 = pi1gen(TrueParameterValues).*CCP1Adjuster;
pi2 = pi2gen(TrueParameterValues).*CCP2Adjuster;


% Transforming CCP vectors into matrices
CCP1Mat = CCP1Transform(CCP1);
CCP2Mat = CCP2Transform(CCP2);

% Step 2: Given initial values of CCP, calculate Exante Value Function

% Calcuating f^P, given a transition matrix, CCP1Matrix, and
% CCP2Matrix. Type "sum(fPsigma,2)" to check whether a summuation
% of each row is indeed equal to 1, if needed.
fPsigma = fP(TransitionMat, CCP1, CCP2);

% Given parameter values, calculate pi1Psigma and pi2Psigma
pi1Psigma = pi1PsigmaGen(pi1, CCP2Mat, NumMaxChoices);
pi2Psigma = pi2PsigmaGen(pi2, CCP1Mat, NumMaxChoices);

% Calculating e^P_i(a,s) -- an 8*1 vector
eP1 = double(eulergamma)*ones(8,3) - CCP1LogTransform(CCP1);
eP2 = double(eulergamma)*ones(8,3) - CCP2LogTransform(CCP2);
%eP1 = - CCP1LogTransform(CCP1);
%eP2 = - CCP2LogTransform(CCP2);

% Obtaining exante value functions through an analytical formula
ExanteV1 = (eye(8)-beta*fPsigma)\sum(CCP1Mat.*(pi1Psigma+eP1),2);
ExanteV2 = (eye(8)-beta*fPsigma)\sum(CCP2Mat.*(pi2Psigma+eP2),2);


% Step 3: Caclulating CCPs based on ExanteVF1 and ExanteVF2

% Calculating fP given a_1 or a_2
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
output = [CCP1Updated, CCP2Updated];
output = [ExanteV1, ExanteV2];


end