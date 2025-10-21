function output = pi1PsigmaGen(pi1, Mat2, ChoiceMaxNum)

pi1_dec = pi1(:,1);
pi1_0   = pi1(:,2);
pi1_inc = pi1(:,3);

pi1m1 = (pi1_dec.*Mat2)*ones(ChoiceMaxNum,1);
pi1_0 = (  pi1_0.*Mat2)*ones(ChoiceMaxNum,1);
pi1i1 = (pi1_inc.*Mat2)*ones(ChoiceMaxNum,1);

output = [ pi1m1 pi1_0 pi1i1 ]; 