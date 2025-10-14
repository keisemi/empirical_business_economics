function output = pi2PsigmaGen(pi2, Mat1, ChoiceMaxNum)


pi2_dec = pi2(:,1);
pi2_0   = pi2(:,2);
pi2_inc = pi2(:,3);


pi2m1 = (pi2_dec.*Mat1)*ones(ChoiceMaxNum,1);
pi2_0 = (  pi2_0.*Mat1)*ones(ChoiceMaxNum,1);
pi2i1 = (pi2_inc.*Mat1)*ones(ChoiceMaxNum,1);


output = [ pi2m1 pi2_0 pi2i1 ]; 