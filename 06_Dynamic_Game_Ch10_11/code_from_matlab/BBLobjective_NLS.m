function objvec = BBLobjective_NLS(theta, NumPerturbations, W1star, W2star, W1_all, W2_all)


Param1 = [ theta(1:5) ; 1 ]; 
Param2 = [ theta(6:10); 1 ]; 

fa1 = repmat(W1star, 1, 1, NumPerturbations);
fa2 = repmat(W2star, 1, 1, NumPerturbations);

diff1 = fa1 - W1_all;
diff2 = fa2 - W2_all; 

% reshapeする。行が各要素になる
temp1 = reshape(diff1, 6, 8*NumPerturbations);
temp2 = reshape(diff2, 6, 8*NumPerturbations);

val1 = temp1'*Param1;
val2 = temp2'*Param2;

val1 = min(val1, 0);
val2 = min(val2, 0);

% すべてStackして返す

objvec = [val1;val2];