function objvalue = BBLobjective(theta, NumSimMarkets, NumPerturbations, W1star, W2star, W1_all, W2_all)


Param1 = [ theta(1:5) ; 1 ]; 
Param2 = [ theta(6:10); 1 ]; 

objvalue = 0; 

for per = 1:NumPerturbations 

    temp1 = min(W1star'*Param1 - W1_all(:,:,per)'*Param1, 0); 
    temp2 = min(W2star'*Param2 - W2_all(:,:,per)'*Param2, 0);

    objvalue = objvalue + temp1'*temp1 + temp2'*temp2; 

end 

objvalue = objvalue;

%/(2*NumSimMarkets*NumPerturbations);
%display(objvalue)
