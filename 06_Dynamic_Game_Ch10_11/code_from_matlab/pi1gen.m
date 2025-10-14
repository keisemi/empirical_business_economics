function output = pi1gen(theta) 


base = [0;
        0;
        theta(1) + theta(3) ;
        theta(1) + theta(2) + theta(3) ;
        0;
        0;
        theta(1) ;
        theta(1) + theta(2)]; 

invdiv = [ theta(4) 0 theta(5)]; 

output = repmat(base,1,3) + repmat(invdiv, 8, 1); 

