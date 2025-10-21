function output = pi2gen(theta) 


base = [0; 
        theta(6) + theta(8) ;
        0; 
        theta(6) + theta(7) + theta(8) ;
        0;
        theta(6) ;
        0; 
        theta(6) + theta(7) ]; 

invdiv = [ theta(9) 0 theta(10)]; 

output = repmat(base,1,3) + repmat(invdiv, 8, 1); 

