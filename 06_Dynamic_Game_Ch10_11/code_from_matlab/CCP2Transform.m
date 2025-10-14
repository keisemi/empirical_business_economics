function output = CCP2Transform(x) 

% x is an 8*1 vector in the example and we make it an 8*3 matrix 
output =    [   0           x(1)    1-x(1); 
                1-x(2) x(2)    0     ; 
                0           x(3)    1-x(3); 
                1-x(4) x(4)    0     ;
                0           x(5)    1-x(5); 
                1-x(6) x(6)    0     ; 
                0           x(7)    1-x(7); 
                1-x(8) x(8)    0          ];


end
