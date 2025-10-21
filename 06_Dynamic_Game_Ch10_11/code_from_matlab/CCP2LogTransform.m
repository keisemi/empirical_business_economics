function output = CCP2LogTransform(x) 

% x is an 8*1 vector in the example and we make it an 8*3 matrix 
output =    [   0           log(x(1))    log(1-x(1)); 
                log(1-x(2)) log(x(2))    0     ; 
                0           log(x(3))    log(1-x(3)); 
                log(1-x(4)) log(x(4))    0     ;
                0           log(x(5))    log(1-x(5)); 
                log(1-x(6)) log(x(6))    0     ; 
                0           log(x(7))    log(1-x(7)); 
                log(1-x(8)) log(x(8))    0          ];
