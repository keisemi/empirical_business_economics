function output = CCP1Transform(x) 

output = [  0           x(1) 1-x(1); 
            0           x(2) 1-x(2); 
            1-x(3) x(3) 0          ; 
            1-x(4) x(4) 0          ;
            0           x(5) 1-x(5); 
            0           x(6) 1-x(6); 
            1-x(7) x(7) 0          ; 
            1-x(8) x(8) 0          ]; 

end 
