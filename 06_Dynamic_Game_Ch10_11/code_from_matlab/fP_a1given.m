function [output1, output2, output3] = fP_a1given(Matrix1, Vec2)


% Exogenous state transition probabilities 
TempMat0 = kron(Matrix1,ones(4));

TempMat2 = [ Vec2(1)    1-Vec2(1)   Vec2(1)     1-Vec2(1)   Vec2(1)     1-Vec2(1)   Vec2(1)     1-Vec2(1) ; 
             1-Vec2(2)  Vec2(2)     1-Vec2(2)   Vec2(2)     1-Vec2(2)   Vec2(2)     1-Vec2(2)   Vec2(2) ;
             Vec2(3)    1-Vec2(3)   Vec2(3)     1-Vec2(3)   Vec2(3)     1-Vec2(3)   Vec2(3)     1-Vec2(3) ;
             1-Vec2(4)  Vec2(4)     1-Vec2(4)   Vec2(4)     1-Vec2(4)   Vec2(4)     1-Vec2(4)   Vec2(4) ;
             Vec2(5)    1-Vec2(5)   Vec2(5)     1-Vec2(5)   Vec2(5)     1-Vec2(5)   Vec2(5)     1-Vec2(5) ;
             1-Vec2(6)  Vec2(6)     1-Vec2(6)   Vec2(6)     1-Vec2(6)   Vec2(6)     1-Vec2(6)   Vec2(6) ;
             Vec2(7)    1-Vec2(7)   Vec2(7)     1-Vec2(7)   Vec2(7)     1-Vec2(7)   Vec2(7)     1-Vec2(7) ;
             1-Vec2(8)  Vec2(8)     1-Vec2(8)   Vec2(8)     1-Vec2(8)   Vec2(8)     1-Vec2(8)   Vec2(8) ];

% When a1 = -1 
MatAdjustMinus1 = [ zeros(1,8) ; zeros(1,8) ; ones(1,8) ; ones(1,8) ; ...
                    zeros(1,8) ; zeros(1,8) ; ones(1,8) ; ones(1,8) ] ; ...
MatAdjustMinus2 = [ ones(8,1)  ones(8,1)  zeros(8,1)  zeros(8,1)  ...
                    ones(8,1)  ones(8,1)  zeros(8,1)  zeros(8,1) ] ; 

output1 = TempMat0.*TempMat2.*MatAdjustMinus1.*MatAdjustMinus2; 

% When a1 = 0 
MatAdjustZero1 = [ 1 0 1 0 ; 0 1 0 1 ; 1 0 1 0 ; 0 1 0 1 ]; 
MatAdjustZero2 = ones(2,2); 
output2 = TempMat0.*TempMat2.*(kron(MatAdjustZero1, MatAdjustZero2)); 

% When a1 = 1 
MatAdjustPlus1 = [ ones(1,8) ; ones(1,8) ; zeros(1,8) ; zeros(1,8) ; ...
                   ones(1,8) ; ones(1,8) ; zeros(1,8) ; zeros(1,8) ]; ... 
MatAdjustPlus2 = [ zeros(8,1)  zeros(8,1) ones(8,1) ones(8,1) ...  
                   zeros(8,1)  zeros(8,1) ones(8,1) ones(8,1)  ] ; 

output3 = TempMat0.*TempMat2.*MatAdjustPlus1.*MatAdjustPlus2; 








            

