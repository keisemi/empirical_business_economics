function [output1, output2, output3] = fP_a2given(Matrix1, Vec1)


TempMat0 = kron(Matrix1,ones(4));

% State transition due to rival's action 
TempMat1 = [ Vec1(1)   Vec1(1)    1-Vec1(1) 1-Vec1(1) Vec1(1) Vec1(1) 1-Vec1(1) 1-Vec1(1) ;
             Vec1(2)   Vec1(2)    1-Vec1(2) 1-Vec1(2) Vec1(2) Vec1(2) 1-Vec1(2) 1-Vec1(2) ;
             1-Vec1(3) 1-Vec1(3)  Vec1(3)   Vec1(3)   1-Vec1(3) 1-Vec1(3)  Vec1(3)   Vec1(3) 
             1-Vec1(4) 1-Vec1(4)  Vec1(4)   Vec1(4)   1-Vec1(4) 1-Vec1(4)  Vec1(4)   Vec1(4) 
             Vec1(5)   Vec1(5)    1-Vec1(5) 1-Vec1(5) Vec1(5) Vec1(5) 1-Vec1(5) 1-Vec1(5) ;
             Vec1(6)   Vec1(6)    1-Vec1(6) 1-Vec1(6) Vec1(6) Vec1(6) 1-Vec1(6) 1-Vec1(6) ;
             1-Vec1(7) 1-Vec1(7)  Vec1(7)   Vec1(7)   1-Vec1(7) 1-Vec1(7)  Vec1(7)   Vec1(7);
             1-Vec1(8) 1-Vec1(8)  Vec1(8)   Vec1(8)   1-Vec1(8) 1-Vec1(8)  Vec1(8)   Vec1(8)];


% When a_2 = -1 
MatAdjustMinus1 = [ zeros(1,8) ; ones(1,8) ; zeros(1,8) ; ones(1,8) ; ...
                    zeros(1,8) ; ones(1,8) ; zeros(1,8) ; ones(1,8) ] ; 
MatAdjustMinus2 = [ ones(8,1) zeros(8,1) ones(8,1) zeros(8,1) ...
                    ones(8,1) zeros(8,1) ones(8,1) zeros(8,1) ]; 
output1 = TempMat0.*TempMat1.*MatAdjustMinus1.*MatAdjustMinus2;  


% When a_2 = 0
MatAdjustZero1 = ones(4,4);
MatAdjustZero2 = eye(2,2); 
output2 = TempMat0.*TempMat1.*(kron(MatAdjustZero1, MatAdjustZero2));  


% When a_2 = 1 
MatAdjustPlus1 = [ ones(1,8) ; zeros(1,8) ; ones(1,8) ; zeros(1,8) ; ...
                   ones(1,8) ; zeros(1,8) ; ones(1,8) ; zeros(1,8) ]; 
MatAdjustPlus2 = [ zeros(8,1) ones(8,1) zeros(8,1) ones(8,1) ... 
                   zeros(8,1) ones(8,1) zeros(8,1) ones(8,1) ]; 

output3 = TempMat0.*TempMat1.*MatAdjustPlus1.*MatAdjustPlus2;  


