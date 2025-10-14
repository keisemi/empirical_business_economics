function output = fP(Matrix1, Vec1, Vec2) 

% Exogenous state transition probabilities 
TempMat0 = kron(Matrix1,ones(4));

% State transition due to their actions 
TempMat1 = [ Vec1(1)   Vec1(1)    1-Vec1(1) 1-Vec1(1) Vec1(1) Vec1(1) 1-Vec1(1) 1-Vec1(1) ;
             Vec1(2)   Vec1(2)    1-Vec1(2) 1-Vec1(2) Vec1(2) Vec1(2) 1-Vec1(2) 1-Vec1(2) ;
             1-Vec1(3) 1-Vec1(3)  Vec1(3)   Vec1(3)   1-Vec1(3) 1-Vec1(3)  Vec1(3)   Vec1(3) ;
             1-Vec1(4) 1-Vec1(4)  Vec1(4)   Vec1(4)   1-Vec1(4) 1-Vec1(4)  Vec1(4)   Vec1(4) ;
             Vec1(5)   Vec1(5)    1-Vec1(5) 1-Vec1(5) Vec1(5) Vec1(5) 1-Vec1(5) 1-Vec1(5) ;
             Vec1(6)   Vec1(6)    1-Vec1(6) 1-Vec1(6) Vec1(6) Vec1(6) 1-Vec1(6) 1-Vec1(6) ;
             1-Vec1(7) 1-Vec1(7)  Vec1(7)   Vec1(7)   1-Vec1(7) 1-Vec1(7)  Vec1(7)   Vec1(7);
             1-Vec1(8) 1-Vec1(8)  Vec1(8)   Vec1(8)   1-Vec1(8) 1-Vec1(8)  Vec1(8)   Vec1(8)];

TempMat2 = [ Vec2(1)    1-Vec2(1)   Vec2(1)     1-Vec2(1)   Vec2(1)     1-Vec2(1)   Vec2(1)     1-Vec2(1) ; 
             1-Vec2(2)  Vec2(2)     1-Vec2(2)   Vec2(2)     1-Vec2(2)   Vec2(2)     1-Vec2(2)   Vec2(2) ;
             Vec2(3)    1-Vec2(3)   Vec2(3)     1-Vec2(3)   Vec2(3)     1-Vec2(3)   Vec2(3)     1-Vec2(3) ;
             1-Vec2(4)  Vec2(4)     1-Vec2(4)   Vec2(4)     1-Vec2(4)   Vec2(4)     1-Vec2(4)   Vec2(4) ;
             Vec2(5)    1-Vec2(5)   Vec2(5)     1-Vec2(5)   Vec2(5)     1-Vec2(5)   Vec2(5)     1-Vec2(5) ;
             1-Vec2(6)  Vec2(6)     1-Vec2(6)   Vec2(6)     1-Vec2(6)   Vec2(6)     1-Vec2(6)   Vec2(6) ;
             Vec2(7)    1-Vec2(7)   Vec2(7)     1-Vec2(7)   Vec2(7)     1-Vec2(7)   Vec2(7)     1-Vec2(7) ;
             1-Vec2(8)  Vec2(8)     1-Vec2(8)   Vec2(8)     1-Vec2(8)   Vec2(8)     1-Vec2(8)   Vec2(8) ];

% The output is a Hadamard product of TempMat0, TempMat1, and TempMat2
output = TempMat0.*TempMat1.*TempMat2; 

end