fP_a2given <- function(Matrix1, Vec1){
  # input: 
  # 2 x 2 exogenous state transition matrix (Matrix1); 
  # 8 x 1 CCP essence vector of firm 1(Vec1)
  # output:
  # list of 3 8x8 state transition matrices
  # some lines are same as in fP funciton; 
  # adjustnment matrix multiplication rules out some impossible transitions
  
  # exogenous state transition probabilities; expansion to 8x8
  TempMat0 <-
    kronecker(Matrix1, matrix(rep(1, 16), ncol = 4))
  
  TempMat1 <- matrix(c(Vec1[1], Vec1[1],1-Vec1[1], 1-Vec1[1],
                       Vec1[2], Vec1[2],1-Vec1[2], 1-Vec1[2],
                       1-Vec1[3], 1-Vec1[3], Vec1[3],Vec1[3], 
                       1-Vec1[4], 1-Vec1[4], Vec1[4],Vec1[4], 
                       Vec1[5], Vec1[5],1-Vec1[5], 1-Vec1[5],
                       Vec1[6], Vec1[6],1-Vec1[6], 1-Vec1[6],
                       1-Vec1[7], 1-Vec1[7], Vec1[7],Vec1[7], 
                       1-Vec1[8], 1-Vec1[8], Vec1[8],Vec1[8] ), ncol = 4, byrow = TRUE)
  TempMat1 <- kronecker(matrix(rep(1, 2), nrow = 1), TempMat1)
  
  # when a_2 = -1
  # 
  MatAdjustMinus1 <- matrix(rep(rep(c(0, 1), each = 8),4), nrow =8 ,byrow = TRUE) 
  MatAdjustMinus2 <- matrix(rep(rep(c(1, 0), each = 8),4), nrow =8 ) 
  output1 <- TempMat0 * TempMat1 * MatAdjustMinus1 * MatAdjustMinus2
  
  # when a_2 = 0
  ForZero = c(1, 0, 1, 0 , 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1)
  MatAdjustZero <- matrix(rep(ForZero, 4) , ncol = 8 ,byrow = TRUE)
  output2 <- TempMat0 * TempMat1 * MatAdjustZero
  
  # when a_2 = 1
  MatAdjustPlus1 <- matrix(rep(rep(c(1, 0), each = 8),4), nrow =8 , byrow = TRUE ) 
  MatAdjustPlus2 <- matrix(rep(rep(c(0, 1), each = 8),4), nrow =8 ) 
  output3 <- TempMat0 * TempMat1 * MatAdjustPlus1 * MatAdjustPlus2
  
  output <- list(output1, output2, output3)
  
  output
}
