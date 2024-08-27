BBLobjective_NLS <- function(theta, NumPerturbations
                             , W1star, W2star, W1_all, W2_all){

  Param1 <- c(theta[1:5] , 1 ) 
  Param2 <- c(theta[6:10], 1 ) 
  
  fa1 <- array(rep(W1star,NumPerturbations),dim = c(6, 8,NumPerturbations))
  fa2 <- array(rep(W2star,NumPerturbations),dim = c(6, 8,NumPerturbations))
  
  diff1 <- fa1 - W1_all
  diff2 <- fa2 - W2_all 
  
  # reshapeする。行が各要素になる
  temp1 <- matrix(diff1,nrow =6, ncol =  8*NumPerturbations)
  temp2 <- matrix(diff2,nrow =6, ncol =  8*NumPerturbations)
  
  val1 <- t(temp1)%*%Param1
  val2 <- t(temp2)%*%Param2
  
  val1 <- pmin(val1, 0)
  val2 <- pmin(val2, 0)
  
  # すべてStackして返す
  objvec <- rbind(val1,val2)
  
  # Rはスカラーのみ最小化できる
  obj <- sum(objvec^2)
  
  obj
}