BBLobjective <- function(theta, NumPerturbations
, W1star, W2star, W1_all, W2_all){
# x <- initial
# theta <- c(x[1], x[3:4], 0, x[5], x[2],x[3:4], 0, x[5])
  
  Param1 <- c(theta[1:5] , 1 ) 
  Param2 <- c(theta[6:10], 1 ) 
  
  objvalue <- 0
  
  for (per in 1:NumPerturbations){
  
    temp1 <- pmin(t(W1star)%*%Param1 - t(W1_all[ , ,per])%*%Param1,0)
    temp2 <- pmin(t(W2star)%*%Param2 - t(W2_all[ , ,per])%*%Param2,0)
  
    objvalue <- objvalue + t(temp1)%*%temp1 + t(temp2)%*%temp2
  }
  
  objvalue
}