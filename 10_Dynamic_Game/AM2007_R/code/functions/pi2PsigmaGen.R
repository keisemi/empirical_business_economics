pi2PsigmaGen <- function(pi2, Mat1){
  # input: 
  # 8 x 3 state-dependent term-wise profit of firm 2 (pi2)
  # 8 x 3 ccp of firm1 in equilibrium
  # ChoiceMaxNumはinputからなくしたほうがよさそう
  # output: 
  # 8 x 3 state-dependent term-wise profit of firm 2 in equilibrium
  
  pi2_dec <- pi2[,1] # 1x8
  pi2_0 <- pi2[,2]
  pi2_inc <- pi2[,3]
  
  pi2_dec <- (pi2_dec * Mat1) %*% matrix(c(1,1,1))
  pi2_0 <- (pi2_0 * Mat1) %*% matrix(c(1,1,1))
  pi2_inc <- (pi2_inc * Mat1) %*% matrix(c(1,1,1))
  
  output <- cbind(pi2_dec, pi2_0)
  output <- cbind(output, pi2_inc)
  
  # output <- pi2_dec %>% cbind(., pi2_0) %>% cbind(.,pi2_inc)
  
  output
  
}