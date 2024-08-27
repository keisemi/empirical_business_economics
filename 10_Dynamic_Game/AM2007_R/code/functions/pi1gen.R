pi1gen <- function(theta){
  # input: 10 x 1 parameter vector
  # output: 8 x 3 term profit of firm 1, dependent on 
  # state of the world (row) and action (col); 
  # some actions are unfeasible and will be adjusted outside this function

  base <- matrix(c(0,
                  0,
                  theta[1] + theta[3],
                  theta[1] + theta[2] + theta[3],
                  0,
                  0,
                  theta[1],
                  theta[1] + theta[2])) # 8 x 1
  
  invdiv <- matrix(c(theta[4], 0, theta[5]), nrow = 1) # 1 x 3
  
  
  base <- rep(base, 3)
  base <- matrix(base, nrow = 8)
  invdiv <- rep(invdiv, 8)
  invdiv <- matrix(invdiv, nrow = 8, byrow = TRUE)
  output <- base + invdiv
  
  # # pipe
  # output <-
  #   base %>% rep(., 3) %>% matrix(., nrow = 8) +
  #   invdiv %>% rep(., 8) %>% matrix(., nrow = 8, byrow = TRUE)
  
  output # return output
}
