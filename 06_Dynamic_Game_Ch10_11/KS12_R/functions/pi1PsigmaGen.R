pi1PsigmaGen <- function(pi1, Mat2) {
  # input:
  # 8 x 3 state-dependent term-wise profit of firm 1 (pi1)
  # 8 x 3 ccp of firm2 in equilibrium
  # output:
  # 8 x 3 state-dependent term-wise profit of firm 1 in equilibrium

  pi1_dec <- pi1[, 1] # 1x8
  pi1_0 <- pi1[, 2]
  pi1_inc <- pi1[, 3]

  pi1_dec <- (pi1_dec * Mat2) %*% matrix(c(1, 1, 1))
  pi1_0 <- (pi1_0 * Mat2) %*% matrix(c(1, 1, 1))
  pi1_inc <- (pi1_inc * Mat2) %*% matrix(c(1, 1, 1))

  output <- cbind(pi1_dec, pi1_0)
  output <- cbind(output, pi1_inc)

  # output <-
  #   pi1_dec %>% cbind(., pi1_0) %>% cbind(., pi1_inc)

  output
}
