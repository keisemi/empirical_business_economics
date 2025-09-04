pi2gen <- function(theta) {
  # input: 10 x 1 parameter matrix
  # output: 8 x 3 term profit of firm 2, dependent on
  # state of the world (row) and action (col);
  # some actions are unfeasible and will be adjusted outside this function

  base <- matrix(c(
    0,
    theta[6] + theta[8],
    0,
    theta[6] + theta[7] + theta[8],
    0,
    theta[6],
    0,
    theta[6] + theta[7]
  )) # 8 x 1

  invdiv <- matrix(c(theta[9], 0, theta[10]), nrow = 1) # 1 x 3

  base <- rep(base, 3)
  base <- matrix(base, nrow = 8)
  invdiv <- rep(invdiv, 8)
  invdiv <- matrix(invdiv, nrow = 8, byrow = TRUE)
  output <- base + invdiv

  # output <-
  #   base |> rep(3) |> matrix(nrow = 8) +
  #   invdiv |> rep(8) |> matrix(nrow = 8, byrow = TRUE)

  output # return output
}
