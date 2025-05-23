fP <- function(Matrix1, Vec1, Vec2) {
  # input:
  # 2 x 2 exogenous state transition matrix (Matrix1);
  # 8 x 1 CCP essence vector of firm 1(Vec1)
  # 8 x 1 CCP essence vector of firm 2(Vec2)
  # output:
  # 8 x 8 state transition matrix given CCPs

  # exogenous state transition probabilities; expansion to 8x8
  TempMat0 <-
    kronecker(Matrix1, matrix(rep(1, 16), ncol = 4))

  # state transition due to their actions
  TempMat1 <- matrix(c(
    Vec1[1], 1 - Vec1[1],
    Vec1[2], 1 - Vec1[2],
    1 - Vec1[3], Vec1[3],
    1 - Vec1[4], Vec1[4],
    Vec1[5], 1 - Vec1[5],
    Vec1[6], 1 - Vec1[6],
    1 - Vec1[7], Vec1[7],
    1 - Vec1[8], Vec1[8]
  ), ncol = 2, byrow = TRUE)
  TempMat1 <- kronecker(TempMat1, matrix(c(1, 1), nrow = 1))
  TempMat1 <- cbind(TempMat1, TempMat1)

  # the above lines effectively conduct the below lines [CONFIRM]:
  # TempMat1 <- matrix(c(Vec1[1], Vec1[1], 1-Vec1[1], 1-Vec1[1], Vec1[1], Vec1[1], 1-Vec1[1], 1-Vec1[1],
  #                      Vec1[2], Vec1[2], 1-Vec1[2], 1-Vec1[2], Vec1[2], Vec1[2], 1-Vec1[2], 1-Vec1[2],
  #                      1-Vec1[3], 1-Vec1[3], Vec1[3], Vec1[3], 1-Vec1[3], 1-Vec1[3], Vec1[3], Vec1[3],
  #                      1-Vec1[4], 1-Vec1[4], Vec1[4], Vec1[4], 1-Vec1[4], 1-Vec1[4], Vec1[4], Vec1[4],
  #                      Vec1[5], Vec1[5], 1-Vec1[5], 1-Vec1[5], Vec1[5], Vec1[5], 1-Vec1[5], 1-Vec1[5],
  #                      Vec1[6], Vec1[6], 1-Vec1[6], 1-Vec1[6], Vec1[6], Vec1[6], 1-Vec1[6], 1-Vec1[6],
  #                      1-Vec1[7], 1-Vec1[7], Vec1[7], Vec1[7], 1-Vec1[7], 1-Vec1[7], Vec1[7], Vec1[7],
  #                      1-Vec1[8], 1-Vec1[8], Vec1[8], Vec1[8], 1-Vec1[8], 1-Vec1[8], Vec1[8], Vec1[8]),
  #                    ncol = 8, byrow = TRUE)

  TempMat2 <- matrix(c(
    Vec2[1], 1 - Vec2[1],
    1 - Vec2[2], Vec2[2],
    Vec2[3], 1 - Vec2[3],
    1 - Vec2[4], Vec2[4],
    Vec2[5], 1 - Vec2[5],
    1 - Vec2[6], Vec2[6],
    Vec2[7], 1 - Vec2[7],
    1 - Vec2[8], Vec2[8]
  ), ncol = 2, byrow = TRUE)
  TempMat2 <- kronecker(matrix(rep(1, 4), nrow = 1), TempMat2)
  # [CONFIRM]

  # The output is the Hadamard product of TempMat0, TempMat1, and TempMat2
  output <- TempMat0 * TempMat1 * TempMat2

  output
}
