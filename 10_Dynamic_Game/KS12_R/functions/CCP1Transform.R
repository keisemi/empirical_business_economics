CCP1Transform <- function(x) {
  # input: 8 x 1 CCP essence vector
  # output: 8 x 3 CCP matrix

  output <- matrix(
    c(
      0, x[1], 1 - x[1],
      0, x[2], 1 - x[2],
      1 - x[3], x[3], 0,
      1 - x[4], x[4], 0,
      0, x[5], 1 - x[5],
      0, x[6], 1 - x[6],
      1 - x[7], x[7], 0,
      1 - x[8], x[8], 0
    ),
    ncol = 3, byrow = TRUE
  )
  output # return output
}
