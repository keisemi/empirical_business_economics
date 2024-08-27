CCP1LogTransform <- function(x){
  # input: 8 x 1 CCP essence vector
  # output: 8 x 3 logged CCP matrix
  
  output <- matrix(c(0, log(x[1]), log(1-x[1]),
                     0, log(x[2]), log(1-x[2]),
                     log(1-x[3]), log(x[3]), 0,
                     log(1-x[4]), log(x[4]), 0,
                     0, log(x[5]), log(1-x[5]),
                     0, log(x[6]), log(1-x[6]),
                     log(1-x[7]), log(x[7]), 0,
                     log(1-x[8]), log(x[8]), 0),
                   ncol = 3, byrow = TRUE)
  output # return output
  
}