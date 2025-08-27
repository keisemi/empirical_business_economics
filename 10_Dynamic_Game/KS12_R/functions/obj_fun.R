obj_fun <- function(param, EstimatedCCP1, EstimatedCCP2, TransitionProb, beta) {
  # パラメタの評価をするための目的関数

  # 1st stageで推定したCCPを元にValueを計算し、そのValue及びパラメタにもとづいてCCPを予測する
  output <- CCP_to_Value_to_prediction(param, EstimatedCCP1, EstimatedCCP2, TransitionProb, beta)

  # 目的関数の値を計算する
  obj <- sum((output[, 1] - EstimatedCCP1)^2) + sum((output[, 2] - EstimatedCCP2)^2)

  obj
}
