check_convergence <- function(oldCCP, newCCP, tol = 1e-12) {
  # CCPが収束しているかどうか確認する関数
  # input
  ## oldCCP: 更新前のCCP
  ## newCCP: 更新後のCCP
  ## tol: 収束判定の閾値
  # output
  ## res1 && res2: どちらも収束していればTRUE、それ以外ではFALSE

  # CCP1の収束判定
  res1 <- t(oldCCP[, 1] - newCCP[, 1]) %*% (oldCCP[, 1] - newCCP[, 1]) < tol

  # CCP2の収束判定
  res2 <- t(oldCCP[, 2] - newCCP[, 2]) %*% (oldCCP[, 2] - newCCP[, 2]) < tol

  # CCP1の差分とCCP2の差分を表示
  print(c(
    t(oldCCP[, 1] - newCCP[, 1]) %*% (oldCCP[, 1] - newCCP[, 1]),
    t(oldCCP[, 2] - newCCP[, 2]) %*% (oldCCP[, 2] - newCCP[, 2])
  ))

  # どちらか一方でも収束していなければFALSEを返す
  return(res1 && res2)
}
