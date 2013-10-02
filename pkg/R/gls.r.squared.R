gls.r.squared <- function(x) {
# based on Judge et al. 1985, eq. 2.3.16
  e = x$resid
  V <- corMatrix(x$modelStruct$corStruct)
  # V = vcv.phylo(tr)
  Y = x$resid + x$fitted
  one <- matrix(1, length(Y), 1)
  a <- as.numeric(solve(t(one) %*% solve(V) %*% one) %*% (t(one) %*% solve(V) %*% Y))
  r.squared= 1-(t(e) %*% solve(V) %*% e) / (t(Y-a) %*% solve(V) %*% (Y-a))
  return(r.squared[1,1])
  }
