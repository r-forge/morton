gls.r.squared <- function(x, tr) {
  e = x$resid
  V = ifelse(class(tr) == "phylo", vcv.phylo(tr), tr)
  Y = x$resid + x$fitted
  one <- matrix(1, length(Y), 1)
  a <- as.numeric(solve(t(one) %*% solve(V) %*% one) %*% (t(one) %*% solve(V) %*% Y))
  r.squared= 1-(t(e) %*% solve(V) %*% e) / (t(Y-a) %*% solve(V) %*% (Y-a))
  return(r.squared[1,1])
  }
  
