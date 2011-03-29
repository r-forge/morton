tTest.slope <- function(obj1, obj2) {
# source: https://stat.ethz.ch/pipermail/r-help/2009-April/195075.html
# obj1, obj2 are gls linear regression objects, each with a single predictor
  s1 = summary(obj1)$tTable
  s2 = summary(obj2)$tTable
  db <- (s2[2,1]-s1[2,1])
  pooled.sd <- sqrt(s2[2,2]^2+s1[2,2]^2)
  df.residual1 <- summary(obj1)$dims$N - summary(obj1)$dims$p
  df.residual2 <- summary(obj2)$dims$N - summary(obj2)$dims$p
  df.full = df.residual1 + df.residual2
  td <- db/pooled.sd
  p <- 2*pt(-abs(td), df.full)
  return(p)
}