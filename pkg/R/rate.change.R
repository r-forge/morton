## new tree scaler -- equivalent to the non-censored O'Meara test
## the censored test has an extra degree of freedom, so the AIC 
## model approach swings away from detecting a change in rate
## Andrew Hipp, 21 feb 2011

require(ape)
require(phylobase)

rate.change <- function(x, ...) UseMethod('rate.change')

rate.change.phylo <- function(x, ...) {
  object = as(x, 'phylo4')
  rate.change.phylo4(object, ...)
  }

rate.change.phylo4 <- function(x, taxVector, dat, limits = c(0.001, 100), rateScale = 'optimize') {
  branches <- sharedBranches(x, taxVector)
  if(class(rateScale) == 'numeric') {
    x@edge.length[branches] <- x@edge.length[branches]*rateScale 
    x.analyzed <- fit.continuous(as(x, 'phylo'), dat)
    return(tr = x, dat.analyzed = dat, results = x.analyzed)
    } # end if
  else {
    scaler <- function(branchScalar, phy = x, branchesToScale = branches, dat.analyzed = dat) {
      phy@edge.length[branchesToScale] <- x@edge.length[branchesToScale]*branchScalar
      phy <- as(phy, 'phylo')
      sig <- sigmaSq(vcv(phy), dat.analyzed[phy$tip.label])
      logLik(sig)
      }
    a <- optimize(scaler, limits, maximum = T)
    x <- as(x, 'phylo')
    b <- logLik(sigmaSq(vcv(x), dat[x$tip.label]))
    a <- c(a[[1]], a[[2]], b)
    names(a) <- c('rate.change.ratio', 'lnL.rateChange', 'lnL.noRateChange')
    return(a)
    } #end else
  }
  
sharedBranches <- function(x, taxVector) {
  scaleNodes <- descendants(x, MRCA(x, taxVector), 'all')
  return(which(x@edge[, 2] %in% scaleNodes))
  }

rootState <- function(vcvMat = vcv(tr), X = x) {
  ## the maximum likelihood estimate of root state based on O'Meara, p. 925
  ## note that this is the same as doing gls(X ~ 1, correlation = corBrownian(phy=tr))$coef, but slower
  one <- matrix(1, length(X), 1) # a matrix of 1s, as tall as the lenght of X
  B0 <- as.numeric(solve(t(one) %*% solve(vcvMat) %*% one) %*% (t(one) %*% solve(vcvMat) %*% X))
  return(B0)
  }

sigmaSq <- function(vcvMat = vcv(tr), X = x) {
  ## following O'Meara 2006, eq. 2
  N <- length(X)
  E.X <- matrix(rootState(vcvMat, X), length(X), 1) # expected value of X is just the rootstate of X
  sigmaSq <- as.numeric((t(X - E.X) %*% solve(vcvMat) %*% (X - E.X)) / N)
  out <- list(V = vcvMat * sigmaSq, sigmaSq = sigmaSq, X = X, E.X. = E.X, N = N)
  class(out) <- "phylogSigmaSq"
  return(out)
  }

logLik.phylogSigmaSq <- function(object, ...) {
  ## following O'Meara 2006, eq. 3
  numer <- exp(-0.5 * t(object$X - object$E.X) %*% solve(object$V) %*% (object$X - object$E.X))
  denom <- sqrt((2*pi)^object$N * det(object$V))
  out <- log(numer / denom)
  attr(out, "nobs") <- object$N
  attr(out, "df") <- 2
  class(out) <- "logLik"
  return(out)
  }

