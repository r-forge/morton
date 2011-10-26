## new tree scaler -- equivalent to the non-censored O'Meara test
## the censored test has an extra degree of freedom, so the AIC 
## model approach swings away from detecting a change in rate
## Andrew Hipp, 21 feb 2011
## 2 Aug 2011: modified to allow for > 2 rates

require(ape)
require(phylobase)

rate.change <- function(x, ...) UseMethod('rate.change')

rate.change.phylo <- function(x, ...) {
  object = as(x, 'phylo4')
  rate.change.phylo4(object, ...)
  }

rate.change.phylo4 <- function(x, taxVector, taxVector2 = NULL, dat, limits = c(0.001, 100), rateScale = 'optimize', plotPhy = TRUE, nB = NULL) {
  assign('counter', counter + 1, .GlobalEnv) # this is a kludge... just to count trees
  branches <- sharedBranches(x, taxVector)
  if(!identical(taxVector2, NULL)) {
    branches2 <- sharedBranches(x, taxVector2)
	branches <- branches[!branches %in% branches2] # this is the portion of branches2 not in branches(1); assumes 2 nested in 1
	}
  if(class(rateScale) == 'numeric') {
    x@edge.length[branches] <- x@edge.length[branches]*rateScale[1]
	if(!identical(taxVector2, NULL)) x@edge.length[branches2] <- x@edge.length[branches2]*rateScale[2]
    x.analyzed <- fit.continuous(as(x, 'phylo'), dat)
    return(tr = x, dat.analyzed = dat, results = x.analyzed)
    } # end if
  if(identical(taxVector2, NULL)) { # this is the univariate rate change case
    scaler <- function(branchScalar, phy = x, branchesToScale = branches, dat.analyzed = dat) {
      phy@edge.length[branchesToScale] <- x@edge.length[branchesToScale]*branchScalar
      phy <- as(phy, 'phylo')
      sig <- sigmaSq(vcv(phy), dat.analyzed[phy$tip.label])
      logLik(sig)
      }
    a <- optimize(scaler, limits, maximum = T)
    x <- as(x, 'phylo')
    # s.sq <- sigmaSq(vcv(x), dat[x$tip.label])
	# b <- logLik(s.sq)
    x.rescaled <- x
	x.rescaled$edge.length[branches] <- x.rescaled$edge.length[branches]*a[[1]]
	s.sq.rescaled <- sigmaSq(vcv(x.rescaled), dat[x.rescaled$tip.label])
    # a <- c(a[[1]],  s.sq.rescaled$sigmaSq, s.sq.rescaled$E.X[1], s.sq$sigmaSq, s.sq$E.X[1], a[[2]], b)
	a <- c(a[[1]],  s.sq.rescaled$sigmaSq, s.sq.rescaled$E.X[1], a[[2]])
	names(a) <- c('rate.change.ratio', 'sigma.sq', 'root', 'lnL')
    } # end if
  else { # this is the multivariate rate change case
    scaler <- function(pars, phy = x, branchesToScale.1 = branches, branchesToScale.2 = branches2, dat.analyzed = dat) {
      phy@edge.length[branchesToScale.1] <- x@edge.length[branchesToScale.1]*pars[1]
      phy@edge.length[branchesToScale.2] <- x@edge.length[branchesToScale.2]*pars[2]	  	
      phy <- as(phy, 'phylo')
      sig <- sigmaSq(vcv(phy), dat.analyzed[phy$tip.label])
      -logLik(sig)
	  }
	a <- optim(c(1,1), scaler, method = "L-BFGS-B", lower = limits[1], upper = limits[2])
    x <- as(x, 'phylo')
    x.rescaled <- x
	x.rescaled$edge.length[branches] <- x.rescaled$edge.length[branches]*a$par[1]
	x.rescaled$edge.length[branches2] <- x.rescaled$edge.length[branches2]*a$par[2]
	s.sq.rescaled <- sigmaSq(vcv(x.rescaled), dat[x.rescaled$tip.label])
	a <- c(a$par[1:2],  s.sq.rescaled$sigmaSq, s.sq.rescaled$E.X[1], -a$value, a$convergence)
	names(a) <- c('rate.change.ratio.1', 'rate.change.ratio.2', 'sigma.sq', 'root', 'lnL', 'convergence')
	}
  if(plotPhy) {
	pdf(file = paste('./rateChangeTrees/', nB, counter, ".pdf", sep = ""))
	plot(x.rescaled)
	dev.off()
	}
  return(a) 
  }
  
sharedBranches <- function(x, taxVector) {
  scaleNodes <- descendants(x, MRCA(x, taxVector), 'all')
  #return(which(x@edge[, 2] %in% scaleNodes))
  a = which(x@edge[, 2] %in% scaleNodes)
  return(c(a, min(a)-1))
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

