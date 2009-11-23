mantelMultiple <- function(Y, X, permutations = 999, strata) {
## Y must be a single distance matrix, the response
## X is a list of distance matrices, the predictors (can be a list of length == 1 for a standard Mantel test)
## The Pearson product moment is the hard-wired correlation coefficient as written (though of course any statistic could be used)
## Andrew Hipp, 22 May 2008 (ahipp@mortonarb.org)

## 0. Format data
  if(class(Y) == 'list' || class(X) != 'list') stop('Y should be a solitary vector or matrix, the reponse; X should be a list of vectors or matrices in the same form; if X is solitary, still wrap it as list(X).')
  if (identical(names(X), NULL)) names(X) <- paste('X', 1:length(X), sep = '')
  Ylist <- list(Y); names(Ylist) <- 'Y'
  matrixList <- c(X, Ylist)
  pcors <- numeric(length(X)); names(pcors) <- names(X)
  pcorsDistTemp <- pcors
  RHsignif <- pcors
  LHsignif <- pcors
  pcorDistribution <- matrix(nrow = permutations, ncol = length(X), dimnames = list(seq(permutations), names(X))) 

## 1. Get the partial and multiple correlation coefficients for all items in X
  pCorMatrix = partialCorrelation(matrixList)
  for(i in names(X)) pcors[i] <- pCorMatrix[i, 'Y']
  multRsquared = Rsquared(Y, X)

## 2. Generate a distribution of partial correlation coefficients by permutation
## The following lines are modified (slightly) from vegan 1.11-2 function 'mantel.partial'
## ---------------------------------------------------------------------------------------
	N <- attributes(Y)$Size
	perm <- rep(0, permutations)
	for (i in 1:permutations) {
	  take <- permutedIndex(N, strata)
	  permvec <- as.dist(as.matrix(Y)[take, take])
	  permVeclist <- list(permvec); names(permVeclist) <- 'Y'
	  permMatrixList <- c(X, permVeclist)
	  pcorsDistTemp = partialCorrelation(permMatrixList)
	  for(j in names(X)) pcorDistribution[i,j] <- pcorsDistTemp[j, 'Y']
	#  rxy <- cor(permvec, ydis, method = method)
	#  rxz <- cor(permvec, zdis, method = method)
	#  perm[i] <- part.cor(rxy, rxz, ryz)
	  }
## ---------------------------------------------------------------------------------------

  
## 3. Calculate RH and LH significance
for(i in names(X)) {
  RHsignif[i] <- sum(pcorDistribution[,i] >= pcors[i]) / permutations
  LHsignif[i] <- sum(pcorDistribution[,i] <= pcors[i]) / permutations
  }  

output = list(pcors = pcors, partialRsquared = pcors^2, multRsquared = multRsquared, pValues = rbind(RHsignif, LHsignif))
return(output)
}


Rsquared <- function (Y, X) {
## multiple R squared for an indefinite number of X (predictor) vectors or matrices on a single Y (response) vector or matrix
    if(class(Y) == 'list' || class(X) != 'list') stop('Y should be a solitary vector or matrix, the reponse; X should be a list of vectors or matrices in the same form; if X is solitary, still wrap it as list(X).')
    Ylist <- list(Y); names(Ylist) <- 'Y'
    matrixList <- c(X, Ylist)
    matrices <- names(matrixList)
    R = matrix(nrow = length(matrices), ncol = length(matrices), dimnames = list(matrices,matrices))
    for (i in matrices) {
      for (j in matrices) {
        R[i,j] <- cor(matrixList[[i]], matrixList[[j]]) }}
    Rxx <- R[1:(length(matrices) - 1), 1:(length(matrices) - 1)]
    Rxy <- R[1:(length(matrices) - 1), length(matrices)]
    Ryx <- R[length(matrices), 1:(length(matrices) - 1)]
    R2 <- Ryx %*% solve(Rxx, Rxy)
    return(as.vector(R2))
}

partialCorrelation <- function(matrixList) {
require(corpcor)
    matrices <- names(matrixList)
    R <- matrix(nrow = length(matrices), ncol = length(matrices), dimnames = list(matrices,matrices))
    for (i in matrices) {
      for (j in matrices) {
        R[i,j] <- cor(matrixList[[i]], matrixList[[j]]) }}
    output <- cor2pcor(R, tol = 0.0001)
    dimnames(output) <- list(matrices,matrices)
    return(output)
    }

permutedIndex <- function (n, strata) 
## This is straight out of vegan 1.11-2 (originally permuted.index, renamed here in case the original function should change)
{
    if (missing(strata) || is.null(strata)) 
        out <- sample(n, n)
    else {
        out <- 1:n
        inds <- names(table(strata))
        for (is in inds) {
            gr <- out[strata == is]
            if (length(gr) > 1) 
                out[gr] <- sample(gr, length(gr))
        }
    }
    out
}
