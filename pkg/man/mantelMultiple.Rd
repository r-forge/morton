\name{mantelMultiple}
\Rdversion{1.1}
\alias{mantelMultiple}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mantelMultiple(Y, X, permutations = 999, strata)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
   This should be a dist object or a matrix... I'm not sure, so try both ways. This is your response variable.
}
  \item{X}{
   This is a list of dist objects or matrices... again, try both ways. I'll fix this in the real documentation. Each item in this list is a predictor of Y.
}
  \item{permutations}{
   The number of permutations of Y that you want to assess the significance of partial correlations.
}
  \item{strata}{
   Not implemented; a holdover from vegan.
}
}
\details{
   It's imperative that the matrices you pass in are ordered in the same way (i.e., samples are in the same order along the rows and the columns). The safest way to do this is to set a vector of individuals that you want to sample and subset every matrix by x[inds, inds], where x is a matrix and inds is a vector of individuals you want to sample.
   
   Watch out! NAs are not supported and may behave in unpredictable ways.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(Y, X, permutations = 999, strata) {
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
