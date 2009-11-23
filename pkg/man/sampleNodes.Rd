\name{sampleNodes}
\Rdversion{1.1}
\alias{sampleNodes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sampleNodes(beastData, ratesVector = kay2006$rate[kay2006$lifeHistory == "H"][1:10], subsample = 5000, CI = 0.95, burnin = round(dim(beastData$dataMatrix)[1] * 0.4))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{beastData}{
%%     ~~Describe \code{beastData} here~~
}
  \item{ratesVector}{
%%     ~~Describe \code{ratesVector} here~~
}
  \item{subsample}{
%%     ~~Describe \code{subsample} here~~
}
  \item{CI}{
%%     ~~Describe \code{CI} here~~
}
  \item{burnin}{
%%     ~~Describe \code{burnin} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
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
function(beastData, ratesVector = kay2006$rate[kay2006$lifeHistory == 'H'][1:10], subsample = 5000, CI = 0.95, burnin = round(dim(beastData$dataMatrix)[1] * 0.4)) {
## Samples depths of selected nodes at random from the Markov chain, returning for each (1) a vector of subsampled parameters, and (2) the same vector * a vector of subsampled rates
## NOTE: This analysis assumes that nodeDepth / rate is a meaningful value (e.g., that rate is in substitutions / site / year, and node depth is in substitution / site)
  if(class(beastData) != 'beastData') message('The data were not formatted using beastIn... I will try to analyze, but who knows what might happen?!')
  postBurnin = dim(beastData$dataMatrix)[1] - burnin
  if(postBurnin < subsample) stop(paste('There are only', postBurnin, 'elements in the post-burnin portion of your run...\nReduce your subsample size and rerun'))
  nodesVector <- beastData$nodes
  nodeDepthSubst <- matrix(nrow = subsample, ncol = length(nodesVector), dimnames = list(seq(subsample), nodesVector))
  nodeDepthYears <- nodeDepthSubst
  CIlistSubst <- vector('list', length(nodesVector)); names(CIlistSubst) <- nodesVector
  CIlistYears <- CIlistSubst
  meanSubst <- numeric(length(nodesVector)); names(meanSubst) <- nodesVector
  meanYear <- meanSubst
  CIelements = c(ceiling(subsample * (1 - CI) / 2), floor(subsample * (1 + CI) / 2))
  for(i in nodesVector) {
    nodeDepthSubst[ ,i] <- sample(beastData$dataMatrix[[i]][burnin : dim(beastData$dataMatrix)[1]], subsample, replace = F)
    nodeDepthYears[ ,i] <- nodeDepthSubst[ ,i] / sample(ratesVector, subsample, replace = T) 
    substSort <- sort(nodeDepthSubst[ ,i]); yearsSort = sort(nodeDepthYears[ ,i])
    CIlistSubst[[i]] <- c(substSort[CIelements[1]], substSort[CIelements[2]])
    CIlistYears[[i]] <- c(yearsSort[CIelements[1]], yearsSort[CIelements[2]])
    meanSubst[i] <- mean(nodeDepthSubst[ ,i])
    meanYear[i] <- mean(nodeDepthYears[ ,i]) }
  outData = list(substMatrix = nodeDepthSubst, yearsMatrix = nodeDepthYears, CI = list(subst = CIlistSubst, years = CIlistYears), mean = list(subst = meanSubst, year = meanYear), nodesVector = nodesVector, CIwidth = CI)
  class(outData) = 'sampleNodes'
  return(outData)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
