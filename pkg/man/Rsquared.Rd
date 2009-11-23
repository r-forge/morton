\name{Rsquared}
\Rdversion{1.1}
\alias{Rsquared}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Rsquared(Y, X)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
%%     ~~Describe \code{Y} here~~
}
  \item{X}{
%%     ~~Describe \code{X} here~~
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
function (Y, X) {
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
    R2 <- Ryx \%*\% solve(Rxx, Rxy)
    return(as.vector(R2))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
