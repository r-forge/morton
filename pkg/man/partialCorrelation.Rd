\name{partialCorrelation}
\Rdversion{1.1}
\alias{partialCorrelation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
partialCorrelation(matrixList)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{matrixList}{
%%     ~~Describe \code{matrixList} here~~
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
function(matrixList) {
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
