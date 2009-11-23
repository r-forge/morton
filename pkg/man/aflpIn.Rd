\name{aflpIn}
\Rdversion{1.1}
\alias{aflpIn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
aflpIn(filename = NULL, interactive = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{
%%     ~~Describe \code{filename} here~~
}
  \item{interactive}{
%%     ~~Describe \code{interactive} here~~
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
function(filename = NULL, interactive = T) {
if (interactive) dat <- read.table(file.choose(), as.is = T, sep = "", fill = T, col.names = c('inds', 'loci'))
  inds <- as.character(dat$inds)
  loci <- as.character(dat$loci)
  outdata <- matrix(data = 0, 
                    nrow = length(unique(inds[!is.na(inds)])), 
                    ncol = length(unique(loci[!is.na(loci)])), 
                    dimnames = list(sort(unique(inds[!is.na(inds)])), sort(unique(loci[!is.na(loci)])))
                    )
  for(i in 1:length(inds)) {
    if (is.na(loci[i])) next
    else outdata[inds[i], loci[i]] <- 1
    }
  return(outdata)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
