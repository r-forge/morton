\name{beastIn}
\Rdversion{1.1}
\alias{beastIn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
beastIn(filename = NULL, nodeList = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{
%%     ~~Describe \code{filename} here~~
}
  \item{nodeList}{
%%     ~~Describe \code{nodeList} here~~
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
function(filename = NULL, nodeList = NULL) {
## Formats Beast log file for subsampling
  filename <- ifelse(identical(filename, NULL), choose.files(multi = F, caption = "Select the BEAST log file with which you want to work"), filename)
  beastData <- read.table(filename, header = T)
  if(identical(nodeList, NULL)) nodeList = select.list(names(beastData), preselect = NULL, multiple = T, title = "Select the nodes that interest you")
  outData <- list(dataMatrix = beastData, nodes = nodeList)
  class(outData) <- 'beastData'
  return(outData)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
