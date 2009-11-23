\name{haversine}
\Rdversion{1.1}
\alias{haversine}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
haversine(lat1, long1, lat2, long2, r = 6372.795)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lat1}{
%%     ~~Describe \code{lat1} here~~
}
  \item{long1}{
%%     ~~Describe \code{long1} here~~
}
  \item{lat2}{
%%     ~~Describe \code{lat2} here~~
}
  \item{long2}{
%%     ~~Describe \code{long2} here~~
}
  \item{r}{
%%     ~~Describe \code{r} here~~
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
function(lat1, long1, lat2, long2, r = 6372.795) {
  lat1 = lat1 * pi / 180
  long1 = long1 * pi / 180
  lat2 = lat2 * pi / 180
  long2 = long2 * pi / 180
  deltaLong = long2 - long1
  deltaLat = lat2 - lat1
  a = sin(deltaLat/2)^2 + cos(lat1)*cos(lat2)*sin(deltaLong/2)^2
  if(length(dim(a) == 2)) a <- apply(a, c(1,2), min, 1)
   else a <- min(a,1)
  c = 2 * asin(sqrt(a))
  d = r * c
  return(d) }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
