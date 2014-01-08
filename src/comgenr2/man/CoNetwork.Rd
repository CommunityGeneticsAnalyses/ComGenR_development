\name{CoNetwork}
\alias{CoNetwork}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
CoNetwork(x, plot.net = TRUE, scalar = 3, min.vsize = 0.1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{plot.net}{
%%     ~~Describe \code{plot.net} here~~
}
  \item{scalar}{
%%     ~~Describe \code{scalar} here~~
}
  \item{min.vsize}{
%%     ~~Describe \code{min.vsize} here~~
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
function (x, plot.net = TRUE, scalar = 3, min.vsize = 0.1) 
{
    bc.d <- as.matrix(vegdist(t(x)))
    prune <- co.net(x)
    bc.d[prune == 0] <- 0
    thresh <- percThreshold(bc.d)$threshold
    pruned.net <- bc.d
    pruned.net[bc.d < thresh] <- 0
    if (plot.net) {
        v.cex <- apply(x, 2, sum)
        v.cex <- (((v.cex/sum(v.cex))/max((v.cex/sum(v.cex)))) * 
            scalar) + min.vsize
        gplot(abs(pruned.net), displaylabels = TRUE, gmode = "graph", 
            pad = 1.5, edge.lwd = (abs(pruned.net)), vertex.cex = v.cex, 
            vertex.col = "grey")
    }
    return(pruned.net)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
