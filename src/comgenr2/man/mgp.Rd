\name{mgp}
\alias{mgp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mgp(scn, com, loc = TRUE, my.coord = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{scn}{
%%     ~~Describe \code{scn} here~~
}
  \item{com}{
%%     ~~Describe \code{com} here~~
}
  \item{loc}{
%%     ~~Describe \code{loc} here~~
}
  \item{my.coord}{
%%     ~~Describe \code{my.coord} here~~
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
function (scn, com, loc = TRUE, my.coord = "") 
{
    e.col <- sign(scn)
    e.col[e.col == 1] <- "grey"
    e.col[e.col == -1] <- "red"
    v.cex <- apply(com, 2, sum)
    v.cex <- log(v.cex, 10)
    v.cex <- v.cex * (1/min(v.cex))
    v.cex <- v.cex/2
    if (length(my.coord) == 1) {
        coord <- gplot(abs(scn), displaylabels = TRUE, gmode = "graph", 
            pad = 1.5, edge.col = e.col, edge.lwd = abs(scn), 
            vertex.cex = v.cex, vertex.col = "darkgrey", vertex.border = "darkgrey")
    }
    else {
        coord <- gplot(abs(scn), displaylabels = TRUE, gmode = "graph", 
            pad = 1.5, edge.col = e.col, edge.lwd = abs(scn), 
            vertex.cex = v.cex, vertex.col = "darkgrey", vertex.border = "darkgrey", 
            coord = my.coord)
    }
    if (loc) {
        return(coord)
    }
    else {
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
