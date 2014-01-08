\name{percThreshold}
\alias{percThreshold}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
percThreshold(x = "network matrix", step.size = 0.01)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{step.size}{
%%     ~~Describe \code{step.size} here~~
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
function (x = "network matrix", step.size = 0.01) 
{
    no.c <- no.clusters(graph.adjacency(x, weighted = TRUE))
    step <- 1
    while (no.c == 1) {
        x[x <= (step * step.size)] <- 0
        no.c <- no.clusters(graph.adjacency(x, weighted = TRUE))
        step <- step + 1
    }
    out <- list(threshold = ((step - 1) * step.size), isolated.nodes = colnames(x)[clusters(graph.adjacency(x, 
        weighted = TRUE))$membership == 2])
    return(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
