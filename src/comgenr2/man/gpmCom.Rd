\name{gpmCom}
\alias{gpmCom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Genotype-phenotype mapping function for input into the cgSim function. 
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
gpmCom(n = "number of species", het.values = c(5, 21), allelic.range = c(0, 3))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
%%     ~~Describe \code{n} here~~
}
  \item{het.values}{
%%     ~~Describe \code{het.values} here~~
}
  \item{allelic.range}{
%%     ~~Describe \code{allelic.range} here~~
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
function (n = "number of species", het.values = c(5, 21), allelic.range = c(0, 
    3)) 
{
    if (n == "number of species") {
        n <- 25
    }
    com <- matrix(NA, nrow = n, ncol = 2)
    com[, 1] <- runif(n, het.values[1], het.values[2])
    com[, 2] <- runif(n, allelic.range[1], allelic.range[2])
    com. <- com
    com.[, 1] <- (com[, 1] - 0.5 * com[, 2])/2
    com.[, 2] <- (com[, 1] + 0.5 * com[, 2])/2
    return(com.)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
