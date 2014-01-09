\name{coSym}
\alias{coSym}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculates the symmetry of species in a network.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
coSym(x = "dependency network", zero.na = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{zero.na}{
%%     ~~Describe \code{zero.na} here~~
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
function (x = "dependency network", zero.na = TRUE) 
{
    out <- x * 0
    for (i in 1:nrow(x)) {
        for (j in 1:ncol(x)) {
            if (max(c(x[i, j], x[j, i])) == 0) {
            }
            else {
                out[i, j] <- abs(x[i, j] - x[j, i])/max(c(x[i, 
                  j], x[j, i]))
            }
        }
    }
    return(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
