\name{co.net}
\alias{co.net}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Co-occurrence network function.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
co.net(x = "species in cols", diag.zero = TRUE, std = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{diag.zero}{
%%     ~~Describe \code{diag.zero} here~~
}
  \item{std}{
%%     ~~Describe \code{std} here~~
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
function (x = "species in cols", diag.zero = TRUE, std = TRUE) 
{
    x[x != 0] <- 1
    out <- matrix(NA, nrow = ncol(x), ncol = ncol(x))
    rownames(out) <- colnames(out) <- colnames(x)
    for (j in 1:ncol(x)) {
        for (k in 1:ncol(x)) {
            out[j, k] <- null.prune(x[, j], x[, k], std = std)
        }
    }
    if (diag.zero) {
        diag(out) <- 0
    }
    else {
    }
    return(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
