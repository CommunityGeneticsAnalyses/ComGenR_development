\name{dep.net}
\alias{dep.net}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
dep.net(x = "species in cols", zero.na = TRUE, prune = TRUE, diag.zero = TRUE, pos = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{zero.na}{
%%     ~~Describe \code{zero.na} here~~
}
  \item{prune}{
%%     ~~Describe \code{prune} here~~
}
  \item{diag.zero}{
%%     ~~Describe \code{diag.zero} here~~
}
  \item{pos}{
%%     ~~Describe \code{pos} here~~
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
function (x = "species in cols", zero.na = TRUE, prune = TRUE, 
    diag.zero = TRUE, pos = TRUE) 
{
    out <- matrix(NA, nrow = ncol(x), ncol = ncol(x))
    for (i in 1:ncol(x)) {
        for (j in 1:ncol(x)) {
            if (pos) {
                out[i, j] <- calcDepend(x[, i], x[, j])
            }
            else {
                out[i, j] <- negDepend(x[, i], x[, j])
            }
        }
    }
    if (prune) {
        out.rm <- co.net(x, diag.zero = diag.zero)
        out[out.rm == 0] <- 0
    }
    else {
    }
    if (diag.zero) {
        diag(out) <- 0
    }
    rownames(out) <- colnames(out) <- colnames(x)
    if (zero.na) {
        out[is.na(out)] <- 0
    }
    return(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
