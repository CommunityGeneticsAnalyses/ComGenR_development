\name{null.prune}
\alias{null.prune}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Pruning function for the Araujo et al. 2011 method.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
null.prune(a, b, std = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a}{
%%     ~~Describe \code{a} here~~
}
  \item{b}{
%%     ~~Describe \code{b} here~~
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
function (a, b, std = TRUE) 
{
    A <- length(a)
    pa <- sum(a)/A
    pb <- sum(b)/A
    paANDb <- pa * pb
    paORb <- pa + pb - 2 * pa * pb
    pNOTaORb <- 1 - pa - pb + pa * pb
    Ei <- A * paANDb
    Eii <- A * paORb
    Eiii <- A * pNOTaORb
    Vi <- A * paANDb * (1 - paANDb)
    Vii <- A * paORb * (1 - paORb)
    Viii <- A * pNOTaORb * (1 - pNOTaORb)
    ci.u <- Ei + 2 * Vi^(1/2)
    ci.l <- Ei - 2 * Vi^(1/2)
    Oi <- length(a[(a + b) == 2])
    if (Oi > ci.u | Oi < ci.l) {
        if (std) {
            out <- (Oi - Ei)/sqrt(Vi)
        }
        else {
            out <- Oi
        }
    }
    else {
        out <- 0
    }
    return(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
