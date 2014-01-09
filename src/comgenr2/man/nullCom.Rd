\name{nullCom}
\alias{nullCom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate a set of null communities based on functions in the vegan
package. 
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
nullCom(com, method = "r1", nits = 5000, burn = 500, thin = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{com}{
%%     ~~Describe \code{com} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
}
  \item{nits}{
%%     ~~Describe \code{nits} here~~
}
  \item{burn}{
%%     ~~Describe \code{burn} here~~
}
  \item{thin}{
%%     ~~Describe \code{thin} here~~
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
function (com, method = "r1", nits = 5000, burn = 500, thin = 10) 
{
    com[com != 0] <- 1
    for (i in 1:burn) {
        post.burn <- commsimulator(x = com, method = method, 
            thin = thin)
    }
    out <- list()
    for (i in 1:nits) {
        out[[i]] <- commsimulator(x = post.burn, method = method, 
            thin = thin)
    }
    return(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
