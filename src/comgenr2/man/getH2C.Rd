\name{getH2C}
\alias{getH2C}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculates the community heritability metric of Shuster et al. 2006. 
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getH2C(x = "NMDS axis", g = "grouping vector", sibs = "sibling proportion (1=clones)")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{g}{
%%     ~~Describe \code{g} here~~
}
  \item{sibs}{
%%     ~~Describe \code{sibs} here~~
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
function (x = "NMDS axis", g = "grouping vector", sibs = "sibling proportion (1=clones)") 
{
    g <- as.character(g)
    aov.tab <- matrix(unlist(summary(aov(x ~ g))), nrow = 2)
    if (sibs == "sibling proportion (1=clones)") {
        sibs <- 1
    }
    S <- length(unique(as.character(g)))
    ni <- table(g)
    n. <- sum(ni)
    if (all(ni == ni[1])) {
        k <- mean(ni)
    }
    else {
        k <- 1/(S - 1) * (n. - (sum(ni^2)/n.))
    }
    s2s <- (aov.tab[1, 3] - aov.tab[2, 3])/k
    s2w <- aov.tab[2, 3]
    s2total <- s2s + s2w
    H2C <- s2s/s2total
    t <- H2C * (1/sibs)
    if (all(ni == ni[1])) {
        SE <- ((2 * (((1 - t)^2) * (1 + (k - 1) * t)^2))/(k * 
            (k - 1) * (S - 1)))^(1/2)
    }
    else {
        SE <- ((2 * (n. - 1) * ((1 - t)^2) * (1 + (k - 1) * t)^2)/((k^2) * 
            (n. - S) * (S - 1)))^(1/2)
    }
    CI <- SE * 1.96
    return(c(H2C = H2C, CI = CI, SE = SE))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
