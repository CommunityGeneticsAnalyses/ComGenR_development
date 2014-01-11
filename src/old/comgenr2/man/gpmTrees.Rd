\name{gpmTrees}
\alias{gpmTrees}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Foundation tree species genotype-phenotype mapping function for cgSim. 
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
gpmTrees(pheno = "phenotype of genotypes", reps = "replication")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pheno}{
%%     ~~Describe \code{pheno} here~~
}
  \item{reps}{
%%     ~~Describe \code{reps} here~~
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
function (pheno = "phenotype of genotypes", reps = "replication") 
{
    if (any(pheno == "phenotype of genotypes")) {
        pheno <- c(11, 12.5, 13.75, 16, 14, 15.25, 17.5, 16.5, 
            18.75, 21)
    }
    if (reps == "replication") {
        reps <- 5
    }
    trees <- list()
    for (i in 1:length(pheno)) {
        trees[[i]] <- rep(pheno[i], reps)
    }
    return(cbind(geno = gl(length(pheno), reps), pheno = unlist(trees)))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
