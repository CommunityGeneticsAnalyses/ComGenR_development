\name{cgSim}
\alias{cgSim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Community Genetics simulator using the methods of Shuster et al 2006.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
cgSim(tree.pheno = "tree phenotype matrix", insect = "community phenotype matrix", reps = 10, GG = 8, YY = 5, VeT = 8, Ve = 0.1, VeN = 15, K = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree.pheno}{
%%     ~~Describe \code{tree.pheno} here~~
}
  \item{insect}{
%%     ~~Describe \code{insect} here~~
}
  \item{reps}{
%%     ~~Describe \code{reps} here~~
}
  \item{GG}{
%%     ~~Describe \code{GG} here~~
}
  \item{YY}{
%%     ~~Describe \code{YY} here~~
}
  \item{VeT}{
%%     ~~Describe \code{VeT} here~~
}
  \item{Ve}{
%%     ~~Describe \code{Ve} here~~
}
  \item{VeN}{
%%     ~~Describe \code{VeN} here~~
}
  \item{K}{
%%     ~~Describe \code{K} here~~
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
function (tree.pheno = "tree phenotype matrix", insect = "community phenotype matrix", 
    reps = 10, GG = 8, YY = 5, VeT = 8, Ve = 0.1, VeN = 15, K = 100) 
{
    if (any(tree.pheno == "tree phenotype matrix")) {
        tree.pheno <- gpmTrees()
    }
    if (any(insect == "community phenotype matrix")) {
        insect <- gpmCom()
    }
    T <- nrow(tree.pheno)
    I <- nrow(insect)
    art_g <- matrix(NA, nrow = T, ncol = I)
    art_z <- matrix(NA, nrow = T, ncol = I)
    art_Vg <- matrix(NA, nrow = T, ncol = I)
    art_Vz <- matrix(NA, nrow = T, ncol = I)
    gen_load <- matrix(NA, nrow = T, ncol = I)
    dem_load <- matrix(NA, nrow = T, ncol = I)
    art_pop <- matrix(NA, nrow = T, ncol = I)
    out <- gg.list <- yy.list <- list()
    tic <- Sys.time()
    for (RR in 1:reps) {
        scores_XX <- matrix(NA, nrow = nrow(tree.pheno), ncol = GG)
        scores_XX[, 1] <- tree.pheno[, 2] + runif(T, 0, 1) * 
            VeT - VeT/2
        for (z in 2:GG) {
            scores_XX[1:T, z] <- scores_XX[1:T, 1]
        }
        trees <- scores_XX
        for (y in 1:YY) {
            for (z in 1:GG) {
                for (i in 1:T) {
                  for (j in 1:I) {
                    if (trees[i, z] < 2 * insect[j, 1]) {
                      art_g[i, j] <- 0
                    }
                    else if (trees[i, z] > 2 * insect[j, 2]) {
                      art_g[i, j] <- 1
                    }
                    else {
                      art_g[i, j] <- (trees[i, z] - 2 * insect[j, 
                        1])/(2 * insect[j, 2] - 2 * insect[j, 
                        1])
                    }
                    art_z[i, j] <- 2 * insect[j, 2] * art_g[i, 
                      j]^2 + 2 * art_g[i, j] * (1 - art_g[i, 
                      j]) * (insect[j, 1] + insect[j, 2]) + 2 * 
                      insect[j, 1] * (1 - art_g[i, j])^2 + runif(1) * 
                      Ve - Ve/2
                    art_Vg[i, j] <- 2 * art_g[i, j] * (1 - art_g[i, 
                      j])
                    art_Vz[i, j] <- art_Vg[i, j] * (insect[j, 
                      2] - insect[j, 1])^2 + Ve
                    gen_load[i, j] <- 0.5 * (7.924e-05 * 2.511886^(z - 
                      1)) * (art_z[i, j] - trees[i, z])^2
                    dem_load[i, j] <- 0.5 * (7.924e-05 * 2.511886^(z - 
                      1)) * (art_Vz[i, j])
                    art_pop[i, j] <- K * (1 - gen_load[i, j] - 
                      dem_load[i, j]) + runif(1) * VeN * (y - 
                      1) - VeN * (y - 1)/2
                    if (art_pop[i, j] < 0) {
                      art_pop[i, j] <- runif(1) * 3
                    }
                    else {
                      art_pop[i, j]
                    }
                  }
                }
                print(paste(RR, y, z, sep = " "))
                gg.list[[z]] <- art_pop
                names(gg.list)[z] <- paste(RR, y, z, sep = "_")
            }
            yy.list[[y]] <- gg.list
        }
        toc <- Sys.time()
        out[[RR]] <- yy.list
    }
    return(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
