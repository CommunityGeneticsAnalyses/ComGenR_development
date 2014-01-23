simSpp <- function (n = "number of species", het.values = c(5, 21), allelic.range = c(0,3)){
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
