neiD <- function(genotypes, n, x = NULL, y = NULL) {
## Unbiased estimate of D based on Nei 1978, Genetics 89, equation (6)
## Arguments:
##  genotypes = a data.frame with column 1 = "locus", column 2 = "allele", and columns 3:length(names(genotypes)) = individuals
##  n = a matrix with rows = loci, columns = individuals
## nomenclature from Nei:
##  p[i], q[i] = frequency of allele i in populations X and Y respectively
##  x[i], y[i] = corresponding sample allele frequences
##  G[x] = mean(sum(p[i]^2)) over all loci in the genome (for the population)
##  G[y] = mean(sum(q[i]^2)) ...
##  G[xy] = mean(sum(p[i]*q[i])) ...
##  J[x], J[y], J[xy] = corresponding sample means
##  D = -log(G[xy] / sqrt(G[x] * G[y])) -- this is the biased D
## to makes this unbiased -- Dhat -- use Ghat in the place of G, calculated as:
##  Ghat[x] = mean((2 * n[x] * J[x] - 1) / (2 * n[x] - 1)) over loci studied
##  Ghat[y] = mean((2 * n[y] * J[y] - 1) / (2 * n[y] - 1)) ...
##  Ghat[xy] = J[xy]
## thus, equation (6) is: Dhat = -log(Ghat[xy] / sqrt(Ghat[x] * Ghat[y]))
## Note that in equation (12), Nei gives an unbiased estimate of single locus genetic distance based on the minimum distance for the kth locus as:
##  d[k] = (2 * n[x] * sum(x[i]^2) - 1) / (2 * (2 * n[x] - 1)) + (2 * n[y] * sum(y[i]^2) - 1) / (2 * (2 * n[y] - 1)) - sum(x[i] * y[i])
##  Dhat-m = mean(d[k])
## Andrew Hipp (ahipp@mortonarb.org), May 2008

GhatXY <- function(x, y, genotypes, loci) {
  locusGhatXY <- numeric(length(loci)); names(locusGhatXY) <- loci
    for(l in loci) {
      Xalleles <- genotypes[[x]][genotypes$locus == l]
      Yalleles <- genotypes[[y]][genotypes$locus == l]
      # following line excludes any loci for which one or both taxa are without alleles
      locusGhatXY[l] <- ifelse((sum(Xalleles) != 0 && sum(Yalleles) != 0), yes = sum(Xalleles * Yalleles), no = NA)
      }
    GhatOut <- mean(locusGhatXY, na.rm = T)
    print(paste("GhatXY for populations",x,"and",y,"=",GhatOut))
    return(GhatOut)}

if(!identical(names(genotypes)[1:2], c("locus", "allele"))) stop("Column 1 must be locus names and column 2 allele names for this function to work properly.")
if(!identical(sort(unique(genotypes$locus)), sort(row.names(n)))) stop("Locus names must be identical in the genotypes and sample-size matrices.")
loci <- sort(unique(genotypes$locus))
populations <- names(genotypes)[3:length(names(genotypes))]
Ghat <- numeric(length(populations)); names(Ghat) <- populations

if(!identical(NULL, x)) {
  for(i in c(x,y)) { 
    popJ <- numeric(length(loci)); names(popJ) <- loci
    popGhat <- numeric(length(loci)); names(popGhat) <- loci
      for (l in loci) {
        popJ[l] <- sum(unlist(lapply(genotypes[[i]][genotypes$locus == l], function(x) x^2))) 
        popGhat[l] <- ifelse(popJ[l] != 0, yes = ((2 * n[l, i] * popJ[l]) - 1) / ((2 * n[l, i]) - 1), no = NA)
        }
      Ghat[i] <- mean(popGhat, na.rm = T) }
  return(-log(GhatXY(x, y, genotypes, loci) / sqrt(Ghat[x] * Ghat[y]))) }
  

## ASSIGN J and Ghat over loci
  # J <- numeric(length(populations)); names(J) <- populations
  for (x in populations) {
    popJ <- numeric(length(loci)); names(popJ) <- loci
    popGhat <- numeric(length(loci)); names(popGhat) <- loci
    for (l in loci) {
      popJ[l] <- sum(unlist(lapply(genotypes[[x]][genotypes$locus == l], function(x) x^2))) 
      popGhat[l] <- ifelse(popJ[l] != 0, yes = ((2 * n[l, x] * popJ[l]) - 1) / ((2 * n[l, x]) - 1), no = NA) }
      #popGhat[l] <- ((2 * n[l, x] * popJ[l]) - 1) / ((2 * n[l, x]) - 1)
    Ghat[x] <- mean(popGhat, na.rm = T) 
    print(paste("Ghat for population", x, "=", Ghat[x]))}
    
## make matrix
  neiMatrix <- matrix(NA, length(populations), length(populations), dimnames = list(populations, populations))
  for (x in 1:length(populations)) {
    for (y in x:length(populations)) {
      neiMatrix[populations[y], populations[x]] <- -log(GhatXY(populations[x], populations[y], genotypes, loci) / sqrt(Ghat[populations[x]] * Ghat[populations[y]]))}}
      neiMatrix[populations[x], populations[y]] <- exp(-neiMatrix[populations[y], populations[x]])
return(neiMatrix) }

