## reads raw ABI output and change to 1/0 file
## andrew hipp, 26 feb 09
## ahipp@mortonarb.org

aflpIn <- function(filename = NULL, interactive = T) {
if (interactive) dat <- read.table(file.choose(), as.is = T, sep = "", fill = T, col.names = c('inds', 'loci'))
  inds <- as.character(dat$inds)
  loci <- as.character(dat$loci)
  outdata <- matrix(data = 0, 
                    nrow = length(unique(inds[!is.na(inds)])), 
                    ncol = length(unique(loci[!is.na(loci)])), 
                    dimnames = list(sort(unique(inds[!is.na(inds)])), sort(unique(loci[!is.na(loci)])))
                    )
  for(i in 1:length(inds)) {
    if (is.na(loci[i])) next
    else outdata[inds[i], loci[i]] <- 1
    }
  return(outdata)
  }