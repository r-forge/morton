## handling pyRAD data
## A Hipp, 2012-09-18
## updated Sept 2012 to accommodate new pyRAD format (GBS data)
## updated Dec 2012 to help with exporting RAD data and blasting

## right now just used for the IUPAC_CODE_MAP vector
## uncomment if not already installed:
#source("http://bioconductor.org/biocLite.R")
#biocLite("Biostrings")


consensus.pyRAD <- function(pyIn, ...) {
## use seqinr to generate a consensus sequence for each pyRAD locus
  require(seqinr)
  }
  
blast.pyRAD <- function(pyConsensus, ...) {}

read.pyRAD <- function(filename, reportInterval = 20000, breakLinesSeparate = FALSE, ...) {
## reads the all.aligned file out of pyRAD, parses into names, loci, sequences
## updated with breakLinesSeparate in Oct 2012 because pyRAD switched to single-line summaries at the end of each aligned file
## updated 2012-11-16 to keep breaklines intact
  message("Reading data...")
  dat <- readLines(filename, ...)
  dat.breakLines <- dat.consensusLines <- grep("//", dat, fixed = TRUE) # this is slow, but only ca. 1 sec for data runs of 10s of thousands
  dat.breakLines.vector <- dat[dat.breakLines] #added 2012-11-16; ignores possibility of separate breakLines
  message("Splitting data...")
  dat.split <- strsplit(dat, " {1,100}")
  dat.names <- as.factor(sapply(dat.split, function(x) x[1]))
  dat.seqs <- as.factor(sapply(dat.split, function(x) x[2]))
  # dat.seqs[dat.breakLines] <- dat.breakLines.vector # shoves the consensus seqs back into the sequence vector, assuming breakLinesSeparate = F
  dat.firstLocusLines <- c(1, (dat.breakLines[1:(length(dat.breakLines)-1)] + 1))
  if(breakLinesSeparate) {
    dat.lastLocusLines <- dat.breakLines - 2
    dat.consensusLines <- dat.breakLines - 1
	}
  else dat.lastLocusLines <- dat.breakLines - 1
  dat.locus.index <- character(length(dat))
  locusCounter <- 1
  message("Assigning locus number...")
  start.time <- Sys.time()
  ## crazy slow! before vectorization:
  #for (i in 1:length(dat)) {
  #  if(i-1 %in% dat.breakLines) locusCounter <- locusCounter + 1
  #	dat.locus.index[i] <- paste("locus.", locusCounter, sep = "")
  #	if(i / reportInterval - i %/% reportInterval == 0) {
  #	   message(paste('...', i, 'of', length(dat.locus.index), 
  #	   '-- Estimated time remaining =', ((Sys.time() - start.time) / i) * (length(dat.locus.index) - i), attr(Sys.time() - start.time, 'units')
  #	   ))
  #	   }
  #	}

  ## after some vectorization:
  for (i in 1:length(dat.firstLocusLines)) {
    dat.locus.index[dat.firstLocusLines[i]:dat.lastLocusLines[i]] <- names(dat.breakLines.vector)[i] <- paste("locus.", i, sep = "")
	if(i / reportInterval - i %/% reportInterval == 0) {
  	   message(paste('...', i, 'of', length(dat.firstLocusLines), 
 	   '-- Estimated time remaining =', ((Sys.time() - start.time) / i) * (length(dat.firstLocusLines) - i), attr(Sys.time() - start.time, 'units')
  	   ))
	}
  }
  dat.locus.index <- as.factor(dat.locus.index) # done only for memory considerations... slows things down in analysis unless you transform back to character
  out = list(tips = dat.names, seqs = dat.seqs, breaks = dat.breakLines, break.vectors = dat.breakLines.vector, cons = dat.consensusLines, locus.index = dat.locus.index)
  class(out) <- 'pyRAD.loci'
  return(out)
  }

summary.pyRAD.loci <- function(object, var.only = FALSE, ...) {
# Arguments:
#  object = a pyRAD.loci object
#  var.only = if T, only includes variable loci; as written, the function assumes a "*" if there is 
  reportInterval <- 2000 # this is just for screen reporting... only matters with really long files
  ## currently, locus.names includes a null (""), b/c the break lines have no locus name
  locus.names <- as.character(unique(object$locus.index)) # this slows things down by a factor of 2 or 3, but it seems to prevent a subscript-out-of-bounds error
  locus.names <- locus.names[locus.names != ""]
  ## REWRITE TO LOOK FOR BREAKLINES THAT HAVE * OR - IN THEM
  variable.loci <- locus.names[!is.na(object$seqs[object$breaks])] # note that this works with the pyRAD output we are currently using... should be checked
  if(var.only) locus.names <- variable.loci
  num.loci <- length(locus.names)
  tip.names <- as.character(unique(object$tips[-c(object$breaks, object$cons)]))
  message("Splitting tips by locus name...")
  tips.per.locus <- split(object$tips, object$locus.index)[locus.names]
  seqs.per.locus <- split(object$seqs, object$locus.index)[locus.names]
  num.inds.per.locus <- sapply(tips.per.locus, length)
  inds.mat <- matrix(NA, nrow = length(tip.names), ncol = num.loci, dimnames = list(tip.names, locus.names))
  message("Making tips matrix...")
  start.time <- Sys.time()
  ## is there some way to vectorize the following:
  for(i in seq(num.loci)) {
    temp <- try(tip.names %in% tips.per.locus[[locus.names[i]]])
	if(class(temp) != "try-error") inds.mat[ , locus.names[i]] <- temp
	else(message(paste("Error occurred with locus", locus.names[i])))
    if(i / reportInterval - i %/% reportInterval == 0) {
  	   message(paste('...', i, 'of', num.loci, 
 	   '-- Estimated time remaining =', ((Sys.time() - start.time) / i) * (num.loci - i), attr(Sys.time() - start.time, 'units')
  	   ))
	   }
	 }

  out <- list(num.loci = num.loci, tips.per.locus = tips.per.locus, break.vectors = object$break.vectors, seqs.per.locus = seqs.per.locus, num.inds.per.locus = num.inds.per.locus, variable.loci = variable.loci, inds.mat = inds.mat)
  class(out) <- 'summary.pyRAD.loci'
  out
  }

overlap.report <- function(dat, repPattern = "2.") {
## reports on how replicate individuals fall out in the loci
  if(class(dat) == 'summary.pyRAD.loci') dat <- dat$inds.mat
  reps <- grep(repPattern, row.names(dat), fixed = TRUE, value = TRUE)
  orig <- sapply(reps, function(x) strsplit(x, repPattern, fixed = TRUE)[[1]][2])
  out <- matrix(NA, nrow = length(orig), ncol = 3, dimnames = list(orig, c("Original", "Replicate", "Both")))
  for(i in 1:length(orig)) {
    message(paste("Doing names", orig[i]))
	out[orig[i], "Original"] <- sum(dat[orig[i],])
	out[orig[i], "Replicate"] <- sum(dat[reps[i],])
	out[orig[i], "Both"] <- sum(colSums(dat[c(orig[i], reps[i]), ]) == 2)
	}
  return(out)
  }

lengths.report <- function(dat, numtodo = 10, reportInterval = 2000, high.mem = TRUE) {
  if(class(dat) != 'pyRAD.loci') stop("This function runs on a pyRAD data object")
  last.lines <- dat$cons - 1
  num.loci <- length(last.lines)
  datSeqs <- as.character(dat$seqs)
  if(high.mem) block.lengths <- sapply(datSeqs[last.lines][1:ifelse(numtodo<1,num.loci,numtodo)], function(x) nchar(as.character(x)))
  else {
    block.lengths = integer(0)
	for(i in 1:num.loci) block.lengths = c(block.lengths, nchar(as.character(datSeqs[last.lines[i]])))
	if(i / reportInterval - i %/% reportInterval == 0) {
  	   message(paste('...', i, 'of', num.loci, 
 	   '-- Estimated time remaining =', ((Sys.time() - start.time) / i) * (num.loci - i), attr(Sys.time() - start.time, 'units')
  	   ))
	   }

	}
  return(block.lengths)
  }

filter.by <- function(dat, taxa) {
  ## returns just loci for which the requested taxa are all present
  if(class(dat) != 'summary.pyRAD.loci') stop("This function only works with summary.pyRAD.loci datatypes")
  dat.mat <- dat$inds.mat[taxa, ]
  return(names(which(apply(dat.mat, 2, sum) == length(taxa))))
  }

hybrid.test <- function(dat = test.18.v2.summary, parents = c('>2830D', '>2893G1'), 
                        f1s = c('>2830Dx2893G1A', '>2830Dx2893G1B','>2830Dx2893G1C','>2830Dx2893G1D', '>2830Dx2893G1E',
						         '>2893Gx2830D1A', '>2893Gx2830D1C'))								 
{
require(Biostrings)
  ## go through all loci and find for each (1) the parent genotypes, 
  ##(2) the expected F1 genotypes are and their ratios, 
  ##(3) the observed F1 genotypes and their ratios,
  ##(4) the percent F1 genotypes that are the expected genotypes
  
  variableSiteCharacters <- c("-", "*") ## change this if Deren rewrites pyRAD to use other characters
  if(class(dat) != 'summary.pyRAD.loci') stop('summary of pyRAD data needed here!')
  parentNumbers <- which(dimnames(dat$inds.mat)[[1]] %in% parents)
  f1Numbers <- which(dimnames(dat$inds.mat)[[1]] %in% f1s)
  loci.to.use <- which(apply(dat$inds.mat, 2, function(x) all(x[parentNumbers] == T) & any(x[f1Numbers] == T))) 
  loci.to.use.names <- dimnames(dat$inds.mat)[[2]][loci.to.use]
  loci.to.use.names <- loci.to.use.names[loci.to.use.names %in% dat$variable.loci]
  rm(loci.to.use)
  matsOut <- structure(vector('list', length(loci.to.use.names)), .Names = loci.to.use.names)
  colNames <- character(0)
  for(locusCounter in loci.to.use.names) {
	message(paste("Doing", locusCounter))
	seqsMat <- t(as.matrix(sapply(as.character(dat$seqs.per.locus[[locusCounter]]), function(x) strsplit(x, split = "")[[1]]), byrow = T))   
	dimnames(seqsMat)[[1]] <- dat$tips.per.locus[[locusCounter]]
	seqLength <- dim(seqsMat)[2]
	seqConsensus <- substr(dat$break.vectors[locusCounter], nchar(dat$break.vectors[locusCounter]) - seqLength + 1, nchar(dat$break.vectors[locusCounter]))
    variable.sites <- which(strsplit(seqConsensus, "")[[1]] %in% variableSiteCharacters)
	message(paste("Found", length(variable.sites), "variable sites"))
	seqsMat <- as.matrix(seqsMat[c(parents, dimnames(seqsMat)[[1]][dimnames(seqsMat)[[1]] %in% f1s]), variable.sites]) # only includes parents and children
	#browser()
	for (i in 1:dim(seqsMat)[2]) {
	  seqPossibilities <- character(0)
	  parent.sites <- strsplit(IUPAC_CODE_MAP[seqsMat[parents, i]], "")
	  for(j in 1:length(parent.sites[[1]])) {
	    for(k in 1:length(parent.sites[[2]])) {
		  seqPossibilities <- c(seqPossibilities, mergeIUPACLetters(paste(parent.sites[[1]][j], parent.sites[[2]][k], sep = "")))
		  }}
	  seqsMat <- cbind(seqsMat, c(NA, NA, seqsMat[3:dim(seqsMat)[1], i] %in% seqPossibilities))
	  }
	matsOut[[locusCounter]] <- seqsMat
	colNames <- c(colNames, paste(locusCounter, "_", seq(dim(seqsMat)[2] / 2), sep = ""))
	}
  message(paste("Columns in summary matrix:", length(colNames)))
  summaryMat <- matrix(NA, nrow = length(c("differ", f1s)), ncol = length(colNames), dimnames = list(c("Parents differ", f1s), colNames))
  colCounter <- 1
  for(i in 1:length(matsOut)) {
    message(paste('DOING MATRIX', i, 'OF', length(matsOut)))
	for(j in ((dim(matsOut[[i]])[2] / 2) + 1):dim(matsOut[[i]])[2]) {
	  for(k in 3:dim(matsOut[[i]])[1]) {
	    message(paste("Doing summary matrix column", colCounter))
		message(paste(" -- working on row", dimnames(matsOut[[i]])[[1]][k]))
		summaryMat[dimnames(matsOut[[i]])[[1]][k], colCounter] <- as.logical(matsOut[[i]][k,j])
		summaryMat[1, colCounter] <- matsOut[[i]][1,(j - dim(matsOut[[i]])[2] / 2)] != matsOut[[i]][2,(j - dim(matsOut[[i]])[2] / 2)]
		}# close k
		colCounter <- colCounter + 1	  
	  }# close j
	}# close i
  parents.differ <- summaryMat[1, ]
  percent.compatible.with.cross <- apply(summaryMat[, parents.differ], 1, mean, na.rm = TRUE)  
  out = list(matsOut = matsOut, summaryMat = summaryMat, percent.compatible.with.cross = percent.compatible.with.cross)
  return(out)
  }
  