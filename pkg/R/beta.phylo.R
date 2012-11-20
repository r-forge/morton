## Analyses of Rooney and Begley data, Fall 2012
## Andrew Hipp and Danielle Begley, 2012-10-15 ff

beta.phylo <- function(samp, tree, depth = c('basal', 'terminal'), sampleType = c('abundance', 'presence'), threshold = 0.001) {
  # data that work in picante will work in this
  # this is the 6th metric of Swenson, 2011, Dpw' if depth = 'basal'
  # ...         4th metric, Dnn' if depth = 'terminal'
  # ...         ignores abundance if sampleType = 'presence'
  if(sampleType[1] == 'presence') samp <- ifelse(samp > threshold, 1, 0) # turns abundances into presence-absence
  else samp <- samp / apply(samp, 1, sum) # makes these actual relative abundances in case they aren't
  patD <- cophenetic.phylo(tree)
  plots <- row.names(samp)
  out <- matrix(NA, length(plots), length(plots), dimnames = list(plots, plots))
  for(i in 1:(length(plots)-1)) {
    for(j in (i+1):length(plots)) {
	  iPlot <- names(which(samp[i, ] > threshold))
	  jPlot <- names(which(samp[j, ] > threshold))
	  iDeltaList <- lapply(iPlot, function(x) samp[i, x] * patD[x, jPlot])
	  jDeltaList <- lapply(jPlot, function(x) samp[j, x] * patD[x, iPlot])
	  if(depth[1] == 'basal') {
	    iDelta <- sum(unlist(sapply(iDeltaList, mean)))
		jDelta <- sum(unlist(sapply(jDeltaList, mean)))
		}
	  if(depth[1] == 'terminal') {	    
	    iDelta <- sum(unlist(sapply(iDeltaList, min)))
		jDelta <- sum(unlist(sapply(jDeltaList, min)))
		}
      out[j, i] <- mean(iDelta, jDelta)
	  }
	}
  return(out)
  }
