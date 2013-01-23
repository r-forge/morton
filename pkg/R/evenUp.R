## make a tree ultrametric by tapping down the tips
## useful if you have manually added tips and need to clean things up a bit
## Andrew Hipp, 20 Nov 2012

evenUp.phylo <- function(tr, absDepth = NA, method = c('median', 'mean', 'min', 'max')) {
  require(phytools)
  if(class(tr) != "phylo") stop('phylo object required')
  tipLabels <- which(tr$edge[, 2] <= length(tr$tip.label))
  tipHeights <- nodeHeights(tr)[tipLabels, 2]
  if(!is.na(absDepth)) tipsOverage <- tipHeights - absDepth
  else tipsOverage <- tipHeights - switch(method[1], 
                                 median = median(tipHeights),
								 mean = mean(tipHeights),
					             min = min(tipHeights),
					             max = max(tipHeights)
					             )
  tr$edge.length[tipLabels] <- tr$edge.length[tipLabels] - tipsOverage
  return(tr)
  }