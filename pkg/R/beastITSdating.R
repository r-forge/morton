## Functions for estimating node ages +/- CI using BEAST output, calibrating using Kay et al. 2006 ITS rates
## To use:
##  1. Call 'beastIn' to input a Beast output file and create a "beastData" object
##  2. call 'sampleNodes' on the 'beastData' object, with subsample = number of trees to sample from your Beast run, CI = confidence interval, and burnin = number of trees to discard at the beginning of the run
##  3. call 'summary' on the 'sampleNodes' object to get useful summary data on your analysis
## Andrew Hipp, 8 Feb 2008
## ahipp@mortonarb.org

beastIn <- function(filename = NULL, nodeList = NULL) {
## Formats Beast log file for subsampling
  filename <- ifelse(identical(filename, NULL), choose.files(multi = F, caption = "Select the BEAST log file with which you want to work"), filename)
  beastData <- read.table(filename, header = T)
  if(identical(nodeList, NULL)) nodeList = select.list(names(beastData), preselect = NULL, multiple = T, title = "Select the nodes that interest you")
  outData <- list(dataMatrix = beastData, nodes = nodeList)
  class(outData) <- 'beastData'
  return(outData)
}
sampleNodes <- function(beastData, ratesVector = kay2006$rate[kay2006$lifeHistory == 'H'][1:10], subsample = 5000, CI = 0.95, burnin = round(dim(beastData$dataMatrix)[1] * 0.4)) {
## Samples depths of selected nodes at random from the Markov chain, returning for each (1) a vector of subsampled parameters, and (2) the same vector * a vector of subsampled rates
## NOTE: This analysis assumes that nodeDepth / rate is a meaningful value (e.g., that rate is in substitutions / site / year, and node depth is in substitution / site)
  if(class(beastData) != 'beastData') message('The data were not formatted using beastIn... I will try to analyze, but who knows what might happen?!')
  postBurnin = dim(beastData$dataMatrix)[1] - burnin
  if(postBurnin < subsample) stop(paste('There are only', postBurnin, 'elements in the post-burnin portion of your run...\nReduce your subsample size and rerun'))
  nodesVector <- beastData$nodes
  nodeDepthSubst <- matrix(nrow = subsample, ncol = length(nodesVector), dimnames = list(seq(subsample), nodesVector))
  nodeDepthYears <- nodeDepthSubst
  CIlistSubst <- vector('list', length(nodesVector)); names(CIlistSubst) <- nodesVector
  CIlistYears <- CIlistSubst
  meanSubst <- numeric(length(nodesVector)); names(meanSubst) <- nodesVector
  meanYear <- meanSubst
  CIelements = c(ceiling(subsample * (1 - CI) / 2), floor(subsample * (1 + CI) / 2))
  for(i in nodesVector) {
    nodeDepthSubst[ ,i] <- sample(beastData$dataMatrix[[i]][burnin : dim(beastData$dataMatrix)[1]], subsample, replace = F)
    nodeDepthYears[ ,i] <- nodeDepthSubst[ ,i] / sample(ratesVector, subsample, replace = T) 
    substSort <- sort(nodeDepthSubst[ ,i]); yearsSort = sort(nodeDepthYears[ ,i])
    CIlistSubst[[i]] <- c(substSort[CIelements[1]], substSort[CIelements[2]])
    CIlistYears[[i]] <- c(yearsSort[CIelements[1]], yearsSort[CIelements[2]])
    meanSubst[i] <- mean(nodeDepthSubst[ ,i])
    meanYear[i] <- mean(nodeDepthYears[ ,i]) }
  outData = list(substMatrix = nodeDepthSubst, yearsMatrix = nodeDepthYears, CI = list(subst = CIlistSubst, years = CIlistYears), mean = list(subst = meanSubst, year = meanYear), nodesVector = nodesVector, CI = CI)
  class(outData) = 'sampleNodes'
  return(outData)
}

summary.sampleNodes <- function(sampleNodesIn) {
  cat("-----------------------------------\n")
  cat("Summary data for sampleNodes object\n")
  cat("-----------------------------------\n\n")
  cat(paste("For each node, age is provided as mean (lower CI, upper CI) on a", CI, "confidence interval\n\n"))
  for(i in sampleNodesIn$nodesVector) {
    cat(paste(i,": ", round(sampleNodesIn$mean$year[i]), " (", round(sampleNodesIn$CI$years[[i]][1]), ", ", round(sampleNodesIn$CI$years[[i]][2]), ")\n", sep = ""))
    }}

print.sampleNodes <- summary.sampleNodes
