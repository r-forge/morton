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

kay2006 <-
## Table from Kay et al. 2006, BMC Evolutionary Biology 2006, 6:36 doi:10.1186/1471-2148-6-36
## Citations not included, but available on original table
structure(list(taxon = structure(c(16L, 29L, 19L, 25L, 4L, 8L, 
3L, 10L, 26L, 2L, 12L, 20L, 17L, 21L, 1L, 11L, 28L, 23L, 18L, 
5L, 6L, 9L, 22L, 14L, 7L, 15L, 24L, 27L, 13L), .Label = c("Adansonia", 
"Aesculus", "Alnus", "AraliaSectDimorphanthus", "Astragalus", 
"Cucurbitoideae", "Dendroseris", "Echium", "Ehrharta", "Empetraceae", 
"Eupatorium", "Gaertnera", "GentianaSectCiminalis", "Gentianella", 
"Gossypium", "Hamamelis", "Inga", "Lupinus", "Nothofagus", "Ormocarpum", 
"Phylica", "Plantago", "RobinioidLegumes", "Robinsonia", "Salicaceae", 
"Saxifraga", "Soldanella", "Tarweeds-HawaiianSilverswords", "Winteraceae"
), class = "factor"), family = structure(c(10L, 21L, 13L, 19L, 
1L, 5L, 3L, 7L, 20L, 11L, 18L, 8L, 8L, 17L, 4L, 2L, 2L, 8L, 8L, 
8L, 6L, 15L, 14L, 9L, 2L, 12L, 2L, 16L, 9L), .Label = c("Araliaceae", 
"Asteraceae", "Betulaceae", "Bombacaceae", "Boraginaceae", "Cucurbitaceae", 
"Empetraceae", "Fabaceae", "Gentianaceae", "Hamamelidaceae", 
"Hippocastanaceae", "Malvaceae", "Nothofagaceae", "Plantaginaceae", 
"Poaceae", "Primulaceae", "Rhamnaceae", "Rubiaceae", "Salicaceae", 
"Saxifragaceae", "Winteraceae"), class = "factor"), lifeHistory = structure(c(2L, 
2L, 2L, 2L, 2L, 2L, 2L, 2L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 1L, 2L, 
2L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 1L, 2L, 1L, 1L), .Label = c("H", 
"W"), class = "factor"), clockTest = structure(c(1L, 1L, 1L, 
1L, 1L, 2L, 3L, 2L, 1L, 2L, 2L, 3L, 3L, 3L, 2L, 2L, 2L, 3L, 1L, 
1L, 1L, 3L, 3L, 2L, 2L, 1L, 2L, 2L, 2L), .Label = c("na", "Passed", 
"Rejected"), class = "factor"), calibrationType = structure(c(4L, 
4L, 2L, 2L, 4L, 4L, 2L, 2L, 9L, 2L, 2L, 4L, 4L, 4L, 2L, 7L, 1L, 
2L, 2L, 2L, 2L, 8L, 4L, 5L, 6L, 3L, 4L, 2L, 4L), .Label = c("C", 
"F", "F&M-cpDNA", "G", "G&F", "G&M-cpDNA", "M-ndhF", "M-rbcL&F", 
"na"), class = "factor"), calibrationDate = c(8.5, 65, 83, 50, 
12, 20, 70, 37, 5.5, 65, 54, 35, 3.5, 2, 47, 14.8, 15, 39.4, 
60, 35, 40, 41, 0.6, 3, 3.3, 8.5, 4, 23.3, 0.1), rate = c(3.8e-10, 
4.5e-10, 5e-10, 6e-10, 1.07e-09, 1.1e-09, 1.1e-09, 1.44e-09, 
1.72e-09, 1.72e-09, 1.99e-09, 2e-09, 2.34e-09, 2.44e-09, 2.48e-09, 
2.51e-09, 3e-09, 3.3e-09, 3.46e-09, 3.5e-09, 3.62e-09, 3.81e-09, 
4.27e-09, 4.52e-09, 5e-09, 5.5e-09, 7.83e-09, 8.34e-09, 1.9e-08
), outlier = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, 
FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, 
FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, 
FALSE, FALSE, FALSE, TRUE)), .Names = c("taxon", "family", "lifeHistory", 
"clockTest", "calibrationType", "calibrationDate", "rate", "outlier"
), class = "data.frame", row.names = c(NA, -29L))