## make a tree from a taxonomic hierarchy
## ahipp@mortonarb.org, 2013-11-25

example.tax <- function() {
  require(ape)
  tree = read.delim(file.choose(), as.is = TRUE)
  trelease.tree <- tax2tree(tree)
  plot(trelease.tree)
  out <- dist.nodes.2(tree)
  return(out)
  }

tax2tree <- function(taxonomy) {
  nEdge = dim(taxonomy)[1]
  ntips <- sum(taxonomy$tip)
  tr <- list(edge = matrix(NA, nrow = nEdge, ncol = 2),
			 edge.length = rep(1, nEdge), 
			 tip.label = taxonomy$child[taxonomy$tip]
			 )
  tr$edge[taxonomy$tip, 2] <- seq(ntips)
  tr$edge[!taxonomy$tip, 2] <- (ntips + 2):(nEdge + 1)
  tr$edge[, 1] <- tr$edge[match(taxonomy$parent, taxonomy$child), 2]
  tr$edge[is.na(tr$edge[, 1]), 1] <- ntips + 1
  tr$node.label = taxonomy$parent[match((ntips + 1):(nEdge + 1), tr$edge[, 1])]
  tr$Nnode <- 1 + (nEdge - ntips)
  class(tr) <- 'phylo'
  return(tr)
  }

dist.nodes.2 <- function(x, drop.pattern = 'n[0123456789]') {
  out <- dist.nodes(x)
  out.names <- c(x$tip.label, x$node.label)
  dimnames(out) <- list(out.names, out.names)
  out <- out[-grep(drop.pattern, out.names), -grep(drop.pattern, out.names)]
  return(out)
  }