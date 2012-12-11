write.nexus = function (..., file = "", translate = TRUE, original.data = TRUE) 
## from ape, but with tree names from the names of the list
{
    obj <- list(...)
    if (length(obj) == 1) {
        if (class(obj[[1]]) == "phylo") 
            ntree <- 1
        else {
            obj <- obj[[1]]
            ntree <- length(obj)
        }
    }
    else ntree <- length(obj)
    cat("#NEXUS\n", file = file)
    cat(paste("[R-package APE, ", date(), "]\n\n", sep = ""), 
        file = file, append = TRUE)
    if (original.data) {
        if (!is.null(attr(obj[[1]], "origin"))) {
            if (!file.exists(attr(obj[[1]], "origin"))) {
                warning(paste("the file", attr(obj[[1]], "origin"), 
                  "cannot be found,\nthe original data won't be written with the tree."))
                original.data <- FALSE
            }
            else {
                ORI <- scan(file = attr(obj[[1]], "origin"), 
                  what = character(), sep = "\n", skip = 1)
                start <- grep("BEGIN TAXA;", ORI)
                ORI <- ORI[-(1:(start - 1))]
                ORI <- gsub("ENDBLOCK;", "END;", ORI)
                endblock <- grep("END;", ORI)
                start <- grep("BEGIN TREES;", ORI)
                end <- endblock[endblock > start][1]
                cat(ORI[1:(start - 1)], file = file, append = TRUE, 
                  sep = "\n")
                ORI <- ORI[-(1:end)]
            }
        }
        else original.data <- FALSE
    }
    N <- length(obj[[1]]$tip.label)
    if (!original.data) {
        cat("BEGIN TAXA;\n", file = file, append = TRUE)
        cat(paste("\tDIMENSIONS NTAX = ", N, ";\n", sep = ""), 
            file = file, append = TRUE)
        cat("\tTAXLABELS\n", file = file, append = TRUE)
        cat(paste("\t\t", obj[[1]]$tip.label, sep = ""), sep = "\n", 
            file = file, append = TRUE)
        cat("\t;\n", file = file, append = TRUE)
        cat("END;\n", file = file, append = TRUE)
    }
    cat("BEGIN TREES;\n", file = file, append = TRUE)
    if (translate) {
        cat("\tTRANSLATE\n", file = file, append = TRUE)
        tmp <- checkLabel(obj[[1]]$tip.label)
        X <- paste("\t\t", 1:N, "\t", tmp, ",", sep = "")
        X[length(X)] <- gsub(",", "", X[length(X)])
        cat(X, file = file, append = TRUE, sep = "\n")
        cat("\t;\n", file = file, append = TRUE)
        token <- as.character(1:N)
        names(token) <- obj[[1]]$tip.label
        obj[[1]]$tip.label <- token
        if (ntree > 1) {
            for (i in 2:ntree) obj[[i]]$tip.label <- token[obj[[i]]$tip.label]
            class(obj) <- NULL
        }
    }
    else {
        for (i in 1:ntree) obj[[i]]$tip.label <- checkLabel(obj[[i]]$tip.label)
    }
    for (i in 1:ntree) {
        if (class(obj[[i]]) != "phylo") 
            next
        if(identical(names(obj)[i], NULL)) treeTitle <- "UNTITLED" else treeTitle = names(obj)[i]
		if (is.rooted(obj[[i]])) 
            cat(paste("\tTREE *", treeTitle, "= [&R] "), file = file, append = TRUE)
        else cat(paste("\tTREE *", treeTitle, "= [&U] "), file = file, append = TRUE)
        cat(write.tree(obj[[i]], file = ""), "\n", sep = "", 
            file = file, append = TRUE)
    }
    cat("END;\n", file = file, append = TRUE)
    if (original.data) 
        cat(ORI, file = file, append = TRUE, sep = "\n")
}
