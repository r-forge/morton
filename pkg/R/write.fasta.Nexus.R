## functions for reading and writing files
## A Hipp, 23 sept 2010 (ahipp@mortonarb.org)

## Functions:
##  write.fasta.nexus

write.fasta.nexus <- function (sequences, otuNames = names(sequences), file.out, open = "w") 
{
    outfile <- file(description = file.out, open = open)
    writeLines("#NEXUS", outfile)
	writeLines(paste("[File written using write.fasta.nexus function in R,", date(), "]"), outfile)
	writeLines("\nBEGIN DATA;", outfile)
    writeLines(paste("   DIMENSIONS NTAX =", length(sequences), 
        " NCHAR =", length(sequences[[1]]), ";"), outfile)
    writeLines("   FORMAT MISSING=? DATATYPE=DNA GAP=- INTERLEAVE;", 
        outfile)
    writeLines("   OPTIONS GAPMODE=MISSING;", outfile)
    writeLines("\nMATRIX\n", outfile)
    write.oneseq <- function(DNAsequence, name) {
        writeLines(paste(name, c2s(DNAsequence), sep = "\t"), outfile)
    }
    if (!is.list(sequences)) {
        write.oneseq(sequence = sequences, name = otuNames)
    }
    else {
        n.seq <- length(sequences)
        sapply(seq_len(n.seq), function(x) write.oneseq(DNAsequence = as.character(sequences[[x]]), 
            name = otuNames[x]))
    }
    writeLines(";", outfile)
    writeLines("end;", outfile)
    close(outfile)
}

