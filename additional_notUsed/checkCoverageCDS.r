## how many codon aare present at 100% in each of the different CDS
## to put on the github

setwd("~/Desktop/HpGlobal/PREPROCESSING/CDS")
for(i in 21:1576) {
    print(i)
    data <- read.table(paste0(i, ".aln"))
    
    strain <- data$V1[seq(from = 1, to = nrow(data), by = 2)]
    cds <- data$V1[seq(from = 2, to = nrow(data), by = 2)]
    
    cds <- strsplit(cds, split = "")
    cds2 <- NULL
    for(j in seq(from = 1, to = length(cds[[1]]), by = 3)) {
        cds2 <- rbind(cds2, sapply(cds, function(x) paste(x[j:(j+2)], collapse = ""))) 
    }
    
    
    ## strains that have the whole gene missing
    imiss <- 0
    jmiss <- NULL
    for(j in 1:length(strain)) {
        if(length(unique(cds[[j]])) == 1) if(unique(cds[[j]]) == '-') {imiss <- imiss + 1; jmiss <- c(jmiss, j)}
    }
    
    ## remove the strains for which the whole gene was missing
    if(length(jmiss) > 0) cds2 <- cds2[,-jmiss]
    
    codons <- apply(cds2, 1, unique)
    ## for each codons, check whether at least one strain has a gap
    codmiss <- 0
    for(j in 1:length(codons)) {
        cod <- codons[[j]]
        cod <- unique(unlist(strsplit(cod, split = "")))
        if('-' %in% cod) codmiss <- codmiss + 1
    }
    
    nbcodon <- data.frame(gene = i, nbCodons = length(codons), nbCodonsWithGap = codmiss, nbStrainMissingGene = imiss)
    if(i == 0) write.table(nbcodon, file = "../coverageCDS.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)
    if(i > 0) write.table(nbcodon, file = "../coverageCDS.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, append = TRUE)
    
}


