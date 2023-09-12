##4#######################################################
## create the file of the comparisons that need to be done 
## for the pairwise comparisons for the dN/dS analysis
#########################################################

setwd("~/Desktop/HpGlobal_noHpGP")

## rmk1: for dN/dS, strain1 vs strain2 and strain2 vs strain1 gives the same result
## do not use the same set of outliers strains as for the HpGP data as they were removed from this new dataset (WHY ?!)
## just use the hpAfrica2 strains present in the dataset (or at least the ones for which a population was assigned...)
pop <- read.table("METADATA/listStrains.txt", header = FALSE)$V1
pop2 <- read.table("METADATA/strains_pop.csv", sep = ",", header = TRUE)
pop2 <- pop2[pop2$Sample %in% pop,]
pop2 <- rbind(pop2, data.frame(Sample = pop[!(pop %in% pop2$Sample)], MergedPop = "NA", pop = "NA"))

outgroup <- pop2$Sample[pop2$pop == "hpAfrica2"]

## do not do the calculations for the outgroup strains
comp <- NULL
for(ipop in unique(pop2$pop)) {
    if(ipop != "hpAfrica2") {
        print(ipop)
        strain <- pop2$Sample[pop2$pop == ipop]
        strain <- strain[!is.na(strain)]
        tmp <- unlist(lapply(1:length(strain), function(x) paste(strain[x], outgroup)))
        comp <- c(comp, tmp)
    }
}
write.table(comp, file = "METADATA/pairwiseComparisons_outgroup.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

## calculate the dN/dS between the Hardy and Ubiquitous strains
## shall I use all of the Ubiquitous/Hardy strains or only a subset of them?
## use the same subset of strains that I were used for the pangenome/haplotype/functional (allele sharing) analyses
sampleStrains <- read.table("METADATA/strains_pangenome.txt", header = TRUE)
hardy <- sampleStrains$strain[sampleStrains$divergence == "high"]
ubi <- sampleStrains$strain[sampleStrains$divergence %in% c("low", "inter")]

comp <- unlist(lapply(hardy, lapply, paste, ubi))

write.table(comp, file = "METADATA/pairwiseComparisons_HardyUbiquitous.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)


#########################################################
#########################################################
