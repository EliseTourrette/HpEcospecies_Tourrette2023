## R script to calculate the FST between high / load load strains 
## and between siberia / indNAm strains
## use BIALLELIC SNPs only

setwd("/Users/elise/Desktop/HpGlobal_noHpGP/FST/")

## use the PopGenome package
library(PopGenome)

meta <- read.table("../METADATA/strains_SiberiaIndAm.csv", sep = ",", header = TRUE)
meta <- meta[meta$fs_pop %in% c("hspIndigenousNAmerica", "hspSiberia"),]

## do the FST on the biallelic core snps (CDS and intergenic)
## same set of snps as for the gwas
data <- readData("fasta_core", FAST = TRUE)
data <- set.populations(data, list(meta$strain[meta$level == "low"], meta$strain[meta$level == "high"]))

# data <- F_ST.stats(data)
# get.F_ST(data)
# data <-neutrality.stats(data, detail=TRUE) 
# data <- F_ST.stats.2(data)

## if want the per site F_ST value
## need to create sliding windows of size one
## (the different stats will be calculated for each window)
slide.data <- sliding.window.transform(data,1,1, type=2)
slide.data <- diversity.stats(slide.data)
str(slide.data@nuc.diversity.within)
x <- slide.data@nuc.diversity.within[,2]

slide.data <- F_ST.stats(slide.data, mode="nucleotide")
str(slide.data@nuc.F_ST.pairwise)
x <- slide.data@nuc.F_ST.pairwise

save(slide.data, file = "FST_highLow.RData")



