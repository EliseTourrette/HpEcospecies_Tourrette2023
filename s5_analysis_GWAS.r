## R script to do the GWAS analysis on the hspSiberia / hspIndigenousNAmerica strains with low/high dNdS (as a binary phenotype)
## use the package bugwas v0.0.0.9000

## the GWAS is done on SNPs

setwd("/Users/elise/Desktop/HpGlobal_noHpGP/GWAS")

###################################
## PREPROCESSING
###################################

## keep only the strains from Siberia and IndigenousNAmerica
## and determine the phenotype based on their dn/ds level, can also be determined by the PCA results 
meta <- read.table("../METADATA/strains_SiberiaIndAm.csv", sep = ",", header = TRUE)
meta$pheno <- 0
meta$pheno[meta$level == "high"] <- 1
pheno <- data.frame(ID = meta$strain[meta$fs_pop %in% c("hspSiberia", "hspIndigenousNAmerica")], pheno = meta$pheno[meta$fs_pop %in% c("hspSiberia", "hspIndigenousNAmerica")])
write.table(pheno, file = "input/pheno.txt", quote = FALSE, row.names = FALSE, sep = "\t")

## look at the phylogenetic tree
library(ggtree)
library(tidyverse)

setwd("/Users/elise/Desktop/HpGlobal/GWAS")

phenotype <- read.table("input/pheno.txt", header = TRUE)
highdNdS <- phenotype$ID[phenotype$pheno == 1]

fs_pop <- read.table("../METADATA/strains_pop.txt", header = TRUE)

tree <- read.tree("input/tree")
#tree <- read.tree("../TREE/tree") ## with IndSAm strains added
#highdNdS <- meta$strain[meta$level == "high"]

dNdS <- data.frame(strain = tree$tip.label, level = "low")
dNdS$level[dNdS$strain %in% highdNdS] <- "high"
dNdS$pop <- "hspSiberia"
dNdS$pop[dNdS$strain %in% fs_pop$Sample[fs_pop$pop == "hspIndigenousNAmerica"]] <- "hspIndigenousNAmerica"
dNdS$pop[dNdS$strain %in% fs_pop$Sample[fs_pop$pop == "hspIndigenousSAmerica"]] <- "hspIndigenousSAmerica"

p <- ggtree(tree, layout = "circular") + 
    geom_treescale() 
    
p %<+% dNdS +
    geom_tiplab(aes(col = interaction(level, pop)), align = TRUE, linesize = .5, size = 2.5)
ggsave("../PLOT/TREE/TREE_GWAS.jpeg", height = 14, width = 14)

###################################
## GWAS proper
###################################

## use the R package bugeas to do the GWAS
## this package make a call to the GEMMA software (v0.93)
## while additionally taking into account the population structure (use of the phylogenetic tree)
setwd("/Users/elise/Desktop/HpGlobal_noHpGP/GWAS/")

library(bugwas)

gen <- "input/geno_biallelic_SNP.txt"
pheno <- "input/pheno.txt"
phylo <- "input/tree"
prefix <- "bugwas"
gem.path <- "/home/elise/testSoftware/bugwas/gemma/gemma.0.93b"

data <- lin_loc(gen = gen, pheno = pheno, phylo = phylo, prefix = prefix, gem.path = gem.path, creatingAllPlots = FALSE)

save(data, file = "data_GWAS.RData")

all_plots(biallelic = data$biallelic, triallelic = data$triallelic, genVars = data$genVars, treeInfo = data$treeInfo, config = data$config)


###################################
## POSTPROCESSING
###################################

setwd("/Users/elise/Desktop/HpGlobal_noHpGP/GWAS/")
library(ggplot2)

## look at the function of the significative SNPs
library(readr)
gff <- data.frame(read_delim("../DATA/26695.gff", delim = "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE, skip = 2))
gffNCDS <- gff[which(gff$X3 != "CDS"),]
gff <- gff[which(gff$X3 == "CDS"),]

## look at the manhattan plot 
## and at the positions of the significant SNPs
## color them based on the PC with which each SNP correlate with
## note that both data0 (GWAS with only biallelic SNPs) and data (GWAS with all SNPs but result only for biallelic SNPs) are the same
data <- read.table("bugwas_biallelic_lmmout_allSNPs.txt", header = TRUE)



## p_lrt is the corrected p-value?
plot(data$ps, data$negLog10)
abline(h = 10, col = "green")
## question: has correction for multiple testing already been done for the p-values?

data$type <- 0
#data2$type <- 1
#data <- data.frame(position = c(data$ps, data2$ps), logp = c(data$negLog10, data2$negLog10), type = c(data$type, data2$type))
data <- data.frame(position = data$ps, logp = data$negLog10, type = 0)

## plot the different regions (50 kb) of the genome separately
## and add the position of the corresponding genes
for(i in seq(from = 0, to = 1700, by = 5)) {
    xmin <- i * 10000
    xmax <- xmin + 50000
    tmp <- data[data$logp != Inf & data$position > xmin & data$position <= xmax,]
    ## only save the figure if the region contains significant SNPs
    if(sum(tmp$logp > 5)) {
        CDS <- gff[which(gff$X4 > xmin & gff$X5 < xmax),]
        ggplot(tmp) +
            geom_point(aes(x = position/1000, y = logp, col = factor(type))) +
            geom_hline(yintercept = 5, col = "green") +
            geom_segment(data = CDS, aes(x = X4/1000, y = rep(-2, nrow(CDS)), xend = X5/1000, yend = rep(-2, nrow(CDS)))) +
            labs(x = "position (kb)", y = "-log(p-val)", col = "") +
            scale_colour_manual(values = c("0" = "black", "2" = "green"), labels = c("biallelic", "significance thereshold"))
            #scale_colour_manual(values = c("0" = "black", "1" = "red", "2" = "green"), labels = c("biallelic", "tri,tetra-allelic", "significance thereshold"))
        ggsave(paste0("../PLOT/GWAS/regions/", xmin, "_", xmax, ".png"), width = 19.6, height = 7.64)
    }
}


sig <- data[data$logp > 5 & data$logp != Inf,]
## if consider that a CDS is 5kb maximum
## remove all SNPs that are not in at most 1kb from their neighbors
## suppose a concentration of SNPs if the whole CDS is significant
## will need to check afterwards with the plots
delta <- c(0, abs(sig$position[-nrow(sig)] - sig$position[-1]))
sig <- sig[delta < 1000,]
sigCDS <- NULL
for(x in sig$position) {
    sigCDS <- rbind(sigCDS, gff[which(gff$X4 < x & gff$X5 > x),])
}
sigCDS <- sigCDS[!duplicated(sigCDS),]

sigCDS$product <- unlist(lapply(sigCDS$X9, function(x) sapply(strsplit(x, split = "product="), function(y) y[2])))

write.table(sigCDS[order(sigCDS$X4),c(4,5,10)], file = "significantCDS.csv", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ",")
write.table(sig, file = "significantSNPs.txt", quote = FALSE, row.names = FALSE, col.names = TRUE)






