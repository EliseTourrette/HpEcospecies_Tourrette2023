## Rscript for the post-processing of the FST and GWAS results
## in particular, apply a criterion based on the FST and GWAS results for the sites
## look only at the sites that were significant in both analysis

setwd("/Users/elise/Desktop/HpGlobal_noHpGP/HAPLOTYPES/")

library(reshape2)
library(ggplot2)
library(cowplot)
library(gridExtra)
library(ggpubr) ## for ggarrange

## FUNCTIONS
#########################################################

## process the dataset with the different alleles for the significant sites
## consider that the data from data_pval10FST0.9.RData has already been loaded into the session

## add the populations to the dataset
addPop <- function(data2) {
    data2$pop <- NA
    pop <- read.table("/Users/elise/Desktop/HpGlobal_noHpGP/METADATA/strains_pop.csv", sep = ",", header = TRUE)
    for(i in unique(pop$MergedPop)) {
        strain_pop <- pop$Sample[pop$MergedPop == i]
        data2$pop[data2$strain %in% strain_pop] <- i
    }
    for(i in unique(pop$pop)) {
        strain_pop <- pop$Sample[pop$pop == i]
        data2$pop[data2$strain %in% strain_pop] <- i
    }
    ## and the load level for some of the pop
    meta <- read.table("/Users/elise/Desktop/HpGlobal_noHpGP/METADATA/strains_SiberiaIndAm.csv", sep = ",", header = TRUE)
    i <- meta$strain[meta$level == "high"]
    data2$pop[data2$strain %in% i] <- paste(data2$pop[data2$strain %in% i],"high", sep = "_")
    i <- meta$strain[meta$level == "low"]
    data2$pop[data2$strain %in% i] <- paste(data2$pop[data2$strain %in% i],"low", sep = "_")
    
    return(data2)
}

## manipulate the populations from the dataset
## assign specific populations for the H. acinonychis and monkey strains
popData <- function(data2) {
    data2$pop[data2$strain %in% c("HEL_AA0408AA_AS", "HEL_BA9262AA_AS")] <- "Hacinonychis"
    data2$pop[data2$strain %in% c("HEL_BA9707AA_AS", "HEL_BA9303AA_AS")] <- "Monkey"
    
    data2$popLevel <- "inter"
    data2$popLevel[data2$pop %in% c("hspSiberia_low", "hspIndigenousSAmerica_low", "hspIndigenousNAmerica_low")] <- "low"
    data2$popLevel[data2$pop %in% c("hspSiberia_high", "hspIndigenousSAmerica_high", "hspIndigenousNAmerica_high")] <- "high"
    data2$popLevel[data2$pop %in% c("Hacinonychis", "Monkey")] <- "high_inter"
    for(i in c("hspSiberia", "hspIndigenousSAmerica", "hspIndigenousNAmerica")) data2$pop[data2$pop %in% paste0(i, c("_low", "_high"))] <- i
    data2 <- data2[order(data2$popLevel, data2$pop, data2$strain, data2$pos),]
    
    data2$strain <- factor(data2$strain, levels = unique(data2$strain))
    
    return(data2)
}

## add the CDS limits
## to plot them on the haplotype plot
addCDS <- function(data2, gff) {
    ## get the limits of the cds
    cdsLim <- c(gff$V4, gff$V5)
    cdsLim <- cdsLim[cdsLim >= min(data2$pos) & cdsLim <= max(data2$pos)]
    
    ## want to only keep the CDS limits
    ## that actually corresponds to sites 
    cdsLim <- cdsLim[round(cdsLim/1000) %in% unique(round(unique(data2$pos/1000)))]
    
    uniqueSet <- data2[data2$pos == min(data2$pos),]
    ## add the limits of the CDS to the dataset
    ## with nuc = NA and allele = limit (color = black)
    for(i in cdsLim) { 
        print(i)
        data2 <- rbind(data2, data.frame(strain = uniqueSet$strain, pos = i, nuc = "NA", pop = uniqueSet$pop, allele = "limit", popLevel = uniqueSet$popLevel))
    }
    
    return(data2)
}

#########################################################


## LOAD THE DATA
#########################################################

## !! need to load the FST results
## from the file ../FST/FST_highLow.RData

setwd("/Users/elise/Desktop/HpGlobal_noHpGP")

load("FST/FST_highLow.RData") ## slide.data
x <- slide.data@nuc.F_ST.pairwise

gwas <- read.table("GWAS/bugwas_biallelic_lmmout_allSNPs.txt", header = TRUE)
gwas <- data.frame(pos = gwas$ps, pval = gwas$negLog10)

tot <- data.frame(pos = gwas$pos, FST = x[gwas$pos])
tot <- tot[order(tot$pos),]
# sum(gwas$pos[order(gwas$pos)] != tot$pos[order(tot$pos)])
tot$pval <- gwas$pval[order(gwas$pos)]
tot$col <- "black"
tot$col[tot$FST > 0.9 & tot$pval > 10] <- "red"
tot$col[tot$FST < 0.5 & tot$pval < 10] <- "green"
plot(tot$pval, tot$FST, col = tot$col)

write.table(tot, file = "HAPLOTYPES/FST_GWASpval_lowHigh.txt", quote = FALSE, row.names = FALSE)

## also save the positions of the SNPs with FST > 0.9 and -log(pval) > 10 (GWAS)
write.table(tot$pos[tot$FST > 0.9 & tot$pval > 10], file = "HAPLOTYPES/pos_pval10FST0.9.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
## SNPs with FST < 0.5 and -log(pval) < 10 (GWAS)
write.table(tot$pos[tot$FST < 0.5 & tot$pval < 10], file = "HAPLOTYPES/pos_nonSigSites_pval10FST0.5.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
## other SNPs
write.table(tot$pos[tot$col == "black"], file = "HAPLOTYPES/pos_otherSNPs.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)

#########################################################

## GET THE ALLELE TYPE AT SIGNIFICANT POSITIONS FOR ALL STRAINS
#########################################################
## three criteria:
## FST > 0.9
## -log(p) > 10 (GWAS)
## and sites with both threshold

set <- "pval10FST0.9"

## the file sites_pval10FST0.9.txt was generated via the python script s7bis_haplotypeSigSitesDivergence.py
# system("python3 s7bis_haplotypeSigSitesDivergence.py")
data <- read.table(paste0("HAPLOTYPES/sites_", set, ".txt"))
pos <- as.numeric(data[1, -1])
strain <- sapply(strsplit(data[-1, 1], split = ">"), function(x) x[2])
data <- data[-1, -1]
data <- as.matrix(data)
colnames(data) <- pos
rownames(data) <- strain
data2 <- melt(data)
colnames(data2) <- c("strain", "pos", "nuc")

## add the population (fineSTRUCTURE)
data2 <- addPop(data2)


## color the sites with three colors
## major allele low load (blue), major high (red), and other (grey)
data2$allele <- "other"
ix <- 1
for(i in unique(data2$pos)) {
    tmp <- data2[data2$pos == i & data2$pop %in% c("hspSiberia_low", "hspIndigenousNAmerica_low"),]
    tmp <- table(tmp$nuc)
    tmp <- attributes(tmp[tmp == max(tmp)])$names
    data2$allele[data2$pos == i & data2$nuc == tmp] <- "low"
    tmp <- data2[data2$pos == i & data2$pop %in% c("hspSiberia_high", "hspIndigenousNAmerica_high"),]
    tmp <- table(tmp$nuc)
    tmp <- attributes(tmp[tmp == max(tmp)])$names
    data2$allele[data2$pos == i & data2$nuc == tmp] <- "high"
    print(ix)
    ix = ix+1
}
data2 <- data2[order(data2$allele, data2$pop, data2$strain, data2$pos),]

## file with the main dataset ie,
## the allele type (low, high, other) for all strain, for all signifcant sites
## with the population, the position and the nucleotide
save(data2, file = "HAPLOTYPES/data_pval10FST0.9.RData")

## count the number of alleles that are the same as the high pop
subdata <- data.frame(strain = unique(data2$strain), nb = NA, pop = NA)
for(i in 1:nrow(subdata)) {
    print(i)
    subdata$pop[i] <- as.character(unique(data2$pop[data2$strain == subdata$strain[i]]))
    subdata$nb[i] <- sum(data2$strain == subdata$strain[i] & data2$allele == "high")
}

## file with the count of "high" allele for each strain
save(subdata, file = "HAPLOTYPES/countHighAllele_pval10FST0.9.RData")

hist(subdata$nb, breaks = 100)
hist(subdata$nb[subdata$nb > 500], breaks = 100)

#########################################################

## VISUALIZE THE HAPLOTYPES FOR ALL THE STRAINS
#########################################################
setwd("/Users/elise/Desktop/HpGlobal_noHpGP/HAPLOTYPES/")

set <- "pval10FST0.9"
load("data_pval10FST0.9.RData")
load("countHighAllele_pval10FST0.9.RData")

gff <- read.delim("~/Desktop/HpGlobal/DATA/26695.gff", header=FALSE, comment.char="#")
gff <- gff[gff$V3 == "CDS",]

colAllele <- c("low" = "blue", "high" = "red", "other" = "grey", "limit" = "black")
colPop <- c("hpAfrica2" = "#000000", "hspAfrica1SAfrica" = "#D9D9D9", "hspAfrica1WAfrica" = "#7F7F7F", "hspCNEAfrica" = "#A3D399", "hpNEAfrica" = "#008041", "hspSWEuropeLatinAmerica" = "#D5A6BD", "hspSWEurope" = "#E6A5D6", "hspSEurope" = "#FF8BD3", "hspEuropeMiddleEast" = "#A64D79", "hspNEurope" = "#D62FAE", "hpEurope" = "#741B47", "hspUral" = "#FAB778", "hspSiberia_low" = "#EB3D4F", "hpNorthAsia" = "#E76710", "hspEAsia" = "#920808", "hpAsia2" = "#FFD966", "hpSahul" = "#D5B7ED", "hspIndigenousNAmerica_low" = "#7B4BA0", "hspIndigenousSAmerica_low" = "#AA8AC3", "Hacinonychis" = "turquoise1", "Monkey" = "green", "hspIndigenousSAmerica_high" = "#ABCDEF", "hspIndigenousNAmerica_high" = "#3E5170", "hspSiberia_high" = "#C7EAFF", "hspSiberia" = "#EB3D4F", "hspIndigenousNAmerica" = "#7B4BA0", "hspIndigenousSAmerica" = "#AA8AC3")

databis <- data2

databis <- popData(databis)

subdata$level <- "low"
subdata$level[subdata$nb > 1900] <- "high"

databis <- addCDS(databis, gff)

p1 <- ggplot(databis) +
    geom_tile(aes(x = as.character(pos), y = strain, fill = allele)) +
    scale_fill_manual(values = colAllele) +
    theme_set(theme_bw()) +
    theme(panel.grid = element_line(colour = NA),axis.title = element_blank(),axis.text = element_blank(),axis.line = element_line())

p2 <- ggplot(databis[databis$pos == min(databis$pos),]) +
    geom_tile(aes(x = 1, y = strain, fill = pop)) +
    scale_fill_manual(values = colPop) +
    theme_set(theme_bw()) +
    theme(panel.grid = element_line(colour = NA),axis.title = element_blank(),axis.text = element_blank(),axis.line = element_line(), legend.title = element_text(size = 6), legend.text = element_text(size = 6)) + 
    guides(fill = guide_legend(override.aes = list(size = 0.5)))

leg1 <- get_legend(p1 + theme(legend.position = "bottom"))
leg2 <- get_legend(p2 + theme(legend.position = "bottom"))

ggarrange(p2 + theme(legend.position = "none"), p1+ theme(legend.position = "none"), NULL, leg1, NULL, leg2, ncol = 2, nrow = 3, widths = c(0.9,10), heights = c(10, 1, 1.1))




#########################################################

## BLOCKS OF HIGH ALLELES - CALCULATION AND VISUALIZATIONS
#########################################################

setwd("~/Desktop/HpGlobal_noHpGP/HAPLOTYPES")

load("data_pval10FST0.9.RData")
databis <- popData(data2)

gff <- read.delim("~/Desktop/HpGlobal_noHpGP/DATA/26695.gff", header=FALSE, comment.char="#")
gff <- gff[gff$V3 == "CDS",]

## calculate the number and size of blocks
## to be considered in the same block, the sites must from the same gene
## use the gff file
## if the two sites are not from the same gene, they are from different block, even if they have the same high/low allele
databis <- databis[databis$allele != "limit",]
POS <- sort(unique(databis$pos))
databis$block <- 0
for(i in 2:length(POS)) {
    print(i)
    pos0 <- POS[i-1]
    pos1 <- POS[i]
    which(gff$V4 <= pos0)
    lim0 <- c(max(gff$V4[gff$V4 <= pos0]), min(gff$V5[gff$V5 >= pos0]))
    lim1 <- c(max(gff$V4[gff$V4 <= pos1]), min(gff$V5[gff$V5 >= pos1]))
    strain0 <- databis[databis$pos == pos0,]
    strain1 <- databis[databis$pos == pos1,]
    if(lim0[1] == lim1[1]) { ## ie the beginning of their gene is the same <=> same gene
        ix <- which(strain1$allele == strain0$allele)
        databis$block[databis$pos == pos1][ix] <- strain0$block[ix]
        if(length(ix) < nrow(strain1)) {
            databis$block[databis$pos == pos1][-ix] <- strain0$block[-ix] + 1
        }
    } else {
        databis$block[databis$pos == pos1] <- strain0$block + 1
    }
}
## file with the data pertaining to the blocks
## for all type of alleles
save(databis, file = "dataBlock_pval10FST0.9.RData")

## extract the blocks of high alleles
highBlock <- NULL
j <- 1
for(i in unique(databis$strain[databis$allele == "high"])) {
    print(j)
    strain <- databis[databis$allele == "high" & databis$strain == i,]
    nbBlock <- length(unique(strain$block))
    tmp <- unlist(lapply(unique(strain$block), function(x) max(strain$pos[strain$block == x]) - min(strain$pos[strain$block == x]) + 1))
    highBlock <- rbind(highBlock, data.frame(strain = i, nbBlock = nbBlock, meanLength = mean(tmp), maxLength = max(tmp), minLength = min(tmp)))
    j <- j+1
}
## file with information about the blocks of high alleles
## number, average, max and min size
save(highBlock, file = "highBlock_pval10FST0.9.RData")


## VISUALIZATIONS
load("dataBlock_pval10FST0.9.RData")
load("highBlock_pval10FST0.9.RData")

colPop <- c("hpAfrica2" = "#000000", "hspAfrica1SAfrica" = "#D9D9D9", "hspAfrica1WAfrica" = "#7F7F7F", "hspCNEAfrica" = "#A3D399", "hpNEAfrica" = "#008041", "hspSWEuropeLatinAmerica" = "#D5A6BD", "hspSWEurope" = "#E6A5D6", "hspSEurope" = "#FF8BD3", "hspEuropeMiddleEast" = "#A64D79", "hspNEurope" = "#D62FAE", "hpEurope" = "#741B47", "hspUral" = "#FAB778", "hspSiberia_low" = "#EB3D4F", "hpNorthAsia" = "#E76710", "hspEAsia" = "#920808", "hpAsia2" = "#FFD966", "hpSahul" = "#D5B7ED", "hspIndigenousNAmerica_low" = "#7B4BA0", "hspIndigenousSAmerica_low" = "#AA8AC3", "Hacinonychis" = "#6996CA", "Monkey" = "#16549D", "hspIndigenousSAmerica_high" = "#AA8AC3", "hspIndigenousNAmerica_high" = "#7B4BA0", "hspSiberia_high" = "#EB3D4F")

## look at the "correlation" between the number of "high" blocks and the number of "high" alleles
load("countHighAllele_pval10FST0.9.RData")
subdata <- addPop(subdata)
subdata$strain <- as.character(subdata$strain)
subdata$pop <- as.character(subdata$pop)
subdata <- subdata[order(subdata$strain),]
highBlock <- highBlock[order(highBlock$strain),]
subdata$pop[subdata$strain %in% c("HEL_AA0408AA_AS", "HEL_BA9262AA_AS")] <- "Hacinonychis"
subdata$pop[subdata$strain %in% c("HEL_BA9707AA_AS", "HEL_BA9303AA_AS")] <- "Monkey"

#sum(subdata$strain != highBlock$strain)
subdata$nbBlock <- highBlock$nbBlock
subdata$meanLength <- highBlock$meanLength

ggplot(subdata) +
    geom_histogram(aes((nb - nbBlock)/nb), bins = 100) +
    facet_wrap(~ pop, scales = "free_y")

subdata$type <- "A"
subdata$type[subdata$nb > 500] <- "B"
fit <- lm(nb ~ nbBlock, subdata[subdata$nb < 1000,])
ggplot(subdata) +
    geom_point(aes(x = nbBlock, y = nb, col = pop)) +
    scale_colour_manual(values = colPop) +
    labs(y = "number of 'high' alleles", x = "number of 'high' blocks") +
    theme(legend.position = "bottom") +
    geom_abline(slope = fit$coefficients[2], intercept = fit$coefficients[1]) +
    facet_wrap(~ type, scales = "free_y")



#########################################################

## NUMBER OF HIGH ALLELES VS GENETIC DISTANCE
#########################################################

## plot nbr of "high" alleles vs genetic distance (nb diff / nb sites total) to the low load strains
## average over the distance for all low load strains

setwd("/Users/elise/Desktop/HpGlobal_noHpGP")

library(ggplot2)

colPop <- c("hpAfrica2" = "#000000", "hspAfrica1SAfrica" = "#D9D9D9", "hspAfrica1WAfrica" = "#7F7F7F", "hspCNEAfrica" = "#A3D399", "hpNEAfrica" = "#008041", "hspSWEuropeLatinAmerica" = "#D5A6BD", "hspSWEurope" = "#E6A5D6", "hspSEurope" = "#FF8BD3", "hspEuropeMiddleEast" = "#A64D79", "hspNEurope" = "#D62FAE", "hpEurope" = "#741B47", "hspUral" = "#FAB778", "hspSiberia_high" = "#EB3D4F", "hpNorthAsia" = "#E76710", "hspEAsia" = "#920808", "hpAsia2" = "#FFD966", "hpSahul" = "#D5B7ED", "hspIndigenousSAmerica_high" = "#AA8AC3", "hspIndigenousNAmerica_high" = "#7B4BA0", "Hacinonychis" = "#6996CA", "Monkey" = "#16549D")

## dataset that contains the count of "high" allele for all strains
## name of the variable: subdata
load("HAPLOTYPES/countHighAllele_pval10FST0.9.RData")
subdata <- addPop(data2 = subdata)

## read the genetic distance to the low load strains
## nbdiff_lowLoadStrains.txt is generated using the python script s7ter_nbDiffSeq.py (in bash, for loop to parallelize)
# system("for i in {0..6874..700}; do python3 nbDiffSeq.py $i &; done")
data0 <- read.table("HAPLOTYPES/nbdiff_lowLoadStrains.txt", header = TRUE, sep = "\t")

## I am interested in the proportion of sites that are different, without gaps
data0$dist <- data0$nbDiffNoGap / (data0$nbSite - data0$nbGaps)
data0$dist2  <- data0$nbDiff / data0$nbSite

data <- aggregate(data0[c("dist", "dist2")], by = list(strain = data0$strain), mean)

subdata <- subdata[subdata$strain %in% data$strain,]
subdata$strain <- as.character(subdata$strain)
subdata <- subdata[order(subdata$strain),]
data <- data[data$strain %in% subdata$strain,]
data <- data[order(data$strain),]

data$nbAll <- subdata$nb
data$pop <- as.character(subdata$pop)

data$pop[data$strain %in% c("HEL_AA0408AA_AS", "HEL_BA9262AA_AS")] <- "Hacinonychis"
data$pop[data$strain %in% c("HEL_BA9707AA_AS", "HEL_BA9303AA_AS")] <- "Monkey"


ggplot(data) +
    geom_point(aes(x = dist, y = nbAll, col = pop)) +
    scale_colour_manual(values = colPop) +
    labs(x = "nb differences with low load strains / nb sites [core genes]", y = "number of high load alleles [high FST & signif GWAS sites]") +
    #facet_wrap(~ type, scale = "free_y") +
    theme(legend.position = "bottom")
ggsave("PLOT/HAPLOTYPES/nbHighAlleleVSdistLow.jpeg", height = 9, width = 10)

## some strains are outliers:
## their number of "high" alleles is higher than expected
## compared to the other strains
## for now, exclude the whole hpSahul, as the whole population is globally a bit higher
ix <- which(data$pop != "hpAfrica2" & data$nbAll < 1000 & ((data$nbAll > 86 & data$dist < 0.042) | (data$nbAll >= 150 & data$dist > 0.042 & data$dist < 0.051) | (data$nbAll >= 190 & data$dist > 0.051 & data$dist < 0.054)))
data[ix,]

## outlier strains
## high number of high alleles compared to the rest of their population
# STRAINS <-  c("26695", "HEL_AA7275AA_AS", "HEL_BA9014AA_AS", "HEL_BA9032AA_AS", "HEL_BA9219AA_AS", "HEL_BA9339AA_AS", "HEL_BA9544AA_AS", "HEL_BA9696AA_AS", "HEL_BA9792AA_AS", "HEL_BA9830AA_AS", "HEL_BA9838AA_AS", "HEL_CA2274AA_AS", "HEL_CA2285AA_AS", "HpGP-CAN-012", "HpGP-CAN-014", "HpGP-CAN-017", "HpGP-CHI-006", "HpGP-CHI-021", "HpGP-CHI-022", "HpGP-CHI-027", "HpGP-CHI-120", "HpGP-CHI-124", "HpGP-CHI-127", "Keto60", "Khanty16", "Khanty27", "Khanty40")


#########################################################