## script for the plots

library(readr)
library(tidyverse)
library(ggplot2)
library(reshape2)
library(cowplot)
library(gridExtra)
library(ggpubr)
library(ggtree)
library(ape)
library(viridis)
library(ggsci)

setwd("/home/elise/Desktop/hp/reviews")

colPop <- c("hpAfrica2" = "red", "hpAfrica1" = "#7F7F7F", "hspAfrica1SAfrica" = "#D9D9D9", "hspAfrica1WAfrica" = "#7F7F7F", "hspCNEAfrica" = "#A3D399", "hpNEAfrica" = "#008041", "hspSWEuropeLatinAmerican" = "#D5A6BD", "hspSWEurope" = "#E6A5D6", "hspSEurope" = "#FF8BD3", "hspEuropeMiddleEast" = "#A64D79", "hspNEurope" = "#D62FAE", "hpEurope" = "#D62FAE" #"#741B47"
            , "hspUral" = "#FAB778", "hspSiberia_low" = "#EB3D4F", "hpNorthAsia" = "#E76710", "hspEAsia" = "#920808", "hpAsia2" = "#FFD966", "hpSahul" = "#D5B7ED", "hspIndigenousNAmerica_low" = "#7B4BA0", "hspIndigenousSAmerica_low" = "#AA8AC3", "Hacinonychis" = "limegreen", "Monkey" = "limegreen", "Primate" = "limegreen", "hspIndigenousSAmerica_high" = "#ABCDEF", "hspIndigenousNAmerica_high" = "#3E5170", "hspSiberia_high" = "#C7EAFF", "hspSiberia" = "purple", "hspIndigenousNAmerica" = "blue", "hspIndigenousSAmerica" = "deepskyblue", "Hcetorum" = "deeppink", "other" = "grey", "NA" = "black")

## use the viridis package to get a colorblind-friendly color palette
colVir <- viridis(14)
colVir = c("#440154FF", "#481D6FFF", "#453581FF", "#3D4D8AFF", "#34618DFF", "#2B748EFF", "#24878EFF", "#1F998AFF", "#25AC82FF", "#40BC72FF", "#67CC5CFF", "#97D83FFF", "#CBE11EFF", "#FDE725FF")
colVir <- turbo(14)
# [1] "#30123BFF" "#424AB3FF" "#467EF4FF" "#31AFF5FF" "#18DAC7FF" "#38F491FF"
# [7] "#83FF52FF" "#BDF534FF" "#E9D539FF" "#FEAA33FF" "#F8721CFF" "#E03F08FF"
#[13] "#B61C02FF" "#7A0403FF"

colPop <- c("Hacinonychis" = colVir[1], "hpAfrica1" = colVir[2], "hpAfrica2" = colVir[3], "hpAsia2" = colVir[4], "hpEurope" = colVir[5], "hpNEAfrica" = colVir[6], "hpNorthAsia" = colVir[7], "hpSahul" = colVir[8], "hspEAsia" = colVir[9], "hspIndigenousNAmerica" = colVir[10], "hspIndigenousSAmerica" = colVir[11], "hspSiberia" = colVir[12], "hspUral" = colVir[13], "Hcetorum" = colVir[14])

strainsToRemove <- c("HEL_BA8803AA_AS", "HEL_BA8804AA_AS", "HEL_BA8805AA_AS", "HEL_BA8806AA_AS", "HEL_BA8807AA_AS")


## use all the populations
getAllPop <- function(data, pop) {
    data$pop <- "NA"
    data$divergence <- "other"
    for(i in unique(pop$pop)) {
        strains <- pop$strain[pop$pop == i]
        data$pop[data$strain %in% strains] <- i
    }
    data$pop[data$strain %in% c("HEL_AA0408AA_AS", "HEL_BA9262AA_AS")] <- "Hacinonychis"
    data$pop[data$strain %in% c("HEL_BA9303AA_AS", "HEL_BA9707AA_AS")] <- "hpSahul"
    data$divergence[data$strain %in% pop$strain[pop$divergence == "high"]] <- "high"
    data$divergence[data$strain %in% pop$strain[pop$divergence == "low"]] <- "low"
    
    data$host <- "humans"
    data$host[data$strain %in% c("HEL_AA0408AA_AS", "HEL_BA9262AA_AS")] <- "Hacinonychis"
    data$host[data$strain %in% c("HEL_AA4372AA_AS", "HEL_BA5627AA_AS", "HEL_BA8243AA_AS", "HEL_BA9186AA_AS", "HEL_BA9638AA_AS", "HEL_BA9303AA_AS", "HEL_BA9707AA_AS")] <- "Primate"
    
    return(data)
}


##### declaration pop / load levels #####

x1 <- data.frame(read_csv("DATA/Enterobase_Hp_070322.csv")) ## information from the enterobase database
x2 <- read.table("DATA/Samples_AvgAncestry_PerSample_PopAndSubpop.csv", sep = ",", header = TRUE) ## population assignment based on the chromosome painting results
x2.1 <- read.table("DATA/Donors.txt") ## donor strains for the chromosome painting 
x3 <- read.table("ANALYSIS_OTHER/fineSTRUCTURE_Roberto/run1/NAM_SIB_fsPops.txt", header = TRUE) ## population assignment based on the fineSTRUCTURE results
x4 <- read.table("METADATA/strains_SiberiaIndAm.csv", sep = ",", header = TRUE) ## load levels for the Siberia and IndN/SAm populations
x5 <- read.table("METADATA/strains_pop.csv", sep = ",", header = TRUE) ## "general" population assignment

## concatenate all the information in one dataset
## could add the length of the genome to compare the length between the high and low strains
pop <- data.frame(strain = c(x2$Sample, x2.1$V1), enterobase_name = NA, divergence = NA, pop = NA, pop_G = NA, pop_CHRP = NA, subpop_CHRP = NA, pop_FS = NA, city = NA, country = NA, continent = NA, latitude = NA, longitude = NA, host_ethnicity = NA, origin = NA)

## add the enterobase name
## and geographical / ethnic information
for(i in 1:nrow(pop)) {
    #print(i)
    if(pop$strain[i] %in% x1$Assembly.barcode...1) {
        pop$continent[i] <- x1$Continent...4[x1$Assembly.barcode...1 == pop$strain[i]]
        pop$country[i] <- x1$Country...2[x1$Assembly.barcode...1 == pop$strain[i]]
        pop$city[i] <- x1$City...2[x1$Assembly.barcode...1 == pop$strain[i]]
        pop$latitude[i] <- x1$Latitude[x1$Assembly.barcode...1 == pop$strain[i]]
        pop$longitude[i] <- x1$Longitude[x1$Assembly.barcode...1 == pop$strain[i]]
        pop$host_ethnicity[i] <- x1$Host.Ethnicity[x1$Assembly.barcode...1 == pop$strain[i]]
        pop$enterobase_name[i] <- x1$Name[x1$Assembly.barcode...1 == pop$strain[i]]
    } 
    else if(pop$strain[i] %in% x1$Name) {
        pop$continent[i] <- x1$Continent...4[x1$Name == pop$strain[i]]
        pop$country[i] <- x1$Country...2[x1$Name == pop$strain[i]]
        pop$enterobase_name[i] <- x1$Assembly.barcode...1[x1$Name == pop$strain[i]]
        pop$city[i] <- x1$City...2[x1$Name == pop$strain[i]]
        pop$latitude[i] <- x1$Latitude[x1$Name == pop$strain[i]]
        pop$longitude[i] <- x1$Longitude[x1$Name == pop$strain[i]]
        pop$host_ethnicity[i] <- x1$Host.Ethnicity[x1$Name == pop$strain[i]]
    } 
    else if(pop$strain[i] %in% x2$Sample) {
        pop$continent[i] <- x2$Continent[x2$Sample == pop$strain[i]]
        pop$country[i] <- x2$Country[x2$Sample == pop$strain[i]]
    }
}

## fill the chromosome painting population
for(i in unique(x2$MergedPop)) {
    strains <- x2$Sample[x2$MergedPop == i]
    pop$pop_CHRP[pop$strain %in% strains | pop$enterobase_name %in% strains] <- i
}
## fill the chromosome painting subpopulation
for(i in unique(x2$Subpopulation)) {
    strains <- x2$Sample[x2$Subpopulation == i]
    pop$subpop_CHRP[pop$strain %in% strains | pop$enterobase_name %in% strains] <- i
}
## fill the fineSTRUCTURE population
for(i in unique(x3$fsPop)) {
    strains <- x3$SampleID[x3$fsPop == i]
    pop$pop_FS[pop$strain %in% strains | pop$enterobase_name %in% strains] <- i
}
for(i in unique(x5$pop)) {
    strains <- x5$Sample[x5$pop == i]
    pop$pop_G[pop$strain %in% strains | pop$enterobase_name %in% strains] <- i
}

## the column pop refers to the population assignment that will be used throughout the analysis
## use the fineSTRUCTURE for the pop assignment
pop$pop <- pop$pop_FS
## for the strains that were not used in the fineSTRUCTURE analysis (6278 / 6874 strains), use the chromosome painting pop assignment
pop$pop[is.na(pop$pop)] <- pop$pop_CHRP[is.na(pop$pop)]
## finally, for the strains that have none (14 strains), use the "general" population assignment
pop$pop[is.na(pop$pop)] <- pop$pop_G[is.na(pop$pop)]

## fill in the divergence levels
pop$divergence = "other"
pop$divergence[pop$strain %in% x4$strain[x4$level == "high"]] <- "high"
pop$divergence[pop$strain %in% x4$strain[x4$level == "low"]] <- "low"

## manual curation for the Monkey and H. acinonychis strains
pop$origin <- "human"
pop$origin[pop$strain %in% c("HEL_BA9186AA_AS", "HEL_BA9638AA_AS")] <- "monkey-Asia2"
pop$origin[pop$strain %in% c("HEL_BA9303AA_AS", "HEL_BA9707AA_AS")] <- "monkey-Sahul"
pop$pop[pop$origin == "monkey-Sahul"] <- "Primate"
pop$origin[pop$strain %in% c("HEL_AA0408AA_AS", "HEL_BA9262AA_AS")] <- "Hacinonychis"
pop$pop[pop$origin == "Hacinonychis"] <- "Hacinonychis"
pop$divergence[pop$origin %in% c("monkey-Sahul", "Hacinonychis")] <- "high"

## HEL_BA9721AA_AS and HpGP-ZAF-008 are assigned to hpAfrica2 based on their subpop_CHRP and to hpEurope based on pop_CHRP, due to their quite "tight" ancestry proportion
## their ancestry profile is closer to hpEurope
pop$pop[pop$strain %in% c("HpGP-ZAF-008", "HEL_BA9721AA_AS")] <- "hpEurope"

## remove the subpopulations that are not of interest in the paper
## and replace them by a lower order population
pop$pop[pop$pop %in% c("hspAfrica1SAfrica", "hspAfrica1WAfrica")] <- "hpAfrica1"
pop$pop[pop$pop %in% c("hspNEurope", "hspEuropeMiddleEast", "hspSEurope", "hspSWEurope", "hspSWEuropeLatinAmerican")] <- "hpEurope"
pop$pop[pop$pop %in% c("hspCNEAfrica")] <- "hpNEAfrica"

##########


#### TO REDO WITH THE WHOLE GENOME
###### 2. tree for all strains ######

tree <- read.nexus("/home/elise/Desktop/hp/reviews/rooted_nexus/whole/coreCDS.nex")
tree <- drop.tip(tree, strainsToRemove) ## remove some of the tips (do not want these strains)
div <- data.frame(strain = tree$tip.label, pop = NA, divergence = "other")
div <- getAllPop(div, pop)

strainsNA <- div$strain[div$pop == "NA"]
tree <- drop.tip(tree, strainsNA)
div <- div[div$pop != "NA",]

p <- ggtree(tree, size = 0.05, layout = "circular") +
    geom_treescale(fontsize = 3) 
p2 <- p %<+% div +
    aes(color = pop) +
    geom_tippoint(aes(color = pop, shape = interaction(divergence, host)), size = 1) +
    scale_color_manual(values = colPop[names(colPop) %in% unique(div$pop)]) +
    scale_shape_manual(values = c(1, 19, 4, NA, 15, 0), labels = c("Hardy (H. acinonychis)", "Hardy", "Ubiquitous (hspSiberia, hspIndigenousNAmerica, hspIndigenousSAmerica)", "Ubiquitous", "Hardy (Primate)", "Ubiquitous (Primate)")) +
    theme(plot.margin = unit(c(-50,-100,-10,-200), "mm")) +
    hexpand(.05, direction = -1) +
    labs(color = "population", fill = "population", shape = "ecotype") +
    theme(text = element_text(size = 18), axis.text = element_text(size = 14), legend.position = c(1.1, 0.35)) +
    guides(colour = guide_legend(override.aes = list(size=6)))
ggsave(p2, filename = paste0("DRAFT/figures/TREE2.pdf"), device = "pdf", width = 16, height = 9, units = "in", dpi = 500)

############

#### OK
###### 3. global PCA ######

pca <- data.frame(read_table("./PCA/onAllGenome/hpglobal_LD_PCA.eigenvec", col_names = FALSE))
eigenval <- scan("./PCA/onAllGenome/hpglobal_LD_PCA.eigenval")

pca <- pca[,-1]
names(pca)[1] <- "strain"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))
pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)
pca <- pca[-1,]
pca <- pca[!(pca$strain %in% strainsToRemove),]

pca <- getAllPop(pca, pop)
pca <- pca[pca$pop != "NA",]

ggplot(pca, aes(PC1, PC2, col = pop, shape = divergence, alpha = divergence)) + 
    geom_point() +
    geom_point(data = pca[pca$host == "Primate",], aes(color = pop, shape = "za"), size = 4) +
    geom_point(data = pca[pca$host == "Hacinonychis",], aes(color = pop, shape = "zz"), size = 4) +
    labs(x = paste0("PC1 (", round(pve$pve[1],1), "%)"), y = paste0("PC2 (", round(pve$pve[2],1), "%)"), col = "population", shape = "ecotype", alpha = "ecotype") +
    scale_color_manual(values = colPop[names(colPop) %in% unique(pca$pop)], guide = guide_legend(override.aes = list(size = rep(4,13), linetype = rep(0,13), shape = rep(19,13)))) +
    scale_shape_manual(values = c(19, 3, 3, 0, 1), labels = c("Hardy", "Ubiquitous", "Ubiquitous", "Primates", "H.acinonychis"), guide = guide_legend(override.aes = list(size = rep(4,5)))) +
    scale_alpha_manual(values = c(1, 1, 0.5), labels = c("Hardy", "Ubiquitous", "Ubiquitous")) +
    theme_set(theme_bw()) +
    theme(text = element_text(size = 16), axis.text = element_text(size = 14)) +
    guides(alpha = "none")
ggsave(filename = paste0("DRAFT/figures/set_allPop/PCA.pdf"), device = "pdf", width = 16, height = 9, units = "in", dpi = 500)

############

#### OK
###### 4. GWAS + FST ######

gwas <- read.table("GWAS/onAllGenome/bugwas_biallelic_lmmout_allSNPs.txt", header = TRUE)
gwas <- data.frame(pos = gwas$ps, pval = gwas$negLog10)

load("FST/onAllGenome/FST_highLow.RData") ## slide.data
x <- slide.data@nuc.F_ST.pairwise

tot <- data.frame(pos = gwas$pos, FST = x[gwas$pos])
tot <- tot[order(tot$pos),]
tot$pval <- gwas$pval[order(gwas$pos)]
tot$col <- "black"
tot$col[tot$FST > 0.9] <- "red"
tot$col[tot$FST < 0.5] <- "grey"

## for the limits of the CDS
gff <- read.delim("DATA/26695.gff", header=FALSE, comment.char="#")
gff <- gff[gff$V3 == "CDS",]

ggplot(tot, aes(x = pos, y = pval, col = col)) +
    geom_point(size = 1) +
    labs(x = "position (bp)", y = "-log10(p-val)", color = "") +
    geom_hline(yintercept = 10, col = "green") +
    theme_set(theme_bw()) +
    scale_color_manual(breaks = c("black", "red", "grey"), values = c("black", "red", "grey"), labels = c("0.5 < FST < 0.9", "FST > 0.9", "FST < 0.5"))
ggsave("DRAFT/figures/GWAS_FST.pdf", device = "pdf", width = 16, height = 9, units = "in", dpi = 500, scale = 0.75)

## zoom
CDS <- gff[which((gff$V4 > 1500000 & gff$V5 < 1550000) | (gff$V4 < 1500000 & gff$V5 > 1500000 & gff$V5 < 1550000)| (gff$V5 > 1550000 & gff$V4 > 1500000 & gff$V4 < 1550000)),]
CDSsig <- data.frame(x = c(1516143, 1526523, 1534664), y = rep(-4,3), label = c("HP1445:1446", "HP1456:1457", "Hp1463:1466"))
ggplot(tot[tot$pos > 1500000 & tot$pos < 1550000,]) +
    geom_point(aes(x = pos, y = pval, color = col), size = 1) +
    labs(x = "position (bp)", y = "-log10(p-val)", color = "") +
    theme_set(theme_bw()) +
    geom_hline(yintercept = 10, color = "green") +
    scale_color_manual(breaks = c("black", "red", "grey"), values = c("black", "red", "grey"), labels = c("0.5 < FST < 0.9", "FST > 0.9", "FST < 0.5")) +
    geom_segment(data = CDS[CDS$V7 == "+",], aes(x = V4, y = rep(-2, nrow(CDS[CDS$V7 == "+",])), xend = V5, yend = rep(-2, nrow(CDS[CDS$V7 == "+",]))), size = 1, arrow = arrow(length = unit(0.01, "npc")), linejoin = "mitre") +
    geom_segment(data = CDS[CDS$V7 == "-",], aes(xend = V4, y = rep(-2, nrow(CDS[CDS$V7 == "-",])), x = V5, yend = rep(-2, nrow(CDS[CDS$V7 == "-",]))), size = 1, arrow = arrow(length = unit(0.01, "npc")), linejoin = "mitre") +
    geom_text(data = CDSsig, aes(x = x, y =y , label = label, hjust = 0))
ggsave("DRAFT/figures/GWAS_FST_zoom.pdf", device = "pdf", width = 16, height = 9, units = "in", dpi = 500, scale = 0.75)


aa <- tot[tot$pos > 1500000 & tot$pos < 1550000 & tot$col == "red" & tot$pval > 10,]

gff <- read.delim("DATA/26695.gff", header=FALSE, comment.char="#")
cds <- gff[gff$V3 == "CDS",]
cds <- unique(unlist(sapply(1:nrow(cds), function(x) cds$V4[x]:cds$V5[x])))
noncds <- gff[gff$V3 != "CDS" & !(is.na(gff$V4)),]
noncds <- unique(unlist(sapply(1:nrow(noncds), function(x) noncds$V4[x]:noncds$V5[x])))

totCDS <- tot[tot$pos %in% cds,]
totNonCDS <- tot[tot$pos %in% noncds,]
totIntergenic <- tot[!(tot$pos %in% cds | tot$pos %in% noncds),]

str(totIntergenic[!((totIntergenic$pval < 10 & totIntergenic$FST < 0.5) | (totIntergenic$pval > 10 & totIntergenic$FST > 0.9)),])

############

#### OK
###### 5. dNdS high div genes ######

gene <- "high"

## dN/dS vs dS
data <- read.table(paste0("DIVERGENT_NO_CDS/lowHighGenes/onAllGenome/PAML/dNdS_between_", gene, "/YN00"), header = TRUE)
colnames(data) <- c("strain", colnames(data)[-1])
data <- data[!(data$strain %in% strainsToRemove),]

data <- getAllPop(data, pop)
data <- aggregate(data[c("omega", "dS")], by = list(strain = data$strain, pop = data$pop, divergence = data$divergence, host = data$host), mean)
data <- data[data$pop != "NA",]

p <- ggplot(data) +
    geom_point(aes(x = dS, y = omega, col = pop, shape = divergence, alpha = divergence)) +
    geom_point(data = data[data$host == "Primate",], aes(x = dS, y = omega, color = pop, shape = "za"), size = 4) +
    geom_point(data = data[data$host == "Hacinonychis",], aes(x = dS, y = omega, color = pop, shape = "zz"), size = 4) +
    scale_color_manual(values = colPop[names(colPop) %in% unique(data$pop)], guide = guide_legend(override.aes = list(size = rep(4,12), linetype = rep(0,12), shape = rep(19,12)))) +
    scale_shape_manual(values = c(19, 3, 3, 0, 1), labels = c("Hardy", "Ubiquitous", "Ubiquitous", "Primates", "H.acinonychis"), guide = guide_legend(override.aes = list(size = rep(4,5)))) +
    scale_alpha_manual(values = c(1, 1, 0.5), labels = c("Hardy", "Ubiquitous", "Ubiquitous")) +
    labs(x = "dS", y = "dN/dS", col = "population", subtitle = "outgroup = hpAfrica2", shape = "ecotype", alpha = "ecotype") +
    theme_set(theme_bw()) +
    theme(text = element_text(size = 18), axis.text = element_text(size = 14)) +
    guides(alpha = "none")
ggsave(p, filename = paste0("DRAFT/figures/set_allPop/DNDS_", gene, ".pdf"), device = "pdf", width = 16, height = 9, units = "in", dpi = 500, scale = 0.75)


##

## tree
tree <- read.nexus(paste0("rooted_nexus/treeHigh.nex"))
tree <- drop.tip(tree, strainsToRemove)
div <- data.frame(strain = tree$tip.label, pop = NA, divergence = "other")

div <- getAllPop(div, pop)

strainsNA <- div$strain[div$pop == "NA"]
tree <- drop.tip(tree, strainsNA)
div <- div[div$pop != "NA",]

p <- ggtree(tree, size = 0.05, layout = "circular") +
    geom_treescale(fontsize = 3) 
p2 <- p %<+% div +
    aes(color = pop) +
    geom_tippoint(aes(color = pop, shape = interaction(divergence, host)), size = 1) +
    scale_color_manual(values = colPop[names(colPop) %in% unique(div$pop)]) +
    scale_shape_manual(values = c(1, 19, 4, NA, 15, 0), labels = c("Hardy (H. acinonychis)", "Hardy", "Ubiquitous (hspSiberia, hspIndigenousNAmerica)", "Ubiquitous", "Hardy (Primate)", "Ubiquitous (Primate)")) +
    theme(plot.margin = unit(c(-150,-100,-40,-200), "mm")) +
    hexpand(.05, direction = -1) +
    labs(color = "population", fill = "population", shape = "ecotype", alpha = "ecotype") +
    theme(text = element_text(size = 18), axis.text = element_text(size = 14), legend.position = c(0.97, 0.35)) +
    guides(colour = guide_legend(override.aes = list(size=6)))
ggsave(p2, filename = paste0("DRAFT/figures/TREE_", gene, ".pdf"), device = "pdf", width = 16, height = 9, units = "in", dpi = 500)
##

## dN/dS Hardy vs Ubiquitous
data <- read.table(paste0("DIVERGENT_NO_CDS/lowHighGenes/onAllGenome/PAML/dNdS_", gene, "_HardyUbiquitous/YN00"), header = TRUE)
colnames(data) <- c(colnames(data)[1], "strain", colnames(data)[-c(1:2)])
data <- getAllPop(data, pop)

## separate the divergence calculated against H. acinonychis
## and the divergence calculated against the other Hardy strains
data$outgroup <- "Hardy"
data$outgroup[data$strain1 %in% c("HEL_AA0408AA_AS", "HEL_BA9262AA_AS")] <- "Hacinonychis"

data <- aggregate(data[c("omega", "dS")], by = list(strain = data$strain, pop = data$pop, outgroup = data$outgroup), mean)
data <- data[data$pop != "NA",]

p <- ggplot(data) +
    geom_point(aes(x = dS, y = omega, col = pop), shape = "+", size = 3, alpha = 0.75) +
    scale_color_manual(values = colPop[names(colPop) %in% unique(data$pop)], guide = guide_legend(override.aes = list(size = rep(4,12), linetype = rep(0,12), shape = rep(19,12)))) +
    labs(x = "dS", y = "dN/dS", col = "population") +
    theme_set(theme_bw()) +
    theme(text = element_text(size = 18), axis.text = element_text(size = 14)) +
    facet_wrap(~ outgroup) +
    guides(alpha = "none")
ggsave(p, filename = paste0("DRAFT/figures/set_allPop/DNDS_", gene, "_HardyUbi.pdf"), device = "pdf", width = 16, height = 9, units = "in", dpi = 500, scale = 0.75)
##

############

#### OK
###### 6. dNdS low div genes ######
## same as $5, gene = "low"
gene <- "low"

## dN/dS vs dS
data <- read.table(paste0("DIVERGENT_NO_CDS/lowHighGenes/onAllGenome/PAML/dNdS_between_", gene, "/YN00"), header = TRUE)
colnames(data) <- c("strain", colnames(data)[-1])
data <- data[!(data$strain %in% strainsToRemove),]

data <- getAllPop(data, pop)
data <- aggregate(data[c("omega", "dS")], by = list(strain = data$strain, pop = data$pop, divergence = data$divergence, host = data$host), mean)
data <- data[data$pop != "NA",]

p <- ggplot(data) +
    geom_point(aes(x = dS, y = omega, col = pop, shape = divergence, alpha = divergence)) +
    geom_point(data = data[data$host == "Primate",], aes(x = dS, y = omega, color = pop, shape = "za"), size = 4) +
    geom_point(data = data[data$host == "Hacinonychis",], aes(x = dS, y = omega, color = pop, shape = "zz"), size = 4) +
    scale_color_manual(values = colPop[names(colPop) %in% unique(data$pop)], guide = guide_legend(override.aes = list(size = rep(4,12), linetype = rep(0,12), shape = rep(19,12)))) +
    scale_shape_manual(values = c(19, 3, 3, 0, 1), labels = c("Hardy", "Ubiquitous", "Ubiquitous", "Primates", "H.acinonychis"), guide = guide_legend(override.aes = list(size = rep(4,5)))) +
    scale_alpha_manual(values = c(1, 1, 0.5), labels = c("Hardy", "Ubiquitous", "Ubiquitous")) +
    labs(x = "dS", y = "dN/dS", col = "population", subtitle = "outgroup = hpAfrica2", shape = "ecotype", alpha = "ecotype") +
    theme_set(theme_bw()) +
    theme(text = element_text(size = 18), axis.text = element_text(size = 14)) +
    guides(alpha = "none")
ggsave(p, filename = paste0("DRAFT/figures/set_allPop/DNDS_", gene, ".pdf"), device = "pdf", width = 16, height = 9, units = "in", dpi = 500, scale = 0.75)

##

## tree
tree <- read.nexus(paste0("rooted_nexus/treeLow.nex"))
tree <- drop.tip(tree, strainsToRemove)
div <- data.frame(strain = tree$tip.label, pop = NA, divergence = "other")

div <- getAllPop(div, pop)

strainsNA <- div$strain[div$pop == "NA"]
tree <- drop.tip(tree, strainsNA)
div <- div[div$pop != "NA",]

p <- ggtree(tree, size = 0.05, layout = "circular") +
    geom_treescale(fontsize = 3) 
p2 <- p %<+% div +
    aes(color = pop) +
    geom_tippoint(aes(color = pop, shape = interaction(host,divergence)), size = 1) +
    scale_color_manual(values = colPop[names(colPop) %in% unique(div$pop)]) +
    scale_shape_manual(values = c(1, 19, 15, 4, NA, 0), labels = c("Hardy (H. acinonychis)", "Hardy", "Hardy (Primate)", "Ubiquitous (hspSiberia, hspIndigenousNAmerica)", "Ubiquitous", "Ubiquitous (Primate)")) +
    theme(plot.margin = unit(c(-50,-100,-10,-200), "mm")) +
    hexpand(.05, direction = -1) +
    labs(color = "population", fill = "population", shape = "ecotype", alpha = "ecotype") +
    theme(text = element_text(size = 18), axis.text = element_text(size = 14), legend.position = c(1.1, 0.35)) +
    guides(colour = guide_legend(override.aes = list(size=6)))
ggsave(p2, filename = paste0("DRAFT/figures/TREE_", gene, ".pdf"), device = "pdf", width = 16, height = 9, units = "in", dpi = 500)
##

## dN/dS Hardy vs Ubiquitous
data <- read.table(paste0("DIVERGENT_NO_CDS/lowHighGenes/onAllGenome/PAML/dNdS_", gene, "_HardyUbiquitous/YN00"), header = TRUE)
colnames(data) <- c(colnames(data)[1], "strain", colnames(data)[-c(1:2)])
data <- getAllPop(data, pop)

## separate the divergence calculated against H. acinonychis
## and the divergence calculated against the other Hardy strains
data$outgroup <- "Hardy"
data$outgroup[data$strain1 %in% c("HEL_AA0408AA_AS", "HEL_BA9262AA_AS")] <- "Hacinonychis"

data <- aggregate(data[c("omega", "dS")], by = list(strain = data$strain, pop = data$pop, outgroup = data$outgroup), mean)
data <- data[data$pop != "NA",]

p <- ggplot(data) +
    geom_point(aes(x = dS, y = omega, col = pop), shape = "+", size = 3, alpha = 0.75) +
    scale_color_manual(values = colPop[names(colPop) %in% unique(data$pop)], guide = guide_legend(override.aes = list(size = rep(4,12), linetype = rep(0,12), shape = rep(19,12)))) +
    labs(x = "dS", y = "dN/dS", col = "population") +
    theme_set(theme_bw()) +
    theme(text = element_text(size = 18), axis.text = element_text(size = 14)) +
    facet_wrap(~ outgroup) +
    guides(alpha = "none")
ggsave(p, filename = paste0("DRAFT/figures/set_allPop/DNDS_", gene, "_HardyUbi.pdf"), device = "pdf", width = 16, height = 9, units = "in", dpi = 500, scale = 0.75)
##

############

#### OK
###### 8. haplotypes sig sites ######

set <- "HAPLOTYPES/onAllGenome/pval10FST0.9"
load("HAPLOTYPES/onAllGenome/data_pval10FST0.9.RData")
load("HAPLOTYPES/onAllGenome/countHighAllele_pval10FST0.9.RData")

gff <- read.delim("DATA/26695.gff", header=FALSE, comment.char="#")
gff <- gff[gff$V3 == "CDS",]

colAllele <- c("low" = "blue", "high" = "red", "other" = "grey100", "limit" = "black")

data2$strain <- as.character(data2$strain)

## read the subset of strains that will be used 
## from the file METADATA/strains_pangenome.txt
strainsSubset <- read.table("METADATA/strains_pangenome.txt", header = TRUE)
strainsSubset <- strainsSubset[!(strainsSubset$Sample %in% strainsToRemove),]
databis <- data2[data2$strain %in% strainsSubset$Sample,]

## if want to look at all the strains
#databis <- data2

databis <- getAllPop(databis, pop)

## add the CDS limits
## get the limits of the cds
cdsLim <- c(gff$V4, gff$V5)
cdsLim <- cdsLim[cdsLim >= min(databis$pos) & cdsLim <= max(databis$pos)]
## want to only keep the CDS limits
## that actually corresponds to sites 
cdsLim <- cdsLim[round(cdsLim/1000) %in% unique(round(unique(databis$pos/1000)))]
uniqueSet <- databis[databis$pos == min(databis$pos),]
## add the limits of the CDS to the dataset
## with nuc = NA and allele = limit (color = black)
j <- cdsLim[1]
for(i in cdsLim[-1]) { 
    print(i)
    if(nrow(databis[databis$pos < i & databis$pos > j,]) > 0) databis <- rbind(databis, data.frame(strain = uniqueSet$strain, pos = i, nuc = "NA", pop = uniqueSet$pop, allele = "limit", divergence = uniqueSet$divergence, host = uniqueSet$host))
    j <- i
}


databis$strain <- factor(databis$strain, levels = unique(c(databis$strain[databis$pop == "hspIndigenousNAmerica" & databis$divergence == "high"], databis$strain[databis$pop == "hspIndigenousSAmerica" & databis$divergence == "high"], databis$strain[databis$pop == "hspSiberia" & databis$divergence == "high"], databis$strain[databis$pop == "Hacinonychis"], databis$strain[databis$pop == "hpSahul" & databis$divergence == "high"], databis$strain[databis$pop == "hpAfrica2"], databis$strain[databis$pop == "hspAfrica1SAfrica"], databis$strain[databis$pop == "hspAfrica1WAfrica"], databis$strain[databis$pop == "hpAfrica1"], databis$strain[databis$pop == "hspCNEAfrica"], databis$strain[databis$pop == "hpNEAfrica"], databis$strain[databis$pop == "hspEuropeMiddleEast"], databis$strain[databis$pop == "hspSEurope"], databis$strain[databis$pop == "hspSWEurope"], databis$strain[databis$pop == "hspNEurope"], databis$strain[databis$pop == "hpEurope"], databis$strain[databis$pop == "hspUral"], databis$strain[databis$pop == "hpNorthAsia"], databis$strain[databis$pop == "hspEAsia"], databis$strain[databis$pop == "hpAsia2"], databis$strain[databis$pop == "hpSahul" & databis$divergence == "low"], databis$strain[databis$pop == "hspIndigenousNAmerica" & databis$divergence == "low"], databis$strain[databis$pop == "hspIndigenousSAmerica" & databis$divergence == "low"], databis$strain[databis$pop == "hspSiberia" & databis$divergence == "low"])))

p1 <- ggplot(databis) +
    geom_tile(aes(x = as.character(pos), y = strain, fill = allele)) +
    scale_fill_manual(values = colAllele) +
    theme_set(theme_bw()) +
    theme(panel.grid = element_line(colour = NA),axis.title = element_blank(),axis.text = element_blank(),axis.line = element_line(), legend.title = element_text(size = 12), legend.text = element_text(size = 12))
p2 <- ggplot(databis[databis$pos == min(databis$pos),]) +
    geom_tile(aes(x = 1, y = strain, fill = pop)) +
    scale_fill_manual(values = colPop) +
    labs(fill = "population") +
    theme_set(theme_bw()) +
    theme(panel.grid = element_line(colour = NA),axis.title = element_blank(),axis.text = element_blank(),axis.line = element_line(), legend.title = element_text(size = 12), legend.text = element_text(size = 12)) +
    guides(fill = guide_legend(override.aes = list(size = 0.5)))
leg1 <- get_legend(p1 + theme(legend.position = "bottom"))
leg2 <- get_legend(p2 + theme(legend.position = "bottom"))
ggarrange(p2 + theme(legend.position = "none"), p1+ theme(legend.position = "none"), NULL, leg1, NULL, leg2, ncol = 2, nrow = 3, widths = c(0.9,10), heights = c(10, 0.7, 1.1))
ggsave(paste0("DRAFT/figures/set_allPop/haplotypes_v2.pdf"), dpi = 300, width = 16.2, height = 15)
ggsave(paste0("DRAFT/figures/set_allPop/haplotypes2bis_v2.png"), dpi = 500, width = 20, height = 15)



############

#### OK
###### 9. number high alleles vs dist low strains ######

load("HAPLOTYPES/onAllGenome/countHighAllele_pval10FST0.9.RData") ## subdata
subdata <- getAllPop(subdata, pop)
subdata <- subdata[subdata$pop != "NA",]

data0 <- read.table("HAPLOTYPES/onAllGenome/nbdiff_lowLoadStrains.txt", header = TRUE, sep = "\t")
data0 <- data0[data0$strain != data0$strainLow,] ## remove the pairs for which the strains are the same
data0$dist <- data0$nbDiffNoGap / (data0$nbSite - data0$nbGaps)
data <- aggregate(data0["dist"], by = list(strain = data0$strain), mean)
data <- getAllPop(data, pop)
data <- data[data$pop != "NA",]

subdata <- subdata[subdata$strain %in% data$strain,]
subdata$strain <- as.character(subdata$strain)
subdata <- subdata[order(subdata$strain),]
data <- data[data$strain %in% subdata$strain,]
data <- data[order(data$strain),]
data$nbAll <- subdata$nb

data <- data[!(data$strain %in% strainsToRemove),]

## check the outliers
lim1 <- 0.043
lim2 <- 140
lim3 <- 100
fit <- lm(nbAll ~ dist, data = data[data$nbAll < 1000 & data$dist > lim1,])
ix <- which(data$pop != "hpAfrica2" & data$nbAll < 1000 & 
                ((data$dist < lim1 & data$nbAll > lim2) | 
                     (data$dist >= lim1 & data$nbAll > (fit$coefficients[2]*data$dist+fit$coefficients[1] + lim3))))
outliers <- data[ix,]
popOutliers <- pop[pop$strain %in% outliers$strain,]

plot(data$dist, data$nbAll, pch = "+")
points(outliers$dist, outliers$nbAll, col = "red")
abline(h = lim2, col = "green")
abline(v = lim1, col = "green")
abline(a = fit$coefficients[1] + lim3, b = fit$coefficients[2], col = "green")

write.table(popOutliers, file = "DRAFT/tables/listOutliers_nbAllDist.csv", sep = ",", quote = FALSE, col.names = TRUE, row.names = FALSE, append = FALSE)

data <- getAllPop(data, pop)
data$outliers <- "0"
data$outliers[ix] <- 1
data$divergence[data$divergence == "low"] <- "other"
ggplot(data) +
    geom_point(aes(x = dist, y = nbAll, col = pop, shape = interaction(divergence, outliers)), size = 2) +
    geom_point(data = data[data$host == "Primate",], aes(x = dist, y = nbAll, color = pop, shape = "za"), size = 4) +
    geom_point(data = data[data$host == "Hacinonychis",], aes(x = dist, y = nbAll, color = pop, shape = "zz"), size = 4) +
    scale_colour_manual(values = colPop[names(colPop) %in% unique(data$pop)]) +
    scale_shape_manual(values = c("high.0" = 19, "other.0" = 3, "other.1" = 4, "za" = 0, "zz" = 1), labels = c("Hardy", "Ubiquitous", "\noutliers\n", "Primate", "H. acinonychis"), guide = guide_legend(override.aes = list(size = rep(4,5)))) +
    labs(x = "# differences with Ubiquitous strains / # sites", y = "number of Hardy alleles", col = "populations", shape = "ecotype") +
    theme_set(theme_bw()) +
    theme(text = element_text(size = 18), axis.text = element_text(size = 14))
ggsave("DRAFT/figures/set_allPop/nbHighAlleleVSdistLow.pdf", device = "pdf", width = 16, height = 9, units = "in", dpi = 500, scale = 0.75)



############

#### OK
###### 10. number high alleles vs number high blocks ######


load("HAPLOTYPES/onAllGenome/dataBlock_pval10FST0.9.RData") ## databis
load("HAPLOTYPES/onAllGenome/highBlock_pval10FST0.9.RData") ## highBlock
load("HAPLOTYPES/onAllGenome/countHighAllele_pval10FST0.9.RData") ## subdata

outliers <- read.table("DRAFT/tables/listOutliers_nbAllDist.csv", sep = ",", header = TRUE)

subdata$outliers <- 0
subdata$outliers[subdata$strain %in% outliers$strain] <- 1

subdata <- getAllPop(subdata, pop)

subdata <- subdata[subdata$strain %in% unique(highBlock$strain),]
highBlock <- highBlock[highBlock$strain %in% unique(as.character(subdata$strain)),]

subdata$strain <- as.character(subdata$strain)
subdata <- subdata[order(subdata$strain),]
highBlock <- highBlock[order(highBlock$strain),]
#sum(subdata$strain != highBlock$strain)
subdata$nbBlock <- highBlock$nbBlock
subdata$meanLength <- highBlock$meanLength

subdata <- subdata[!(subdata$strain %in% strainsToRemove),]
subdata <- subdata[subdata$pop != "NA",]


subdata$divergence[subdata$divergence == "other"] <- "low"
subdata$divergence[subdata$outliers == 1] <- "outliers"
ggplot(subdata) +
    geom_point(aes(x = nbBlock, y = nb, col = pop, shape = divergence), size = 2) +
    geom_point(data = subdata[subdata$host == "Primate",], aes(x = nbBlock, y = nb, color = pop, shape = "za"), size = 4) +
    geom_point(data = subdata[subdata$host == "Hacinonychis",], aes(x = nbBlock, y = nb, color = pop, shape = "zz"), size = 4) +
    scale_colour_manual(values = colPop[names(colPop) %in% unique(subdata$pop)]) +
    scale_shape_manual(values = c("high" = 19, "low" = 3, "outliers" = 4, "za" = 0, "zz" = 1), labels = c("Hardy", "Ubiquitous", "\noutliers\n", "Primate", "H. acinonychis"), guide = guide_legend(override.aes = list(size = rep(4,5)))) +
    labs(y = "number of Hardy alleles", x = "number of Hardy blocks", col = "populations", shape = "ecotype") +
    theme_set(theme_bw()) +
    theme(text = element_text(size = 18), axis.text = element_text(size = 14))
ggsave("DRAFT/figures/set_allPop/highBlockAllele.pdf", device = "pdf", width = 16, height = 9, units = "in", dpi = 500, scale = 0.75)


############

#### OK
###### 11. tree VacA ######

## add the VacA alleles manually?

tree <- read.nexus("rooted_nexus/HpGP_VacA.nex")
div <- data.frame(straini = tree$tip.label, pop = NA, divergence = "other")
div$strain <- sapply(strsplit(div$straini, split = ".", fixed = TRUE), function(x) x[1])

div <- getAllPop(div, pop)
div$pop[div$strain == "AFI05407"] <- "Hcetorum"

p <- ggtree(tree, layout = "circular") +
    geom_treescale() 
p2 <- p %<+% div +
    aes(color = pop) +
    geom_tippoint(aes(color = pop, shape = divergence), size = 1) +
    geom_tiplab(geom = "text", align = TRUE) +
    scale_color_manual(values = colPop[names(colPop) %in% unique(div$pop)]) +
    scale_shape_manual(values = c(19, NA, NA), labels = c("Hardy", "Ubiquitous", "Ubiquitous")) +
    theme(plot.margin = unit(c(20,0,20,0), "mm")) +
    hexpand(.05, direction = -1) +
    labs(color = "population", fill = "population", shape = "ecotype", alpha = "ecotype") +
    theme(text = element_text(size = 18), axis.text = element_text(size = 14), legend.position = c(1.35,0.5)) +
    guides(colour = guide_legend(override.aes = list(size=6)))
ggsave(p2, filename = paste0("DRAFT/figures/TREE_VacA.pdf"), device = "pdf", width = 16, height = 9, units = "in", dpi = 500)

############

#### OK
###### 12. enrichment analysis ######
data <- data.frame(read_delim("DIVERGENT_NO_CDS/lowHighGenes/functionalEnrichment_high/david_results/functional_annotation_chart.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE))
data$significant <- " "
data$significant[data$Benjamini < 0.05] <- "*"

data$Term <- factor(data$Term, levels = data$Term[order(data$Count, decreasing = TRUE)])
data$Term2 <- paste0(data$Term, " \n(n = ", data$Count, ")")
data$Term2 <- sapply(strsplit(as.character(data$Term2), split = "~|:|,"), paste, collapse = "\n")
data$Term2 <- factor(data$Term2, levels = data$Term2[order(data$Fold.Enrichment, decreasing = TRUE)])

ggplot(data[data$significant == "*",], aes(x = Term2, y = Fold.Enrichment)) +
    geom_bar(stat = 'identity') +
    #geom_text(aes(label = significant), size = 8) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    labs(x = "term", y = "fold enrichment") +
    theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
    theme(text = element_text(size = 16), axis.text = element_text(size = 14))
ggsave(filename = paste0("DRAFT/figures/enrichment.pdf"), device = "pdf", width = 10, height = 10, units = "in", dpi = 500)
############

#### OK
###### 13. tree UreA ######

gene <- "UreA"
Hcetorum <- c("AFI05732", "AFI05559")

seqType <- "_nt"
tree <- read.nexus(paste0("rooted_nexus/", gene, seqType, ".nex"))
tree <- drop.tip(tree, strainsToRemove)
div <- data.frame(straini = tree$tip.label, pop = NA, divergence = "other")
div$strain <- sapply(strsplit(div$straini, split = ".", fixed = TRUE), function(x) x[1])

div <- getAllPop(div, pop)
div$pop[div$strain %in% Hcetorum] <- "Hcetorum"
p <- ggtree(tree, layout = "circular") +
    geom_treescale() 
p2 <- p %<+% div +
    aes(color = pop) +
    geom_tippoint(aes(color = pop, shape = interaction(divergence, host)), size = 1) +
    #geom_tiplab(geom = "text", align = TRUE) +
    scale_color_manual(values = colPop[names(colPop) %in% unique(div$pop)]) +
    scale_shape_manual(values = c(1, 19, NA, NA, 15, 0), labels = c("Hardy (H. acinonychis)", "Hardy", "Ubiquitous", "Ubiquitous", "Hardy (Primate)", "Ubiquitous (Primate)")) +
    theme(plot.margin = unit(c(-150,-200,-50,-200), "mm")) +
    hexpand(.05, direction = -1) +
    labs(color = "population", fill = "population", shape = "ecotype", alpha = "ecotype") +
    theme(text = element_text(size = 18), axis.text = element_text(size = 14), legend.position = c(0.8,0.3)) +
    guides(colour = guide_legend(override.aes = list(size=6)))
ggsave(p2, filename = paste0("DRAFT/figures/TREE_", gene, seqType, ".pdf"), device = "pdf", width = 16, height = 9, units = "in", dpi = 500)

seqType <- "_nt1"
tree <- read.nexus(paste0("rooted_nexus/", gene, seqType, ".nex"))
tree <- drop.tip(tree, strainsToRemove)
div <- data.frame(straini = tree$tip.label, pop = NA, divergence = "other")
div$strain <- sapply(strsplit(div$straini, split = ".", fixed = TRUE), function(x) x[1])

div <- getAllPop(div, pop)
div$pop[div$strain %in% Hcetorum] <- "Hcetorum"
p <- ggtree(tree, layout = "circular") +
    geom_treescale() 
p2 <- p %<+% div +
    aes(color = pop) +
    geom_tippoint(aes(color = pop, shape = interaction(divergence, host)), size = 1) +
    #geom_tiplab(geom = "text", align = TRUE) +
    scale_color_manual(values = colPop[names(colPop) %in% unique(div$pop)]) +
    scale_shape_manual(values = c(1, 19, NA, NA, 15, 0), labels = c("Hardy (H. acinonychis)", "Hardy", "Ubiquitous", "Ubiquitous", "Hardy (Primate)", "Ubiquitous (Primate)")) +
    theme(plot.margin = unit(c(-120,-270,-120,-270), "mm")) +
    hexpand(.05, direction = -1) +
    labs(color = "population", fill = "population", shape = "ecotype", alpha = "ecotype") +
    theme(text = element_text(size = 18), axis.text = element_text(size = 14), legend.position = c(0.85,0.5)) +
    guides(colour = guide_legend(override.aes = list(size=6)))
ggsave(p2, filename = paste0("DRAFT/figures/TREE_", gene, seqType, ".pdf"), device = "pdf", width = 16, height = 9, units = "in", dpi = 500)

seqType <- "_nt2"
tree <- read.nexus(paste0("rooted_nexus/", gene, seqType, ".nex"))
tree <- drop.tip(tree, strainsToRemove)
div <- data.frame(straini = tree$tip.label, pop = NA, divergence = "other")
div$strain <- sapply(strsplit(div$straini, split = ".", fixed = TRUE), function(x) x[1])

div <- getAllPop(div, pop)
div$pop[div$strain %in% Hcetorum] <- "Hcetorum"
p <- ggtree(tree, layout = "circular") +
    geom_treescale() 
p2 <- p %<+% div +
    aes(color = pop) +
    geom_tippoint(aes(color = pop, shape = interaction(divergence, host)), size = 1) +
    #geom_tiplab(geom = "text", align = TRUE) +
    scale_color_manual(values = colPop[names(colPop) %in% unique(div$pop)]) +
    scale_shape_manual(values = c(1, 19, NA, NA, 15, 0), labels = c("Hardy (H. acinonychis)", "Hardy", "Ubiquitous", "Ubiquitous", "Hardy (Primate)", "Ubiquitous (Primate)")) +
    theme(plot.margin = unit(c(-100,-270,-120,-270), "mm")) +
    hexpand(.05, direction = -1) +
    labs(color = "population", fill = "population", shape = "ecotype", alpha = "ecotype") +
    theme(text = element_text(size = 18), axis.text = element_text(size = 14), legend.position = c(0.85,0.5)) +
    guides(colour = guide_legend(override.aes = list(size=6)))
ggsave(p2, filename = paste0("DRAFT/figures/TREE_", gene, seqType, ".pdf"), device = "pdf", width = 16, height = 9, units = "in", dpi = 500)

############

#### OK
###### 14. tree UreB ######

gene <- "UreB"
Hcetorum <- c("AFI05560", "AFI05733")

seqType <- "_nt"
tree <- read.nexus(paste0("rooted_nexus/", gene, seqType, ".nex"))
tree <- drop.tip(tree, strainsToRemove)
div <- data.frame(straini = tree$tip.label, pop = NA, divergence = "other")
div$strain <- sapply(strsplit(div$straini, split = ".", fixed = TRUE), function(x) x[1])

div <- getAllPop(div, pop)
div$pop[div$strain %in% Hcetorum] <- "Hcetorum"
p <- ggtree(tree, layout = "circular") +
    geom_treescale() 
p2 <- p %<+% div +
    aes(color = pop) +
    geom_tippoint(aes(color = pop, shape = interaction(divergence, host)), size = 1) +
    #geom_tiplab(geom = "text", align = TRUE) +
    scale_color_manual(values = colPop[names(colPop) %in% unique(div$pop)]) +
    scale_shape_manual(values = c(1, 19, NA, NA, 15, 0), labels = c("Hardy (H. acinonychis)", "Hardy", "Ubiquitous", "Ubiquitous", "Hardy (Primate)", "Ubiquitous (Primate)")) +
    theme(plot.margin = unit(c(-150,-200,-50,-400), "mm")) +
    hexpand(.05, direction = -1) +
    labs(color = "population", fill = "population", shape = "ecotype", alpha = "ecotype") +
    theme(text = element_text(size = 18), axis.text = element_text(size = 14), legend.position = c(1,0.45)) +
    guides(colour = guide_legend(override.aes = list(size=6)))
ggsave(p2, filename = paste0("DRAFT/figures/TREE_", gene, seqType, ".pdf"), device = "pdf", width = 16, height = 9, units = "in", dpi = 500)


seqType <- "_nt1"
tree <- read.nexus(paste0("rooted_nexus/", gene, seqType, ".nex"))
tree <- drop.tip(tree, strainsToRemove)
div <- data.frame(straini = tree$tip.label, pop = NA, divergence = "other")
div$strain <- sapply(strsplit(div$straini, split = ".", fixed = TRUE), function(x) x[1])

div <- getAllPop(div, pop)
div$pop[div$strain %in% Hcetorum] <- "Hcetorum"
p <- ggtree(tree, layout = "circular") +
    geom_treescale() 
p2 <- p %<+% div +
    aes(color = pop) +
    geom_tippoint(aes(color = pop, shape = interaction(divergence, host)), size = 1) +
    #geom_tiplab(geom = "text", align = TRUE) +
    scale_color_manual(values = colPop[names(colPop) %in% unique(div$pop)]) +
    scale_shape_manual(values = c(1, 19, NA, NA, 15, 0), labels = c("Hardy (H. acinonychis)", "Hardy", "Ubiquitous", "Ubiquitous", "Hardy (Primate)", "Ubiquitous (Primate)")) +
    theme(plot.margin = unit(c(-150,-270,-50,-270), "mm")) +
    hexpand(.05, direction = -1) +
    labs(color = "population", fill = "population", shape = "ecotype", alpha = "ecotype") +
    theme(text = element_text(size = 18), axis.text = element_text(size = 14), legend.position = c(0.9,0.35)) +
    guides(colour = guide_legend(override.aes = list(size=6)))
ggsave(p2, filename = paste0("DRAFT/figures/TREE_", gene, seqType, ".pdf"), device = "pdf", width = 16, height = 9, units = "in", dpi = 500)


seqType <- "_nt2"
tree <- read.nexus(paste0("rooted_nexus/", gene, seqType, ".nex"))
tree <- drop.tip(tree, strainsToRemove)
div <- data.frame(straini = tree$tip.label, pop = NA, divergence = "other")
div$strain <- sapply(strsplit(div$straini, split = ".", fixed = TRUE), function(x) x[1])

div <- getAllPop(div, pop)
div$pop[div$strain %in% Hcetorum] <- "Hcetorum"
p <- ggtree(tree, layout = "circular") +
    geom_treescale() 
p2 <- p %<+% div +
    aes(color = pop) +
    geom_tippoint(aes(color = pop, shape = interaction(divergence, host)), size = 1) +
    #geom_tiplab(geom = "text", align = TRUE) +
    scale_color_manual(values = colPop[names(colPop) %in% unique(div$pop)]) +
    scale_shape_manual(values = c(1, 19, NA, NA, 15, 0), labels = c("Hardy (H. acinonychis)", "Hardy", "Ubiquitous", "Ubiquitous", "Hardy (Primate)", "Ubiquitous (Primate)")) +
    theme(plot.margin = unit(c(-50,-270,-50,-270), "mm")) +
    hexpand(.05, direction = -1) +
    labs(color = "population", fill = "population", shape = "ecotype", alpha = "ecotype") +
    theme(text = element_text(size = 18), axis.text = element_text(size = 14), legend.position = c(0.8,0.3)) +
    guides(colour = guide_legend(override.aes = list(size=6)))
ggsave(p2, filename = paste0("DRAFT/figures/TREE_", gene, seqType, ".pdf"), device = "pdf", width = 16, height = 9, units = "in", dpi = 500)

############

#### OK
###### 15. map origin ######
library(maps)
library(rworldmap)
library(ggplot2)

data0 <- data.frame(read_csv("DATA/Enterobase_Hp_070322.csv"))
data0 <- data0[,c(1, 2, 3, 4, 8, 15, 28, 29)]
colnames(data0) <- c("Assembly.barcode", "Country", "City", "Continent", "Name", "Host.Ethnicity", "Latitude", "Longitude")
data <- data0
## for the strains with a country but without coordinates, add them based on the country / city given
## need to make a distinction between the strains for which:
## - only the country is given
## - latitude and longitude are given 
## - host ethnicity is also given

divpop <- read.table("METADATA/strains_SiberiaIndAm.csv", sep = ",", header = TRUE)
outliers <- read.table("DRAFT/tables/listOutliers_nbAllDist.csv", sep = ",", header = TRUE)

data <- data[data$Assembly.barcode %in% divpop$strain | data$Name %in% divpop$strain | data$Assembly.barcode %in% outliers$strain | data$Name %in% outliers$strain,]
data$level <- "low"
data$level[data$Assembly.barcode %in% divpop$strain[divpop$level == "high"] | data$Name %in% divpop$strain[divpop$level == "high"]] <- "high"
data$level[data$Assembly.barcode %in% outliers$strain | data$Name %in% outliers$strain] <- "outliers"

ix <- which(!(divpop$strain %in% data$Assembly.barcode) & !(divpop$strain %in% data$Name))
data <- rbind(data, data.frame(Assembly.barcode = divpop$strain[ix], Country = divpop$country[ix], City = NA, Continent = divpop$continent[ix], Name = NA, Host.Ethnicity = divpop$host_ethnicity[ix], Latitude = NA, Longitude = NA, level = divpop$level[ix]))
ix <- which(!(outliers$strain %in% data$Assembly.barcode) & !(outliers$strain %in% data$Name))
data <- rbind(data, data.frame(Assembly.barcode = outliers$strain[ix], Country = outliers$country[ix], City = NA, Continent = outliers$continent[ix], Name = NA, Host.Ethnicity = outliers$host_ethnicity[ix], Latitude = NA, Longitude = NA, level = "outliers"))

## fill the missing latitude / longitude manually
## if the host ethnicity is given, use the coordinates of other strains with same ethnicity and country and known coordinates
## otherwise, if the host ethnicity is given use it to approximately locate the strain origin
## otherwise, use the coordinates of the country (center of the country? same coordinates as other strains from the same country but without known host ethnicity?)
data$Latitude[data$Assembly.barcode %in% c("HEL_AA0365AA_AS", "HEL_BA9433AA_AS")] <- 7.14 ## Venezuela -> center
data$Longitude[data$Assembly.barcode %in% c("HEL_AA0365AA_AS", "HEL_BA9433AA_AS")] <- -66.31
data$Latitude[data$Assembly.barcode %in% c("HEL_AA6482AA_AS")] <- 70 ## Inuit, Canada -> center of Nunavut
data$Longitude[data$Assembly.barcode %in% c("HEL_AA6482AA_AS")] <- -90
data$Latitude[data$Assembly.barcode %in% c("HEL_AA0336AA_AS", "HEL_AA0337AA_AS")] <- 68.22 ## Amerind, Canada -> Aklavik
data$Longitude[data$Assembly.barcode %in% c("HEL_AA0336AA_AS", "HEL_AA0337AA_AS")] <- -135.01
data$Latitude[data$Assembly.barcode %in% c("HEL_AA0369AA_AS", "HEL_AA0375AA_AS")] <- -9.19 ## Peru
data$Longitude[data$Assembly.barcode %in% c("HEL_AA0369AA_AS", "HEL_AA0375AA_AS")] <- -75.01
data$Latitude[data$Assembly.barcode %in% c("HEL_AA0344AA_AS", "HEL_AA0345AA_AS", "HEL_AA0346AA_AS", "HEL_AA0406AA_AS")] <- -10.75 ## Amerind, Peru -> Shima
data$Longitude[data$Assembly.barcode %in% c("HEL_AA0344AA_AS", "HEL_AA0345AA_AS", "HEL_AA0346AA_AS", "HEL_AA0406AA_AS")] <- -73.74
data$Latitude[data$Assembly.barcode %in% c("HEL_AA0357AA_AS", "HEL_AA0358AA_AS")] <- -15.84 ## Amerind, Peru -> Puno
data$Longitude[data$Assembly.barcode %in% c("HEL_AA0357AA_AS", "HEL_AA0358AA_AS")] <- -70.02
data$Latitude[data$Assembly.barcode %in% c("HEL_AA0370AA_AS")] <- -16.40 ## Amerind, Peru -> Monte Carmelo
data$Longitude[data$Assembly.barcode %in% c("HEL_AA0370AA_AS")] <- -71.52
data$Latitude[data$Assembly.barcode %in% c("HEL_AA6759AA_AS")] <- 27.51 ## Tarahumara, Mexico -> Copper Canyon, Chihuahua
data$Longitude[data$Assembly.barcode %in% c("HEL_AA6759AA_AS")] <- -107.71
data$Latitude[data$Assembly.barcode %in% c("HEL_AA6511AA_AS")] <- 23 ## Otomi, Mexico -> Central Plateau (Mexican)
data$Longitude[data$Assembly.barcode %in% c("HEL_AA6511AA_AS")] <- -102
data$Latitude[data$Assembly.barcode %in% c("HEL_BA9042AA_AS", "HEL_BA9750AA_AS")] <- 59.05 ## Koryak, Russia
data$Longitude[data$Assembly.barcode %in% c("HEL_BA9042AA_AS", "HEL_BA9750AA_AS")] <- 159.56
data$Latitude[data$Assembly.barcode %in% c("HEL_BA9047AA_AS")] <- 64.16 ## Evenk (Evenky), Russia
data$Longitude[data$Assembly.barcode %in% c("HEL_BA9047AA_AS")] <- 100.12
data$Latitude[data$Assembly.barcode %in% c("HEL_BA9244AA_AS")] <- 56.06 ## Evenk (Even), Russia
data$Longitude[data$Assembly.barcode %in% c("HEL_BA9244AA_AS")] <- 158.52
data$Latitude[data$Assembly.barcode %in% c("HEL_BA9070AA_AS", "HEL_BA9696AA_AS")] <- 61.36 ## Ket, Russia
data$Longitude[data$Assembly.barcode %in% c("HEL_BA9070AA_AS", "HEL_BA9696AA_AS")] <- 91.11
data$Latitude[data$Assembly.barcode %in% c("HEL_BA9234AA_AS")] <- 52.2 ## Orok, Russia
data$Longitude[data$Assembly.barcode %in% c("HEL_BA9234AA_AS")] <- 143.03
data$Latitude[data$Assembly.barcode %in% c("HEL_BA9437AA_AS")] <- 53.2 ## Tubalar, Russia
data$Longitude[data$Assembly.barcode %in% c("HEL_BA9437AA_AS")] <- 83.46
data$Latitude[data$Assembly.barcode %in% c("HEL_AA7275AA_AS")] <- 68.03 ## Chukch, Russia
data$Longitude[data$Assembly.barcode %in% c("HEL_AA7275AA_AS")] <- 166.26
data$Latitude[data$Assembly.barcode %in% c("HEL_BA9514AA_AS", "HEL_BA9559AA_AS")] <- 51.46 ## Nivkh, Russia
data$Longitude[data$Assembly.barcode %in% c("HEL_BA9514AA_AS", "HEL_BA9559AA_AS")] <- 143.07
data$Latitude[data$Assembly.barcode %in% c("HEL_BA9566AA_AS", "HEL_BA9807AA_AS")] <- 52.25 ## Tuvan (Todzha), Russia
data$Longitude[data$Assembly.barcode %in% c("HEL_BA9566AA_AS", "HEL_BA9807AA_AS")] <- 96.36
data$Latitude[data$Assembly.barcode %in% c("HEL_BA9184AA_AS")] <- 41.607 ## Kyrgystan -> center
data$Longitude[data$Assembly.barcode %in% c("HEL_BA9184AA_AS")] <- 74.751
data$Latitude[data$Assembly.barcode %in% c("HEL_BA9354AA_AS", "HEL_BA9653AA_AS")] <- 41.607 ## Kyrgyz, Kyrgystan -> center of Kyrgyzstan
data$Longitude[data$Assembly.barcode %in% c("HEL_BA9354AA_AS", "HEL_BA9653AA_AS")] <- 74.751
data$Latitude[data$Assembly.barcode %in% c("HEL_CA0045AA_AS", "HEL_CA0046AA_AS", "HEL_CA0047AA_AS", "HEL_CA0057AA_AS")] <- 39.78 ## USA -> center of the US
data$Longitude[data$Assembly.barcode %in% c("HEL_CA0045AA_AS", "HEL_CA0046AA_AS", "HEL_CA0047AA_AS", "HEL_CA0057AA_AS")] <- -100.45
data$Latitude[data$Assembly.barcode %in% c("HEL_CA0050AA_AS", "HEL_CA0051AA_AS")] <- 3.73 ## Colombia -> center
data$Longitude[data$Assembly.barcode %in% c("HEL_CA0050AA_AS", "HEL_CA0051AA_AS")] <- -73.125
data$Latitude[data$Assembly.barcode %in% c("HEL_CA0052AA_AS", "HEL_CA0053AA_AS", "HEL_CA0054AA_AS")] <- 43.25 ## Japan -> Ashoro
data$Longitude[data$Assembly.barcode %in% c("HEL_CA0052AA_AS", "HEL_CA0053AA_AS", "HEL_CA0054AA_AS")] <- 143.58

## also need to add the HpGP strains and the other absent strains from the enterobase table
## -> in their case, only have the country, based on their name
ix <- which(is.na(divpop$enterobase_name) & is.na(divpop$host_ethnicity))
data <- rbind(data, data.frame(Assembly.barcode = divpop$strain[ix], Country = divpop$country[ix], City = NA, Continent = divpop$continent[ix], Name = divpop$strain[ix], Host.Ethnicity = NA, Latitude = NA, Longitude  = NA, level = divpop$level[ix]))
data$Country[data$Assembly.barcode %in% c(paste0("CAN_", c("019", "095", "136", "172", "243", "258", "290", "219", "196")), "Ca_97_34")] <- "Canada"
data$Continent[data$Assembly.barcode %in% c(paste0("CAN_", c("019", "095", "136", "172", "243", "258", "290", "219", "196")), "Ca_97_34")] <- "North America"
data$Latitude[is.na(data$Latitude) & data$Country == "Canada"] <- 56 ## center
data$Longitude[is.na(data$Longitude) & data$Country == "Canada"] <- -101
data$Latitude[is.na(data$Latitude) & data$Country == "Chile"] <- -33.724 ## center
data$Longitude[is.na(data$Longitude) & data$Country == "Chile"] <- -71.015
data$Latitude[is.na(data$Latitude) & data$Country == "Colombia"] <- 3.73 ## center
data$Longitude[is.na(data$Longitude) & data$Country == "Colombia"] <- -73.125
data$Latitude[is.na(data$Latitude) & data$Country == "Mexico"] <- 23 ## center
data$Longitude[is.na(data$Longitude) & data$Country == "Mexico"] <- -102
data$Latitude[is.na(data$Latitude) & data$Country == "Peru"] <- -9.5 ## center
data$Longitude[is.na(data$Longitude) & data$Country == "Peru"] <- -74.62 
data$Latitude[is.na(data$Latitude) & data$Country == "Spain"] <- 40.45 ## center
data$Longitude[is.na(data$Longitude) & data$Country == "Spain"] <- -3.25
data$Latitude[is.na(data$Latitude) & data$Country == "USA"] <- 39.78 ## center
data$Longitude[is.na(data$Longitude) & data$Country == "USA"] <- -100.45
data$Latitude[is.na(data$Latitude) & data$Country == "Kazakhstan"] <- 48 ## center
data$Longitude[is.na(data$Longitude) & data$Country == "Kazakhstan"] <- 68
data$Latitude[is.na(data$Latitude) & data$Country == "Kyrgyzstan"] <- 41.607 ## center
data$Longitude[is.na(data$Longitude) & data$Country == "Kyrgyzstan"] <- 74.751

data$pop <- NA
for(i in unique(pop$pop)) {
    strains <- pop$strain[pop$pop == i]
    data$pop[data$Assembly.barcode %in% strains | data$Name %in% strains] <- i
}
data$pop[data$strain %in% c("HEL_AA0408AA_AS", "HEL_BA9262AA_AS")] <- "Hacinonychis"
data$pop[data$strain %in% c("HEL_BA9303AA_AS", "HEL_BA9707AA_AS")] <- "Primate"

data <- rbind(data, cbind(data0[data0$Assembly.barcode %in% pop$strain[pop$pop == "hpSahul"] | data0$Name %in% pop$strain[pop$pop == "hpSahul"],], level = "low", pop = "hpSahul"))

data$Latitude[is.na(data$Latitude) & data$pop == "hpSahul" & data$Country == "Papua New Guinea" & data$Host.Ethnicity == "Highlander"] <- -6.05 ## Papua New Guinea, Highlander
data$Longitude[is.na(data$Longitude) & data$pop == "hpSahul" & data$Country == "Papua New Guinea" & data$Host.Ethnicity == "Highlander"] <- 145.23
data$Latitude[is.na(data$Latitude) & data$pop == "hpSahul" & data$Country == "Australia" & data$City == "Perth" & is.na(data$Host.Ethnicity)] <- -31.9554 ## Australia, Perth
data$Longitude[is.na(data$Longitude) & data$pop == "hpSahul" & data$Country == "Australia" & data$City == "Perth" & is.na(data$Host.Ethnicity)] <- 115.8586

#write.table(data, file = "DRAFT/tables/listMap.csv", sep = ",", quote = FALSE, row.names = FALSE)

dataMonkey <- data.frame(Latitude = 38.53990, Longitude = -121.80483, Label = "Davis Primate \nCenter")


world <- map_data("world") #https://rdrr.io/cran/maps/man/world.html
baseMap <- ggplot() + 
    geom_polygon(data = world, aes(x = long, y = lat, group = group), color = "grey", fill = "white") +
    ylim(c(-55, 83)) +
    coord_fixed(1.3)

baseMap +
    geom_text(data = data[data$level == "high",], aes(x = Longitude, y = Latitude, label = Host.Ethnicity, color = pop), check_overlap = TRUE, hjust = 0, nudge_x = 3) +
    geom_point(data = data, aes(x = Longitude, y = Latitude, color = pop, shape = level, size = level), stroke = 0.5) +
    scale_shape_manual(values = c("high" = 1, "low" = 0, "outliers" = 18), labels = c("Hardy", "Ubiquitous", "outliers")) +
    scale_size_manual(values = c("high" = 6, "low" = 3, "outliers" = 2), labels = c("Hardy", "Ubiquitous", "outliers")) +
    scale_color_manual(values = colPop[names(colPop) %in% unique(data$pop)]) +
    geom_text(data = dataMonkey, aes(x = Longitude, y = Latitude, label = Label), check_overlap = TRUE, col = "black") +
    geom_point(data = dataMonkey, aes(x = Longitude, y = Latitude), col = "black") +
    labs(x = "longitude", y = "latitude", col = "population", shape = "ecotype", size = "ecotype") +
    theme(text = element_text(size = 18), axis.text = element_text(size = 18))
ggsave(filename = paste0("DRAFT/figures/map.pdf"), device = "pdf", width = 16, height = 9, units = "in", dpi = 500)
############

## TO CHECK
##### generate the pop supp table TS2 #####
## use what was already done for the map creation

listPop <- pop[c("strain", "enterobase_name", "divergence", "pop", "latitude", "longitude", "country", "continent", "host_ethnicity")]
colnames(listPop)[3] <- "ecotype"

listPop$ecotype[listPop$ecotype == "high"] <- "hardy"
listPop$ecotype[listPop$ecotype != "hardy"] <- "ubiquitous"
for(i in 1:nrow(data)) {
    straini <- data$Assembly.barcode[i]
    listPop$latitude[listPop$strain == straini | listPop$enterobase_name == straini] <- data$Latitude[i]
    listPop$longitude[listPop$strain == straini | listPop$enterobase_name == straini] <- data$Longitude[i]
}

listPop <- listPop[!(listPop$strain %in% strainsToRemove),]

listPop$analysis <- ""

## analysis done
## 1A - tree 3 pop
# tree <- read.nexus("TREE/tree")
# strain <- tree$tip.label
# listPop$analysis[listPop$strain %in% strain | listPop$enterobase_name %in% strain] <- sapply(listPop$analysis[listPop$strain %in% strain | listPop$enterobase_name %in% strain], function(x) paste(x, "1A", sep = ";"))
# strains in fineSTRUCTURE but not in the dataset:
# data[data$SampleID %in% strain[!(strain %in% listPop$strain) & !(strain %in% listPop$enterobase_name)],]
# SampleID                    fsPop
# 74     BURYAT49               hspSiberia
# 214    Nentsy82                  hspUral
# 454      SA172C                hpAfrica2
# 456      SA253C                hpAfrica2
# 483 MX_2004_102        hspAfrica1WAfrica
# 520     HN_G262 hspSWEuropeLatinAmerican
# 524     HN_G266 hspSWEuropeLatinAmerican
# 525     HN_G257 hspSWEuropeLatinAmerican
# 526     HN_G274 hspSWEuropeLatinAmerican
# 530       ES_10              hspSWEurope
# 531       ES_11              hspSWEurope
# 532       ES_15              hspSWEurope
# 535     Pt_7854              hspSWEurope
# 536        ES_7              hspSWEurope
# 537     Pt_7739              hspSWEurope
# 538       ES_13              hspSWEurope
# 540     Pt_8186              hspSWEurope
# 581     IL_3227      hspEuropeMiddleEast
# 585     IL_3199      hspEuropeMiddleEast
## 1A - GWAS
data <- read.table("GWAS/onAllGenome/input/pheno.txt", header = TRUE)
strain <- data$ID
listPop$analysis[listPop$strain %in% strain | listPop$enterobase_name %in% strain] <- sapply(listPop$analysis[listPop$strain %in% strain | listPop$enterobase_name %in% strain], function(x) paste(x, "1A", sep = ";"))
## 1B,C - tree high/low genes
strain <- read.nexus(paste0("DIVERGENT_NO_CDS/lowHighGenes/onAllGenome/trees/treeHigh"))$tip.label
listPop$analysis[listPop$strain %in% strain | listPop$enterobase_name %in% strain] <- sapply(listPop$analysis[listPop$strain %in% strain | listPop$enterobase_name %in% strain], function(x) paste(x, "1B;1C", sep = ";"))
## 3A - haplotype
strain <- read.table("METADATA/strains_pangenome.txt", header = TRUE)$Sample
listPop$analysis[listPop$strain %in% strain | listPop$enterobase_name %in% strain] <- sapply(listPop$analysis[listPop$strain %in% strain | listPop$enterobase_name %in% strain], function(x) paste(x, "3A", sep = ";"))
## 3B - nb high allele vs dist
load("HAPLOTYPES/onAllGenome/countHighAllele_pval10FST0.9.RData") ## subdata
data0 <- read.table("HAPLOTYPES/onAllGenome/nbdiff_lowLoadStrains.txt", header = TRUE, sep = "\t")
strain1 <- subdata$strain
strain2 <- data0$strain
strain <- as.character(strain1[strain1 %in% strain2])
strain <- unique(c(strain, data0$strainLow))
listPop$analysis[listPop$strain %in% strain | listPop$enterobase_name %in% strain] <- sapply(listPop$analysis[listPop$strain %in% strain | listPop$enterobase_name %in% strain], function(x) paste(x, "3B", sep = ";"))
## 3C - map
divpop <- read.table("METADATA/strains_SiberiaIndAm.csv", sep = ",", header = TRUE)$strain
outliers <- read.table("DRAFT/tables/listOutliers_nbAllDist.csv", sep = ",", header = TRUE)$strain
strain <- unique(c(divpop, outliers))
listPop$analysis[listPop$strain %in% strain | listPop$enterobase_name %in% strain] <- sapply(listPop$analysis[listPop$strain %in% strain | listPop$enterobase_name %in% strain], function(x) paste(x, "3C", sep = ";"))
## 4A - pangenome
strain <- read.table("METADATA/strains_pangenome.txt", header = TRUE)$Sample
listPop$analysis[listPop$strain %in% strain | listPop$enterobase_name %in% strain] <- sapply(listPop$analysis[listPop$strain %in% strain | listPop$enterobase_name %in% strain], function(x) paste(x, "4A", sep = ";"))
## 4B - Cag+/VacA/UreAB presence/absence
## HpGP strains used for the UreAB trees + 3 Cag+ strains + monkey strains + Hcetorum + Hacinonychis + 26695
strain <- read.table("GENES/UreB/UreB_pop.csv", header = TRUE, sep = ",")$strain.fragment
strain <- unique(sapply(strsplit(strain, split = ".", fixed = TRUE), function(x) x[1]))[c(512:565,5,312,413,84,305,296,395,271)]
listPop$analysis[listPop$strain %in% strain | listPop$enterobase_name %in% strain] <- sapply(listPop$analysis[listPop$strain %in% strain | listPop$enterobase_name %in% strain], function(x) paste(x, "4B", sep = ";"))
## S1 - tree whole pop
tree <- read.nexus("TREE/coreCDS")
strain <- tree$tip.label
listPop$analysis[listPop$strain %in% strain | listPop$enterobase_name %in% strain] <- sapply(listPop$analysis[listPop$strain %in% strain | listPop$enterobase_name %in% strain], function(x) paste(x, "S1", sep = ";"))
## S2 - PCA
strain <- data.frame(read_table("./PCA/onAllGenome/hpglobal_LD_PCA.eigenvec", col_names = FALSE))$X1
listPop$analysis[listPop$strain %in% strain | listPop$enterobase_name %in% strain] <- sapply(listPop$analysis[listPop$strain %in% strain | listPop$enterobase_name %in% strain], function(x) paste(x, "S2", sep = ";"))
## S3 - fineSTRUCTURE
data <- read.table("ANALYSIS_OTHER/fineSTRUCTURE_Roberto/run1/NAM_SIB_fsPops.txt", header = TRUE)
strain <- data$SampleID
listPop$analysis[listPop$strain %in% strain | listPop$enterobase_name %in% strain] <- sapply(listPop$analysis[listPop$strain %in% strain | listPop$enterobase_name %in% strain], function(x) paste(x, "S3", sep = ";"))
## S4 - GWAS
data <- read.table("GWAS/onAllGenome/input/pheno.txt", header = TRUE)
strain <- data$ID
listPop$analysis[listPop$strain %in% strain | listPop$enterobase_name %in% strain] <- sapply(listPop$analysis[listPop$strain %in% strain | listPop$enterobase_name %in% strain], function(x) paste(x, "S4", sep = ";"))
## S5 - ancestry profiles monkey / Sahul
listPop$analysis[listPop$pop %in% c("hpSahul", "Primate")] <- sapply(listPop$analysis[listPop$pop %in% c("hpSahul", "Primate")], function(x) paste(x, "S5", sep = ";"))
## S6 - dN/dS genes
data <- read.table(paste0("DIVERGENT_NO_CDS/lowHighGenes/onAllGenome/PAML/dNdS_between_high/YN00"), header = TRUE)
strain <- unique(c(data$strain1, data$strain2))
listPop$analysis[listPop$strain %in% strain | listPop$enterobase_name %in% strain] <- sapply(listPop$analysis[listPop$strain %in% strain | listPop$enterobase_name %in% strain], function(x) paste(x, "S6", sep = ";"))
## S7 - nb alleles vs nb blocks
load("HAPLOTYPES/onAllGenome/highBlock_pval10FST0.9.RData") ## highBlock
load("HAPLOTYPES/onAllGenome/countHighAllele_pval10FST0.9.RData") ## subdata
strain <- unique(highBlock$strain, subdata$strain)
listPop$analysis[listPop$strain %in% strain | listPop$enterobase_name %in% strain] <- sapply(listPop$analysis[listPop$strain %in% strain | listPop$enterobase_name %in% strain], function(x) paste(x, "S7", sep = ";"))
## S8 dotplots
strain <- c("HEL_AA0336AA_AS", "HEL_CA0045AA_AS", "HEL_AA0408AA_AS", "HEL_CA0054AA_AS", "HEL_AA0370AA_AS", "HEL_AA0367AA_AS", "HEL_AA4428AA_AS", "HEL_AA0364AA_AS", "HpGP-CAN-019", "HpGP-CHI-116", "HpGP-KAZ-002", "HpGP-ZAF-001")
listPop$analysis[listPop$strain %in% strain | listPop$enterobase_name %in% strain] <- sapply(listPop$analysis[listPop$strain %in% strain | listPop$enterobase_name %in% strain], function(x) paste(x, "S8", sep = ";"))
## S9 - tree H.cetorum
tree <- read.nexus("DIVERGENT_NO_CDS/Hcetorum_MIT99-5656/tree_Hcetorum/nwkFile/1_nt.nex")
strain <- unique(sapply(strsplit(tree$tip.label, split = ".", fixed = TRUE), function(x) x[1]))
listPop$analysis[listPop$strain %in% strain | listPop$enterobase_name %in% strain] <- sapply(listPop$analysis[listPop$strain %in% strain | listPop$enterobase_name %in% strain], function(x) paste(x, "S9", sep = ";"))
## S10 - tree VacA
tree <- read.nexus("GENES/VacA/Nucleotide_1fragHcetorum//HpGP_VacA.nex")
strain <- unique(sapply(strsplit(tree$tip.label, split = ".", fixed = TRUE), function(x) x[1]))
listPop$analysis[listPop$strain %in% strain | listPop$enterobase_name %in% strain] <- sapply(listPop$analysis[listPop$strain %in% strain | listPop$enterobase_name %in% strain], function(x) paste(x, "S10", sep = ";"))
## S11 - tree UreA
tree <- read.nexus("GENES/UreA/Nucleotide_sample/UreA_nt.nex")
strain <- unique(sapply(strsplit(tree$tip.label, split = ".", fixed = TRUE), function(x) x[1]))
listPop$analysis[listPop$strain %in% strain | listPop$enterobase_name %in% strain] <- sapply(listPop$analysis[listPop$strain %in% strain | listPop$enterobase_name %in% strain], function(x) paste(x, "S11", sep = ";"))
## S12 - tree UreB
tree <- read.nexus("GENES/UreB/Nucleotide_sample/UreB_nt.nex")
strain <- unique(sapply(strsplit(tree$tip.label, split = ".", fixed = TRUE), function(x) x[1]))
listPop$analysis[listPop$strain %in% strain | listPop$enterobase_name %in% strain] <- sapply(listPop$analysis[listPop$strain %in% strain | listPop$enterobase_name %in% strain], function(x) paste(x, "S12", sep = ";"))


## at the end of the table, add
## MIT 99-5656,NA,NA,Hcetorum,NA,NA,;4B;S9;S10;S11;S8,"genome assembly ASM25927v1 NCBI RefSeq sequence GCF_000259275.1 https://www.ncbi.nlm.nih.gov/data-hub/genome/GCF_000259275.1/ downloaded the 10/02/2023"

write.table(listPop, "DRAFT/tables/listStrains.csv", sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE)

## add the H. cetorum strain manually
## + the additional information necessary (which strain etc...)

##########

#### OK
###### 16. dS against the Monkey strains ######

## dN/dS vs dS
data <- read.table(paste0("PAML/betweenMonkey/monkeySahul/YN00"), header = TRUE)
colnames(data) <- c("strain", colnames(data)[-1])

data <- data[!(data$strain %in% strainsToRemove),]

data <- getAllPop(data, pop)
data <- aggregate(data[c("omega", "dS")], by = list(strain = data$strain, pop = data$pop, divergence = data$divergence), mean)
data$divergence[data$divergence == "other"] <- "low"
dataM <- aggregate(data[c("omega", "dS")], by = list(pop = data$pop, divergence = data$divergence), mean)
data$pop <- factor(data$pop, levels = unique(dataM$pop[order(dataM$dS)]))
ggplot(data) +
    geom_boxplot(aes(y = dS, x = pop, col = divergence)) +
    scale_color_manual(values = c("red", "black"), labels = c("Hardy", "Ubiquitous")) +
    labs(x = "population", y = "dS", col = "ecotype") +
    theme_set(theme_bw()) +
    theme(text = element_text(size = 18), axis.text = element_text(size = 14), axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))
ggsave(filename = paste0("DRAFT/figures/set_allPop/dS_monkey.pdf"), device = "pdf", width = 16, height = 9, units = "in", dpi = 500, scale = 0.75)


##

############

#### OK
###### 17. tree H. cetorum ######

for(gene in 1:81) {
    geneName <- strsplit(system(paste0("grep 'GENE ", gene, " =====>' outTreeHits.out"), intern = TRUE), split = "=====> ")[[1]][2]
    tree <- read.nexus(paste0("rooted_nexus/genes/", gene))
    tree <- drop.tip(tree, strainsToRemove)
    div <- data.frame(straini = tree$tip.label, pop = NA, divergence = "other")
    div$strain <- sapply(strsplit(div$straini, split = ".", fixed = TRUE), function(x) x[1])
    div <- getAllPop(div, pop)
    div$pop[div$pop == "NA"] <- "Hcetorum"
    div$divergence[div$divergence == "other"] <- "low"
    div$divergence[div$pop == "Hcetorum"] <- "Hcetorum"
    div$divergence[div$pop == "host"] <- "whales"
    p <- ggtree(tree, layout = "rectangular") +
        geom_treescale() 
    p2 <- p %<+% div +
        aes(color = pop) +
        geom_tippoint(aes(color = pop, shape = interaction(divergence, host)), size = 1) +
        #geom_tiplab(geom = "text", align = TRUE) +
        scale_color_manual(values = colPop[names(colPop) %in% unique(div$pop)]) +
        scale_shape_manual(values = c(1, 15, 19, NA, 0), labels = c("Hardy (H. acinonychis)", "H. cetorum",  "Hardy", "Ubiquitous", "Hardy (Primates)")) +
        #theme(plot.margin = unit(c(-150,-50,-60,-200), "mm")) +
        # hexpand(.05, direction = -1) +
        labs(color = "population", fill = "population", shape = "ecotype", subtitle = geneName) +
        theme(text = element_text(size = 12), legend.position = "none")  
    ggsave(p2, filename = paste0("DRAFT/figures/TREE/TREE_", gene, ".pdf"), device = "pdf", width = 16, height = 9, units = "in", dpi = 500)
}

############















