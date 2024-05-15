
## to put in the general plot script!!!

## 1) calculate the average distance between each population that has both ecotypes


library(ape)

listStrainHistory <- read.csv("/media/eliset/PortableSSD/2021_2023_Postdoc_IPS/HpGlobal/wHpGP/TREE/evolutionaryHistory/listStrainEvol.csv")

calcDistPop <- function(popNames, strainsPop, data) {
    distMat <- dist.dna(data, as.matrix = TRUE, model = "raw")
    distPop <- matrix(0, nrow = length(popNames), ncol = length(popNames))
    colnames(distPop) <- rownames(distPop) <- popNames
    
    for(i in 1:length(popNames)) {
        for(j in 1:length(popNames)) {
            if(i != j) {
                tmp <- distMat[rownames(distMat) %in% strainsPop[[i]], colnames(distMat) %in% strainsPop[[j]]]
                distPop[i,j] <- mean(tmp)
            }
        }
    } 
    return(distPop)
}

histDistPop <- function(data, listStrainHistory) {
    listStrainHistory$pop[listStrainHistory$pop == "Hacinonychis"] <- "hpAfrica2"
    names(data) <- unlist(lapply(names(data), function(x) listStrainHistory$pop[listStrainHistory$strain == x]))
    distMat <- dist.dna(data, as.matrix = TRUE, model = "raw")
    distPop <- data.frame(pop = NULL, dist = NULL)
    for(i in 1:(nrow(distMat)-1)) {
        for(j in (i+1):ncol(distMat)) {
            compPop <- c(rownames(distMat)[i], colnames(distMat)[j])
            distPop <- rbind(distPop, data.frame(pop = paste(sort(compPop), collapse = "\n"), dist = distMat[i,j]))
        }
    }
    return(distPop)
}


strainsPop <- list(afr2 = listStrainHistory$strain[listStrainHistory$pop %in% c("Hacinonychis", "hpAfrica2")],
                   sahul = listStrainHistory$strain[listStrainHistory$pop %in% c("Primate", "hpSahul")],
                   sib = listStrainHistory$strain[listStrainHistory$pop == "hspSiberia"],
                   nam = listStrainHistory$strain[listStrainHistory$pop == "hspIndigenousNAmerica"],
                   sam = listStrainHistory$strain[listStrainHistory$pop == "hspIndigenousSAmerica"])
popNames <- c("Africa", "Sahul", "Siberia", "Indigenous NAmerica", "Indigenous SAmerica")

# distMat <- read.table("hardyDiff.dist", skip = 1)
# strains <- distMat$V1
# distMat <- as.matrix(distMat[,-1])
# colnames(distMat) <- rownames(distMat) <- strains
data <- read.FASTA("/media/eliset/PortableSSD/2021_2023_Postdoc_IPS/HpGlobal/wHpGP/TREE/evolutionaryHistory/hardyDiff.aln")
distPop <- calcDistPop(popNames, strainsPop, data)
tree1_0 <- nj(as.dist(distPop)) 
## root the tree: outgroup = africa
tree1 <- root(tree1_0, outgroup = "Africa", resolve.root = TRUE)
a1 <- histDistPop(data, listStrainHistory)

data <- read.FASTA("/media/eliset/PortableSSD/2021_2023_Postdoc_IPS/HpGlobal/wHpGP/TREE/evolutionaryHistory/ubiDiff.aln")
distPop <- calcDistPop(popNames, strainsPop, data)
tree2_0 <- nj(as.dist(distPop)) 
## root the tree: outgroup = africa
tree2 <- root(tree2_0, outgroup = "Africa", resolve.root = TRUE)
a2 <- histDistPop(data, listStrainHistory)

data <- read.FASTA("/media/eliset/PortableSSD/2021_2023_Postdoc_IPS/HpGlobal/wHpGP/TREE/evolutionaryHistory/undiff.aln")
distPop <- calcDistPop(popNames, strainsPop, data)
tree3_0 <- nj(as.dist(distPop)) 
## root the tree: outgroup = africa
tree3 <- root(tree3_0, outgroup = "Africa", resolve.root = TRUE)
a3 <- histDistPop(data, listStrainHistory)

library(ggtree)
ggtree(tree1) + geom_tiplab()
ggsave(file = "hardyDiff.pdf")
ggtree(tree2) + geom_tiplab()
ggsave(file = "ubiDiff.pdf")
ggtree(tree3) + geom_tiplab()
ggsave(file = "undiff.pdf")

par(mfrow = c(2,2))
pdf("hardyDiff.pdf")
plot(tree1)
dev.off()
pdf("ubiDiff.pdf")
plot(tree2)
dev.off()
pdf("undiff.pdf")
plot(tree3)
dev.off()

library(ggplot2)
distPop <- rbind(cbind(a1, set = "hardy_diff"), cbind(a2, set = "ubi_diff"), cbind(a3, set = "undiff"))
ggplot(distPop) +
    geom_boxplot(aes(x = pop, y = dist, col = set))
ggsave(file = "distribDist.pdf", height = 15, width = 25, units = "in", scale = 0.75)




