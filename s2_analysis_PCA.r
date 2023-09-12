## plot the PCA results

setwd("~/Desktop/HpGlobal_noHpGP")

library(tidyverse)
pca <- read_table("./PCA/fullAll/hpglobal_LD_PCA.eigenvec", col_names = FALSE)
eigenval <- scan("./PCA/fulAll/hpglobal_LD_PCA.eigenval")

# sort out the pca data
# remove nuisance column
pca <- pca[,-1]
# set names
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))


pop <- read.table("./METADATA/strains_pop.csv", sep = ",", header = TRUE)

pca$pop <- NA
pca$fspop <- NA
for(i in unique(pop$MergedPop)) {
  strains <- pop$Sample[pop$MergedPop == i]
  pca$pop[pca$ind %in% strains] <- i
}
for(i in unique(pop$pop)) {
  strains <- pop$Sample[pop$pop == i]
  pca$fspop[pca$ind %in% strains] <- i
}
pca <- pca[!is.na(pca$pop),]

# first convert to percentage variance explained
pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)

a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_light()

# plot pca
ggplot(pca, aes(PC1, PC2, col = fspop, fill = fspop)) + 
  geom_point() +
  xlab(paste0("PC1 (", round(pve$pve[1],1), "%)")) +
  ylab(paste0("PC2 (", round(pve$pve[2],1), "%)"))
ggsave("PLOT/PCA/PCAcoreCDS_africa2.jpg", height = 14, width = 14)

ggplot(pca, aes(PC3, PC4, col = fspop, fill = fspop)) + 
  geom_point() +
  xlab(paste0("PC3 (", round(pve$pve[3],1), "%)")) +
  ylab(paste0("PC4 (", round(pve$pve[4],1), "%)"))
ggsave("PLOT/PCA/PCA2coreCDS_africa2.jpg", height = 14, width = 14)

pca2 <- pivot_longer(pca, cols = 3:6, names_to = "PC", values_to = "vec")
ggplot(pca2[pca2$fspop %in% c("hspIndigenousNAmerica", "hspIndigenousSAmerica", "hspSiberia"),], aes(x = PC1, y = vec, col = fspop)) +
geom_point() +
facet_wrap(~ PC, scales = "free")


