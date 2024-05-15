library(readxl)

setwd("/media/eliset/PortableSSD/2021_2023_Postdoc_IPS/HpGlobal/wHpGP")

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


## get the dN/dS values

gene <- "high"
data <- read.table(paste0("DIVERGENT_NO_CDS/lowHighGenes/onAllGenome/PAML/dNdS_", gene, "_HardyUbiquitous/YN00"), header = TRUE)
colnames(data) <- c("strain", colnames(data)[-1])

data <- getAllPop(data, pop)
data <- data[data$pop != "NA",]

## ANI
#data <- read_excel("/media/eliset/PortableSSD/HP_lastVersions/fastANI_Hardy_vs_all.xlsx", sheet = "fastANI_undifferentiated_Hardy_")
#mean(data$ANI[which(data$Population == "hpAfrica2" & data$Ecospecies == "Hardy")])

ani <- read_excel("/media/eliset/PortableSSD/HP_lastVersions/fastANI_Hardy_vs_all.xlsx", sheet = "differentiated sites")


outgroups <- unique(pop$strain[pop$pop == "hpAfrica2" & pop$divergence != "high"], pop$enterobase_name[pop$pop == "hpAfrica2" & pop$divergence != "high"])

data$ANI <- NA
for(i in unique(ani$Query)) {
    print(i)
    for(j in unique(ani$Reference)) {
        if(length(which(data$strain == i & data$strain2 == j)) != 0) data$ANI[which(data$strain == i & data$strain2 == j)] <- ani$ANI[which(ani$Query == i & ani$Reference == j)]    
        if(length(which(data$strain2 == i & data$strain == j)) != 0) data$ANI[which(data$strain2 == i & data$strain == j)] <- ani$ANI[which(ani$Query == i & ani$Reference == j)]
    }
}

data$outgroup <- "Hardy"
data$outgroup[data$strain %in% c("HEL_AA0408AA_AS", "HEL_BA9262AA_AS")] <- "Hacinonychis"

ggplot(data) +
    geom_point(aes(x = dS, y = ANI, col = pop), size = 3, alpha = 0.75) +
    labs(x = "dS", y = "ANI", col = "population") +
    theme_set(theme_bw()) +
    theme(text = element_text(size = 18), axis.text = element_text(size = 14)) +
    facet_wrap(~ outgroup) +
    guides(alpha = "none")


    





