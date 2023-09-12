# HpEcospecies_Tourrette2023
repository for the analysis done for the hpEcospecies paper (preprint version 2023)

An ancient ecospecies of Helicobacter pylori found in Indigenous populations and animal adapted lineages
Elise Tourrette, Roberto C. Torres, Sarah L. Svensson, Takashi Matsumoto, Muhammad Miftahussurur, Kartika Afrida Fauzia, Ricky Indra Alfaray, Ratha-Korn Vilaichone, Vo Phuoc Tuan, HelicobacterGenomicsConsortium, Difei Wang, Abbas Yadegar, Lisa M. Olsson, Zhemin Zhou, Yoshio Yamaoka, Kaisa Thorell, Daniel Falush
bioRxiv 2023.04.28.538659; doi: https://doi.org/10.1101/2023.04.28.538659

For the data analysis:

s0_workflow.sh: main script that describe the whole workflow (order of the scripts used) and the different command lines used

## REQUIREMENTS
###############
### SOFTWARES
## PAML 4.9 [calculation of dN/dS and dS values -> yn00 module]
## snp_sites 2.5.1 (installed via conda) [transformation of a multi-FASTA file into a VCF file]
## PLINK 1.9 (on the cluster) [LD-pruning and PCA]
## gemma 0.93 (executable file given in the github bugwas repository) [GWAS]
## FastTree 2.1.10 (cluster) [phylogenetic tree]
## MAFFT 7.505 [multi-sequence alignment]
## BLAST 2.11.0+ (blastn)

### PYTHON LIBRARIES
## python 3.10
## numpy
## Bio (SeqIO)
## subprocess
## sys

### R PACKAGES
## R 4.1.1
## bugwas 0.0.0.9000 [GWAS; https://www.nature.com/articles/nmicrobiol201641]
## PopGenome
## ggplot2 3.3.6
## ape 5.7-1
## ggtree 3.2.1
## tidyverse 1.3.2
## reshape2 1.4.4
## cowplot 1.1.1
## gridExtra 2.3
## ggpubr 0.4.0
###############



