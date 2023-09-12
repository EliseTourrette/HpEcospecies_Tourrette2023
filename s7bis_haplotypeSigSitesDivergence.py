#!/usr/bin/python

## in HpGlobal_noHpGP

## look at the haplotypes of subset of sites
## based on a file that contains the position that must be checked
## three criteria:
## FST > 0.9
## -log(p) > 10 (GWAS)
## and sites with both threshold

## first, need all these specific sites in one file
## one file for each criterion
import numpy as np

setthr = "pval10FST0.9"

## these are the true positions
## in the global alignment
pos = list(np.loadtxt("HAPLOTYPES/pos_" + setthr + ".txt", dtype = "int"))


## for now, we only want to look at the differences between the low and high load strains
## based on their alleles
## no need to take the reverse complement (for the visualization) as we are only interested in the differences (if the nucleotide is different from the "low" allele - for example -, then it will also be the case for the reverse complement)
## I will pay attention to taking the reverse complement when the CDS is read on the reverse strand when looking at the CDS sequence in detail (for example, the type of nucleotides that are different)
strain = []
sequence = []
with open("DATA/woHpGP_dnds.cvg70.aln", 'r') as f: ## if using the whole aln
#with open("PREPROCESSING/coreCDSclean.aln", 'r') as f: ## if using the core aln
    for l in f:
        iline = l.rstrip('\n').split()
        if(iline[0][0] == '>'):
            strain.append(iline[0])
        else:
            ## index -1: python: i0 = 0; while pos0 = 1
            sequence.append([iline[0][i - 1] for i in pos]) ## if using the whole aln
            #sequence.append([iline[0][i] for i in posCoreix]) ## if using the core aln


pos = list(map(str, pos))
with open("HAPLOTYPES/sites_" + setthr + ".txt", 'w') as f:
    f.write("strain" + '\t' + '\t'.join(pos) + '\n')
    for i in range(0, len(sequence)):
        f.write(strain[i] + '\t' + '\t'.join(sequence[i]) + '\n')


