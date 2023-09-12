#!/usr/bin/python

import numpy as np
import sys

ixfocus = sys.argv[1]

listStrain = np.loadtxt("METADATA/listStrains.txt", dtype = "str")
listStrain = listStrain[int(ixfocus):(int(ixfocus)+700)]

strainData = np.loadtxt("METADATA/listStrains.txt", dtype = "str")

meta = np.loadtxt("METADATA/strains_SiberiaIndAm.csv", dtype = "str", delimiter = ",")

low = meta[np.where(meta[:,8] == "low")[0],0]
low = np.array([i for i in low if np.count_nonzero(strainData == i)  > 0])


## read the alignment file
strain = []
sequence = []
with open("DATA/woHpGP_dnds.cvg70.aln", 'r') as f:
    for l in f:
        iline = l.rstrip('\n').split()
        if(iline[0][0] == '>'):
            strain.append(iline[0][1:])
        else:
            sequence.append(iline[0])

strain = np.array(strain)

if(ixfocus == '0'):
    with open("HAPLOTYPES/nbDiff_" + str(ixfocus) + ".txt", 'a') as f:
            f.write('strain\tstrainLow\tnbDiff\tnbDiffGaps\tnbDiffNoGap\tnbSite\tnbGaps\n')

for focus in listStrain:
    focusIn = np.where(low == focus)[0]
    ## only calculate the genetic distance if the strain doesn't have a low load (ie focusIn.size == 0)
    ## or not...
    if focusIn.size >= 0:
        ## load the sequence of the "focus" strain
        if np.count_nonzero(strainData == focus) > 0:
            seqFocus = sequence[int(np.where(strain == focus)[0])]
            for i in range(0, low.size):
                seqLow = sequence[int(np.where(strain == low[i])[0])]
                nbdiff = [focus, low[i], str(sum(x != y for x,y in zip(seqFocus, seqLow)))]
                ## when calculating the number of sites that are different between the two strains
                ## do not take into account the gaps
                ## PAML, GWAS and FST also remove the sites with gaps
                ## thus, when calculating the percentage of difference, we need to divide the number of different sites by the total number of sites WITHOUT gaps
                aa = [ [seqFocus[j], seqLow[j]] for j in range(0, len(seqFocus)) if seqFocus[j] != seqLow[j]]
                bb = [j for j in aa if j[0] == '-' or j[1] == '-']
                cc = [j for j in range(0, len(seqFocus)) if seqFocus[j] == '-' or seqLow[j] == '-']
                with open("HAPLOTYPES/nbDiff_" + str(ixfocus) + ".txt", 'a') as f:
                    f.write('\t'.join(nbdiff) + '\t' + str(len(bb)) + '\t' + str(len(aa) - len(bb)) + '\t' + str(len(seqFocus)) + '\t' + str(len(cc)) + '\n')

   

