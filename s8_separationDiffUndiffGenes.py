## script that separate the genes into high and low divergent genes
## based on the FST and GWAS p-values of the segregating sites
## low divergent genes: no sites with high FST / significant sites
## high divergent genes: contains more than x% of the sites with high FST and significant for the GWAS

## create two separate multi FASTA alignment files
## these files will be used to generate two phylogenetic trees (FastTree software) based on the subset of core genes 

## note that the SNPs are not necessarily from CDS but can also come from the intergenic region

import numpy as np

## read the positions of the sites with -log(pval) > 10 and FST > 0.9
pos = list(np.loadtxt("HAPLOTYPES/pos_pval10FST0.9.txt", dtype = "int"))
## as well as the positions of the non-significant SNPs, with -log(pval) < 10 and FST <0.5
posNonSig = list(np.loadtxt("HAPLOTYPES/pos_nonSigSites_pval10FST0.5.txt", dtype = "int"))
## and the other sites
posInter = list(np.loadtxt("HAPLOTYPES/pos_otherSNPs.txt", dtype = "int"))

## read the alignment file
strain = []
sequence = []
with open("DATA/woHpGP_dnds.cvg70.aln", 'r') as f:
    for l in f:
        iline = l.rstrip('\n').split()
        if(iline[0][0] == '>'):
            strain.append(iline[0])
        else:
            sequence.append(iline[0])

nbSeq = len(sequence)

## separate the CDS into high divergence and low divergence ones
## based on the number of (non-)significant sites each CDS has (more than 5 sig SNPs -> high divergence; all non-sig SNPs (ie no sig and intermediate sites) -> low divergence; all the other genes aree put in the intermediate genes)
## note that we take the reverse complement when the CDS is on the reverse strand (strand = '-'), necessary for the dN/dS analysis
## take the positions (i1 and i2) minus 1 to transform the positions into index
highCDS = []
lowCDS = []
interCDS = []
## number of significant (GWAS, FST) sites for each gene
lenSig = []
## annotation for low/high divergent genes
geneHigh = []
geneLow = []
geneInter = []
## sequence, for all strains, of the different CDS
highSeq = [[] for i in range(0, nbSeq)]
lowSeq = [[] for i in range(0, nbSeq)]
interSeq = [[] for i in range(0, nbSeq)]
thr = 5
ix = 0 
with open("PREPROCESSING/posCDS.csv", "r") as f:
    for l in f:
        print(ix); ix += 1
        iline = l.rstrip('\n').split("\t")
        i1 = int(iline[3]) - 1
        i2 = int(iline[4]) - 1
        if iline[6] == '+':
            ixCDS = list(range(i1, i2+1))
            sigPos = [i for i in ixCDS if (i + 1) in pos]
            nonSigPos = [i for i in ixCDS if (i + 1) in posNonSig]
            otherSites = [i for i in ixCDS if (i + 1) in posInter]
            if len(sigPos) == 0:
                if len(otherSites) == 0:
                    lowCDS.append(ixCDS)
                    geneLow.append(iline)
                    for i in range(0, nbSeq):
                        seq = sequence[i][i1:(i2+1)]
                        lowSeq[i].append(seq)
                else:
                    interCDS.append(ixCDS)
                    geneInter.append(iline)
                    for i in range(0, nbSeq):
                        seq = sequence[i][i1:(i2+1)]
                        interSeq[i].append(seq)
            elif len(sigPos) < thr:
                interCDS.append(ixCDS)
                geneInter.append(iline)
                for i in range(0, nbSeq):
                    seq = sequence[i][i1:(i2+1)]
                    interSeq[i].append(seq)
            elif len(sigPos) >= thr:
                highCDS.append(ixCDS)
                lenSig.append(len(sigPos))
                geneHigh.append(iline)
                for i in range(0, nbSeq):
                    seq = sequence[i][i1:(i2+1)]
                    highSeq[i].append(seq)
        if iline[6] == '-':
            ixCDS = list(range(i2, i1-1, -1))
            sigPos = [i for i in ixCDS if (i + 1) in pos]
            nonSigPos = [i for i in ixCDS if (i + 1) in posNonSig]
            otherSites = [i for i in ixCDS if (i + 1) in posInter]
            if len(sigPos) == 0:
                if len(otherSites) == 0:
                    lowCDS.append(ixCDS)
                    geneLow.append(iline)
                    for i in range(0, nbSeq):
                        seq = sequence[i][(i2):(i1-1):-1]
                        seqRev = ["T" if seq[iseq] == "A" else "A" if seq[iseq] == "T" else "C" if seq[iseq] == "G" else "G" if seq[iseq] == "C" else "-" for iseq in range(0, len(seq))]
                        lowSeq[i].append(''.join(seqRev))
                else:
                    interCDS.append(ixCDS)
                    geneInter.append(iline)
                    for i in range(0, nbSeq):
                        seq = sequence[i][(i2):(i1-1):-1]
                        seqRev = ["T" if seq[iseq] == "A" else "A" if seq[iseq] == "T" else "C" if seq[iseq] == "G" else "G" if seq[iseq] == "C" else "-" for iseq in range(0, len(seq))]
                        interSeq[i].append(''.join(seqRev))
            elif len(sigPos) < thr:
                interCDS.append(ixCDS)
                geneInter.append(iline)
                for i in range(0, nbSeq):
                    seq = sequence[i][(i2):(i1-1):-1]
                    seqRev = ["T" if seq[iseq] == "A" else "A" if seq[iseq] == "T" else "C" if seq[iseq] == "G" else "G" if seq[iseq] == "C" else "-" for iseq in range(0, len(seq))]
                    interSeq[i].append(''.join(seqRev))
            elif len(sigPos) >= thr:
                highCDS.append(ixCDS)
                lenSig.append(len(sigPos))
                geneHigh.append(iline)
                for i in range(0, nbSeq):
                    seq = sequence[i][(i2):(i1-1):-1]
                    seqRev = ["T" if seq[iseq] == "A" else "A" if seq[iseq] == "T" else "C" if seq[iseq] == "G" else "G" if seq[iseq] == "C" else "-" for iseq in range(0, len(seq))]
                    highSeq[i].append(''.join(seqRev))

###





###
## remove the strains with only gaps (do not have the gene)
missing = []
with open("DIVERGENT_NO_CDS/highCoreCDS.aln", 'w') as f:
    for i in range(0, len(strain)):
        tmp = [k for j in highSeq[i] for k in j]
        if list(set(tmp)) != ['-']:
            f.write(strain[i] + '\n')
            f.write(''.join(tmp) + '\n')
        else:
            missing.append(strain[i])


missing = []
with open("DIVERGENT_NO_CDS/lowCoreCDS.aln", 'w') as f:
    for i in range(0, len(strain)):
        tmp = [k for j in lowSeq[i] for k in j]
        if list(set(tmp)) != ['-']:
            f.write(strain[i] + '\n')
            f.write(''.join(tmp) + '\n')
        else:
            missing.append(strain[i])


missing = []
with open("DIVERGENT_NO_CDS/interCoreCDS.aln", 'w') as f:
    for i in range(0, len(strain)):
        tmp = [k for j in interSeq[i] for k in j]
        if list(set(tmp)) != ['-']:
            f.write(strain[i] + '\n')
            f.write(''.join(tmp) + '\n')
        else:
            missing.append(strain[i])

###

## for the high divergence genes
## write one multi-FASTA file per gene (name = gene_beginning_end)
for i in range(0, len(highSeq[0])):
    with open("DIVERGENT_NO_CDS/highGenesSeq/gene_" + geneHigh[i][3] + "-" + geneHigh[i][4] + ".aln", "w") as f:
        for j in range(0, len(highSeq)):
            f.write(strain[j] + '\n')
            f.write(highSeq[j][i] + '\n')

###        
## for PAML, we also need to have one fasta file per sequence
## plus, need to remove the stop (TAG, TGA, TAA) codons from the alignments
highSeq0 = highSeq
highSeq = [[k for j in i for k in j] for i in highSeq0]
CDSfilter = [[] for i in range(0,len(highSeq))]
for i in range(0, len(highSeq[0]), 3):
    seq = [j[i] + j[i+1] + j[i+2] for j in highSeq]
    if ((seq.count('TAG') == 0) & (seq.count('TAA') == 0) & (seq.count('TGA') == 0)):
        for j in range(0,len(highSeq)):
            CDSfilter[j].append(seq[j])

for i in range(0, len(strain)) :
    with open("DIVERGENT_NO_CDS/highDivStrains/" + strain[i][1:] + ".fasta", 'w') as f:
        f.write(strain[i] + '\n')
        f.write(''.join(CDSfilter[i]) + '\n')


lowSeq0 = lowSeq
lowSeq = [[k for j in i for k in j] for i in lowSeq0]
CDSfilter = [[] for i in range(0,len(lowSeq))]
for i in range(0, len(lowSeq[0]), 3):
    print(i)
    seq = [j[i] + j[i+1] + j[i+2] for j in lowSeq]
    if ((seq.count('TAG') == 0) & (seq.count('TAA') == 0) & (seq.count('TGA') == 0)):
        for j in range(0,len(lowSeq)):
            CDSfilter[j].append(seq[j])

for i in range(0, len(strain)) :
    with open("DIVERGENT_NO_CDS/lowDivStrains/" + strain[i][1:] + ".fasta", 'w') as f:
        f.write(strain[i] + '\n')
        f.write(''.join(CDSfilter[i]) + '\n')


interSeq0 = interSeq
interSeq = [[k for j in i for k in j] for i in interSeq0]
CDSfilter = [[] for i in range(0,len(interSeq))]
for i in range(0, len(interSeq[0]), 3):
    print(i)
    seq = [j[i] + j[i+1] + j[i+2] for j in interSeq]
    if ((seq.count('TAG') == 0) & (seq.count('TAA') == 0) & (seq.count('TGA') == 0)):
        for j in range(0,len(interSeq)):
            CDSfilter[j].append(seq[j])

for i in range(0, len(strain)) :
    with open("DIVERGENT_NO_CDS/interDivStrains/" + strain[i][1:] + ".fasta", 'w') as f:
        f.write(strain[i] + '\n')
        f.write(''.join(CDSfilter[i]) + '\n')

###    
        
## print the list of genes 
with open("DIVERGENT_NO_CDS/geneHigh.txt", "w") as f:
    for i in range(0, len(geneHigh)):
        f.write('\t'.join(geneHigh[i]) + '\n')

with open("DIVERGENT_NO_CDS/geneLow.txt", "w") as f:
    for i in range(0, len(geneLow)):
        f.write('\t'.join(geneLow[i]) + '\n')

with open("DIVERGENT_NO_CDS/geneInter.txt", "w") as f:
    for i in range(0, len(geneInter)):
        f.write('\t'.join(geneInter[i]) + '\n')

## and the number of significant sites per genes
with open("DIVERGENT_NO_CDS/nbSigSitesHigh.txt", "w") as f:
    for i in range(0, len(geneHigh)):
        f.write('\t'.join(geneHigh[i]) + ' \t' + str(lenSig[i]) + '\n')


