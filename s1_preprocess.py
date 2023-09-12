#!/usr/bin/python3

## python 3.10.0
## script to get only the CDS from the aligned genomes
## plus only the CDS that are in the core genome

## read the alignment file
strain = []
sequence = []
CDS = []
with open("../DATA/woHpGP_dnds.cvg70.aln", 'r') as f:
    for l in f:
        iline = l.rstrip('\n').split()
        if(iline[0][0] == '>'):
            strain.append(iline[0])
        else:
            sequence.append(iline[0])
            CDS.append([])


nbSeq = len(sequence)

## only keep the sites that are part of the CDS
## pos - 1 as in python, the first index is 0 while the position of the first site should be 1
## note that I also put the codon in the "correct" reading direction
## ie if strand (column 6) = -, inverse the direction of the codon (ALSO NEED TO USE THE COMPLEMENTARY NUCLEOTIDE)
## as it will impact the amino acid the result from the codon, it will directly impact the dN and dS values
## also need to pay attention the phase of the CDS (column 8); if not CDS, phase = '.'
## For features of type "CDS", the phase indicates where the feature begins with reference to the reading frame. The phase is one of the integers 0, 1, or 2, indicating the number of bases that should be removed from the beginning of this feature to reach the first base of the next codon.
## note that in HpGlobal, all CDS have phase = 0
posCDS = []
ix = 0
with open("posCDS.csv", "r") as f:
    for l in f:
        print(ix); ix += 1
        iline = l.rstrip('\n').split("\t")
        i1 = int(iline[3]) - 1
        i2 = int(iline[4]) - 1
        if iline[6] == '+':
            posCDS.append(list(range(i1, i2+1)))
            for i in range(0, nbSeq):
                seq = sequence[i][i1:(i2+1)]
                CDS[i].append(seq)
        if iline[6] == '-':
            posCDS.append(list(range(i2, i1-1, -1)))
            for i in range(0, nbSeq):
                seq = sequence[i][(i2):(i1-1):-1]
                seqRev = ["T" if seq[iseq] == "A" else "A" if seq[iseq] == "T" else "C" if seq[iseq] == "G" else "G" if seq[iseq] == "C" else "-" for iseq in range(0, len(seq))]
                CDS[i].append(''.join(seqRev))


## save the gene sequence of each CDS separately
for i in range(0, len(CDS[0])):
    with open("CDS/" + str(i) + ".aln", "w") as f:
        for j in list(range(0,nbSeq)):
            f.write(strain[j] + "\n")
            [f.write(k) for k in CDS[j][i]]
            f.write("\n")

CDS = [[k for j in i for k in j] for i in CDS]
posCDS = [j for i in posCDS for j in i]
nbCDS = len(CDS[0])

## remove the codons that are not present in at least 99% of the sequences
## calculate the minimum number of sequences needed
th = int(nbSeq/100)
CDSfilter = [[] for i in range(0,nbSeq)]
posCDScore = []
for i in range(0, nbCDS, 3):
    ## group the sites by codons
    ## get the same codon for all sequences
    seq = [j[i] + j[i+1] + j[i+2] for j in CDS] 
    print(i)
    if (seq.count('---') <= th) and (seq.count('TAG') == 0) and (seq.count('TAA') == 0) and (seq.count('TGA') == 0):
        ## if the given codon respect all the above conditions, keep it
        ## only keep the codons that have less missing sequence than the threshold
        ## consider a codon as missing if all its nucleotide are missing
        ## should I also add the case when only some nucleotides are missing?
        ## also need to remove the codons when there is a STOP codon in at least one strain
        ## (alternative: could replace the STOP codons with gaps and then treat it like any other missing codon ie remove the codon only if more than 1% of the strains are missing this codon)
        ## STOP codon = TAG, TAA, TGA 
        ## also, save the true position of the codon     
        posCDScore.append([posCDS[i:(i+3)]])
        for j in range(0,nbSeq):
            CDSfilter[j].append(seq[j])



posCDScore = [k for i in posCDScore for j in i for k in j]


## rewrite the new alignment in a different file
with open("coreCDS.aln", 'a') as f:
    for i in list(range(0,nbSeq)):
        f.write(strain[i] + "\n")
        [f.write(j) for j in CDSfilter[i]]
        f.write("\n")

## also save the sequences
## in one file per strain
## rename the first strain of the dataset as reference (same as 26695 but avoid double)
strain[0] = "Reference"
for i in list(range(0, nbSeq)):
    with open("strains/" + strain[i][1:] + ".fasta", 'w') as f:
        f.write(strain[i] + "\n")
        [f.write(j) for j in CDSfilter[i]]
        f.write("\n")

## as well as the true position of the different core CDS kept
## ordered following the new alignment file
## ATTENTION: indexation in python (posCDS is based on index and first index = 0 -> there is a decalage on 1 between the position given here and the position in the gff file)
## !! need to take that into account when look at the alleles for the signifcant sites in the GWAS 
with open("posCoreCDS.txt", 'w') as f:
    [f.write(str(i + 1) + '\n') for i in posCDScore]

