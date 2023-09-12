## extract the wanted strains from the global alignment file
## for a specific subpop

## in HpGlobal_noHpGP/GWAS/input
 
# read the list of wanted strains
strainsubpop = []
with open("strainsSubpop.txt", 'r') as f:
    for l in f:
        iline = l.rstrip('\n').split()
        strainsubpop.append(iline[0])


strain = []
sequence = []
saveSeq = "No"
with open("../../DATA/woHpGP_dnds.cvg70.aln", 'r') as f:
    for l in f:
        iline = l.rstrip('\n').split()
        if(iline[0][0] == '>') & (strainsubpop.count(iline[0][1:]) > 0):
            strain.append(iline[0][1:])
            saveSeq = "Yes"
        else:
            if saveSeq == "Yes":
                sequence.append(iline[0])
                saveSeq = "No"
            else:
                saveSeq = "No"


with open("alignmentSub.aln", 'a') as f:
    for i in list(range(0,len(strain))):
        f.write(">" + strain[i] + "\n")
        [f.write(j) for j in sequence[i]]
        f.write("\n")
