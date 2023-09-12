## file to generate the input files for the bugwas analysis
## in particular, the genotype file
## for which I replace 0/1/2/3 by A/T/C/G

## also provide the option to remove the variants that are not biallelic (if biallelic = True)

biallelic = True
header = []
pos = []
ref = []
alt = []
snp = []
i = 0
with open("tmp2.txt", 'r') as f:
    for l in f:
        iline = l.rstrip('\n').split()
        i += 1
        print(i)
        if i == 1:
            header.append(iline[0])
            header.append(iline[3:])
        if i > 1:
            if not biallelic:
                pos.append(iline[0])
                ref = iline[1]
                alt = iline[2].split(",")
                seq = iline[3:]
                seq = [ref if seq[j] == '0' else alt[int(seq[j]) - 1] for j in range(0,len(seq))]
                snp.append(seq)
            if biallelic:
                ref = iline[1]
                alt = iline[2].split(",")
                if len(alt) == 1:
                    pos.append(iline[0])
                    seq = iline[3:]
                    seq = [ref if seq[j] == '0' else alt[int(seq[j]) - 1] for j in range(0,len(seq))]
                    snp.append(seq)


if biallelic:
    nameFile = "geno_biallelic_SNP.txt"
else:
    nameFile = "geno_multiallelic_SNP.txt"

with open(nameFile, 'w') as f:
    f.write(header[0] + '\t' + '\t'.join(header[1]) + '\n')
    for i in list(range(0,len(pos))):
        f.write(pos[i] + '\t' + '\t'.join(snp[i]) + '\n')

