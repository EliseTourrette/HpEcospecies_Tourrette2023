## python3 script 

## in a vcf file generated from coreCDS*.aln
## replace the position of the variants by their true position
## the positions in the position are the position in the coreCDS alignment
## which only contains the core CDS
## hence the position in the generated vcf file are not the true positions
## use the file posCoreCDS.txt
## that contains the true positions (ie in the original alignment file)
## the positions are in the same order as in coreCDS*.aln

## in GWAS/input

pos = []
with open("../PREPROCESSING/posCoreCDS.txt", 'r') as f:
    for l in f:
        iline = l.rstrip('\n').split()
        pos.append(iline[0])

vcf = []
ix = 0
## should give the path to the vcf file for which we want to change the variant positions
## as an argument
with open("tmp1.vcf", 'r') as f:
    for l in f:
        ix += 1
        iline = l.rstrip('\n').split()
        if ix > 4:
             ## index begins at 0 while the position begins at 1
            iline[1] = pos[int(iline[1]) - 1]
        vcf.append(iline)

with open("tmp2.vcf", 'w') as f:
    for i in range(0, len(vcf)):
        f.write('\t'.join(vcf[i]) + '\n')  
