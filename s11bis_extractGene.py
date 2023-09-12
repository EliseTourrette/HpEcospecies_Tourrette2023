#!/usr/bin/python3

## want to extract the nucleotide sequence of a given gene for the different strains from the annotation files (ffn or fna)
## note that the gene can be changed by changing the gene in grep (to get the tag of the gene)
## could also put the gene name as an argument
## use the product of the gene as an argument, however one product can correspond to multiple genes (in particular with grep)
## thus, need, for each strain, the list of the gene(s) that we actually want (with the gene product, the gene name and the gene tag specific to the strains)

from Bio import SeqIO
import subprocess
import sys

nbArg = len(sys.argv)

strain = sys.argv[1]
gene = sys.argv[2]
seqType = sys.argv[3]
saveFolder = sys.argv[4]
thr = int(float(sys.argv[5])) ## threshold in the number of site/aa above which we consider a gene fragment as more or less complete
## if the length threshold is a float string, first need to convert it as a float (int either expect a int string or a float)
FILE = strain + "/" + strain + "." + seqType

## in some cases, want to give a list of gene as an argument (the 6th argument)
## in particular, when there are multiple genes with the same product
## need to do a distinction with the tag
## note that listGene is a file that contains the tag of the wanted gene(s) in the first column
if nbArg > 6:
    listGeneFile = sys.argv[6]
    ## the files in listGene contains, for each strains the list of the significant gene with their product, gene name and gene tag (specific for the strain)
    listGene = [] 
    with open(listGeneFile, 'r') as f:
        for i in f:
            iline = i.split("\n")[0].split("\t")
            listGene.append(iline[0])
        
    ## get the id of the VacA gene/fragments for a given strain
    tag = subprocess.run('grep "' + gene + '" ' + FILE, shell = True, capture_output = True, check = True, text = True)
    tag = tag.stdout.split("\n")
    tag = [i.split(" ")[0][1:] for i in tag if len(i) > 0]    
    ## check which tag actually correspond to the wanted gene
    indexTag = [i.split("_")[-1] for i in tag if listGene.count(i) > 0]
    tag = [i for i in tag if listGene.count(i) > 0]
else:
    indexTag = [i.split("_")[-1] for i in tag]


genomeSubset = []
genome = list(SeqIO.parse(FILE, "fasta"))
genomeID = [record.id for record in genome]
for igene in tag:
    ix = genomeID.index(igene)
    genomeSubset.append(genome[ix])

genomeSubset = [list(i.seq) for i in genomeSubset]

## concatenate together the different fragment 
## IF their individual length is less than a given threshold (fraction of the total length of the reference version of the gene)
## IF their locus tag are adjacent
## IF the length of the previous sequence is also less than the threshold length
tmp = []
ix = 0
for i in range(0, len(indexTag)):
    if i == 0:
        tmp.append(genomeSubset[i])
    else:
        if (len(genomeSubset[i]) < thr) & (int(indexTag[i]) - int(indexTag[i-1]) == 1) & (len(tmp[ix]) < thr):
            tmp[ix] = tmp[ix] + genomeSubset[i]
        else:
            tmp.append(genomeSubset[i])
            ix += 1

geneName = gene.split(" ")[-1].split("/") ## if have a '/' in the name file, it will be a problem when want to read it via bash
with open(saveFolder + strain + "_" + geneName[-1] + "." + seqType, "w") as f:
    for i in range(0, len(tmp)):
        f.write(">" + strain + "." + str(i+1) + "\n")
        f.write("".join(tmp[i]) + "\n")

