#!/bin/bash

## for each gene that returned a hit in H. cetorum
## extract its sequence for a subset of strain (strains for which I have the annotation)
## align them together using mafft
## and generate a ML tree using FastTree

mkdir ntSeq
mkdir ffnFile
mkdir alnFile
mkdir nwkFile
mkdir outFile

geneNb=0
while IFS= read -r line; do
    geneNb=$((geneNb + 1))
    geneName=$(echo $line | cut -d "," -f1)
    thrLength=$(echo $line | rev | cut -d "," -f1 | rev)
    echo "GENE $geneNb =====> $geneName"
    cd Annotation
    for strain in *; do
        inData=$(grep $strain ../listStrains.txt)
        inData=${#inData}
        if [[ $inData -gt 0 ]]; then
            echo "strain $geneNb ---------------------> $strain"
            python3 ../s11bis_extractGene.py "$strain" "$geneName" "ffn" "../ntSeq/" "$thrLength" "../listGene/$strain.txt"
        fi
    done
    cd ..
    cat ntSeq/* >> ffnFile/"$geneNb"_nt.ffn
    rm ntSeq/*
    
    seqProt=$(grep "$geneName" HcetorumSequence/seqHitsHcetorum.txt | rev | cut -d " " -f1 | rev)
    if [ ${#seqProt} -ge 1 ];then
        echo ">Hcetorum" >> ffnFile/"$geneNb"_nt.ffn
        echo "$seqProt" >> ffnFile/"$geneNb"_nt.ffn
    fi
    
    mafft --auto ffnFile/"$geneNb"_nt.ffn > alnFile/"$geneNb"_nt.aln 2> outFile/out_"$geneNb"nt.out
    FastTree -nt alnFile/"$geneNb"_nt.aln > nwkFile/"$geneNb"_nt.nwk 2> outFile/tree_"$geneNb"nt.error
done < listHitsHcetorum.csv

rm -r ntSeq
