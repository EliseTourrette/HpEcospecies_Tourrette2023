#!/bin/bash

## instead of saying whether we are calculating the within or between pop dN/dS
## give, as an argument, the name of the file containing the name of the strains for the comparisons to be done
## separate the file containing the comparisons into different section (PAML is run sequentially for each section while the different section are run in parallel)
## arg1 = where (line) the section begin, arg2 = the number of lines in the section, arg3 = file that contains the different comparisons, arg4 = folder where the fasta files are put
## need the path where the executable of yn00 (PAML) is saved

echo $1
echo $3
echo $4
listStrain=$(seq $1 $(($1 + $2)))
for i in $listStrain; do
    echo $i
    comp=$(head -n $i $3 | tail -n 1)
    strain1=$(echo $comp | cut -d " " -f1)
    strain2=$(echo $comp | cut -d " " -f2)
    res="res_$strain1-$strain2"
    folder="out_$strain1-$strain2"
    mkdir $folder
    cat $4/$strain1.fasta $4/$strain2.fasta > $folder/ali.fasta
    cp /Users/elise/Desktop/HpGlobal/PAML/yn00.ctl $folder/
    cd $folder 
    echo $res
    /Users/elise/Desktop/HpGP1k-Genome/ANALYSIS/paml_hpgp/paml4.9j/bin/yn00 yn00.ctl > out
    val=$(grep -A 9 '(A) Nei-Gojobori' res | tail -n 1)
    echo "$strain1 $val" >> ../res/NG86_$1
    val=$(grep -A 8 '(B) Yang & Nielsen' res | tail -n 1 | cut -d " " -f 9-)
    echo "$strain1 $strain2 $val" >> ../res/YN00_$1
    val=$(head -n 3 res | tail -n 1)
    echo "$strain1 $strain2 $val" >> ../res/nbCodons_$1
    cd ..
    rm -r $folder
done
