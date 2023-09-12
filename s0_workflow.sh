## Analysis of the global H. pylori dataset
## some strains of theoriginal dataset have been removed due to an embargo on them (not published) 
## the results of the analysis are mostly the same even without them
## exception: hspIndigenousSAmerica is not a "focus" population anymore - only hspIndigenousNAmerica and hspSiberia

## folder in which all the data and results are saved
cd /Users/elise/Desktop/HpGlobal_noHpGP

## folder that contains the data
mkdir DATA
## folder for the plots
mkdir PLOT

## clean up the names of the strains
## initial data: strains aligned to the 26695 reference genome -> global alignment
cd DATA
gzip -d woHpGP_dnds.cvg70.aln.gz
sed 's/.fasta.delta.3//' woHpGP_dnds.cvg70.aln > woHpGP_dnds.cvg70_2.aln
rm woHpGP_dnds.cvg70.aln
mv woHpGP_dnds.cvg70_2.aln woHpGP_dnds.cvg70.aln
cd ..

## folder for some information about the data
mkdir METADATA

## list the strains that are included in the dataset
grep ">" DATA/woHpGP_dnds.cvg70.aln | cut -b2- > METADATA/listStrains.txt

## when make change on the data
## but not the full analysis
mkdir PREPROCESSING
cd PREPROCESSING
mkdir strains
mkdir CDS

## get the list of CDS
grep "26695\tProdigal:2.6\tCDS" ../DATA/26695.gff > posCDS.csv

## prepare the sequences for the dN/dS analysis:
## remove the STOP codons and the codons that are completely absent (---) in more than 1% of the sequences
## input: DATA/woHpGP_dnds.cvg70.aln [global alignment], posCDS.csv [position of the CDS]
## output: strains/strain_i.aln [one file per strain - all CDS in each file -> cleaned up strain sequence], CDS/CDS_i.aln [one file per CDS - all strain in each file -> cleaned up CDS sequence], coreCDS.aln [all CDS, all strain], posCoreCDS.txt [position of the kept CDS - in practice no CDS have been removed]
python3 s1_preprocess.py

grep ">" coreCDS.aln | sort | uniq -d
mv coreCDS.aln coreCDSclean.aln

#######
### PCA
#######
cd ..
mkdir PCA
mkdir PLOT/PCA
cd PCA

## did the PCA on different set of SNPs
## to check whether the different set of SNPs change the results 
## (not the case)
mkdir coreAll
mkdir fullAll
mkdir cds

## here: do the PCA on the core (100%) SNPS + on all the SNPs, both for CDS and intergenic sequences
## use the variants that were previously called by Roberto
## compare the PCA results to what is obtained when doing the same things as for the whole (HpGP included) dataset
#snp-sites -v -o cds/res.vcf ../PREPROCESSING/coreCDSclean.aln > out
#snp-sites -c -v -o coreAll/res.vcf ../DATA/woHpGP_dnds.cvg70.aln
snp-sites -v -o fullAll/res.vcf ../DATA/woHpGP_dnds.cvg70.aln  
rm out
# 797586 SNPs

## PCA done with plink
## command line with all the SNPs
## same command line for the other set of SNPs
cd fullAll
## first: LD pruning in order to remove linked SNPs
plink --vcf res.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 10 0.1 --out hpglobal > plink.prune.out
## then: PCA + SNPs loadings
plink --vcf res.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# --extract hpglobal.prune.in --make-bed --pca var-wts --out hpglobal_LD_PCA >> plink.out
## exit cluster

# processing of the results
## input: PCA/fullAll/hpglobal_LD_PCA.eigenvec [eigenvector], PCA/fullAll/hpglobal_LD_PCA.eigenval [eigenvalue], METADATA/strains_pop.csv [pop assignment for each strain]
## output: PLOT/PCA/PCAcoreCDS_africa2.jpg [plot PC1 vs PC2], PLOT/PCA/PCA2coreCDS_africa2.jpg [plot PC3 vs PC4]
## here did the plots with all the SNPs, but can change the name of the files (folder) in the script
Rscript s2_analysis_PCA.r

#######
### GWAS
#######
cd ..
mkdir GWAS
mkdir PLOT/GWAS
mkdir PLOT/GWAS/regions
cd GWAS
mkdir input
cd input

## use the bugwas R package to do the GWAS
## generation of the necessary input files: phenotype file, genotype file and tree of the dataset
## only use the strains from the two focus populations (ie the populations with the two ecospecies- as seen in the PCA plot)

## extract the strains from hspSiberia and hspIndigenousNAmerica ("focus" populations)
grep 'hspSiberia\|hspIndigenousNAmerica' ../../METADATA/strains_pop.csv | cut -d "," -f1  > strainsSubpop.txt

## extract the sequences of the strains from the focus populations from the global alignment
## input: strainsSubpop.txt [list of the strains of the focus populations], DATA/woHpGP_dnds.cvg70.aln [global alignment]
## output: alignmentSub.aln [global alignment but with only the strains from the focus populations]
python3 s3_extractStrainsSubpop.py

## concatenate the cleaned-up sequence for all the strains from the two focus populations
## based on the list of strains previously generated
#while IFS= read -r line; do
#    file="../../PREPROCESSING/strains/$line.fasta"
#    if test -f "$file"; then
#        cat ../../PREPROCESSING/strains/$line.fasta >> subpopSiberiaIndNAm_coreCDS.aln
#    else
#        echo $line
#    fi
#done < strainsSubpop.txt

## generate the phylogenetic tree for the gwas needed as an input for the gwas
## did it with the global alignment and the cleaned-up CDS
nohup FastTree -nt alignmentSub.aln > treeFull.nwk 2> treeFull.error &
#nohup FastTree -nt subpopSiberiaIndNAm_coreCDS.aln > treeCDS.nwk 2> treeCDS.error &
## also generate the tree for the whole dataset (saved in the folder TREE)
nohup FastTree -nt woHpGP_dnds.cvg70.aln > wholePop.nwk 2> wholePop.error &

## generate the genotype file (use all the SNPs - both CDS and intergenic)
snp-sites -c -v -o tmp1.vcf alignmentSub.aln
# 278019 SNPs: used the full sequence (CDS and intergenic regions), BIALLELIC SNPs
sed "s/#CHROM/CHROM/" tmp1.vcf | grep -v "#" | cut -f 2,4-5,10-219  | sed "1s/POS/ps/" > tmp2.txt
## input: tmp2.txt [temporary file that contains the alleles for each SNP]
## output: geno_biallelic_SNP.txt [genotype file, use only BIALLELIC SNPs]
python3 s4_input_bugwas.py
rm tmp*
## used the tree done with all (CDS and intergenic) sites ie treeFull
cd ..
## generation of the input + GWAS +  postprocessing + plot of the trees
## input: 
##          1. generation input: METADATA/strains_SiberiaIndAm.csv [only contains the strains from the focus pop, extract the phenotype; this file was generated manually based on the metadata of t he different strains - such as their location, ecospecies appartenance etc] -> input/pheno.txt [phenotype file], strains_pop.txt [pop information of all the strains], input/tree [tree of the data used for the GWAS]
##          2. GWAS: input/geno_biallelic_SNP.txt [genotype file], input/pheno.txt, input/tree, path to the gemma exe
##          3. analyses: DATA/26695.gff [annotation file of the reference genome - downloaded from the NCBI], bugwas_biallelic_lmmout_allSNPs.txt [bugwas output]
## output: 
##          1. generation input: input/pheno.txt [phenotype file], PLOT/TREE/TREE_GWAS.jpeg [plot of the tree used as input fro the GWAS]
##          2. GWAS: bugwas output + data_GWAS.RData [GWAS output but as a RData object]
##          3. analyses: PLOT/GWAS/regions/* [Manhattan plot with a zoom on the different regions of the genome - ie zoom], significantCDS.csv[list of significant CDS], significantSNPs.txt [list of significant SNPs]
Rscript s5_analysis_GWAS.r

#######
### FST
#######
cd ..
mkdir FST 
cd FST
mkdir fasta_core
## use the same alignment as for the GWAS
## want to look at the same SNPs
mv ../GWAS/input/alignmentSub.aln  fasta_core/

## r script to calculate the FST between high and low strains
## input: METADATA/strains_SiberiaIndAm.csv, fasta_core [folder tha contains the alignment file that will be sued as input]
## output: FST_highLow.RData [RData object with the FST values for each SNPs]
Rscript s6_analysis_FST.r


#######
### separate the genes into differentiated and undifferentiated, based on the GWAS and FST results
#######

### haplotype visualization at the differentiated genes
cd ..
mkdir HAPLOTYPES

## determine the biallelic SNPs that are significant for the GWAS and the FST
## input: 
##              1. METADATA/strains_pop.csv, METADATA/strains_SiberiaIndAm.csv, DATA/26695.gff [annotation of 26695], FST/FST_highLow.RData, GWAS/bugwas_biallelic_lmmout_allSNPs.txt
##              2. HAPLOTYPES/sites_pval10FST0.9.txt [haplotypes at the differentiated SNPs, for all strains, extracted from the global alignment]
##              3. HAPLOTYPES/nbdiff_lowLoadStrains.txt [genetic distance to the low load strains from the focus populations]
## output: 
##              1. summary: GWAS/ FST: HAPLOTYPES/FST_GWASpval_lowHigh.txt [file that summarize the GWAS p-value and FST result for each SNP], HAPLOTYPES/pos_pval10FST0.9.txt [positions of the differentiated SNPs], HAPLOTYPES/pos_nonSigSites_pval10FST0.5.txt [positions of the undifferentiated SNPs], HAPLOTYPES/pos_otherSNPs.txt [positions of the intermediate SNPs]
##              2. visualization: HAPLOTYPES/data_pval10FST0.9.RData [haplotype at the differentiated SNPs, but also include additional information concerning the strains]
##              3. count of the high alleles: HAPLOTYPES/countHighAllele_pval10FST0.9.RData [number of "high" allele at the differentiated SNPs, for all strains]
##              4. blocks: dataBlock_pval10FST0.9.RData / highBlock_pval10FST0.9.RData [blocks of "high" alleles"]
##              5. number of high alleles vs genetic distance: PLOT/HAPLOTYPES/nbHighAlleleVSdistLow.jpeg
Rscript s7_analysis_GWAS-FST.r

## need to be called in the middle of  s7_analysis_GWAS-FST.r
## generate the haplotypes at the significant sites
## as need the position of the differentiated SNPs as input
## input: HAPLOTYPES/pos_pval10FST0.9.txt, DATA/woHpGP_dnds.cvg70.aln
## output: HAPLOTYPES/sites_pval10FST0.9.txt [haplotypes at the differentiated SNPs, for all strains, extracted from the global alignment]
python3 s7bis_haplotypeSigSitesDivergence.py

## need to be called in the middle of  s7_analysis_GWAS-FST.r
## although could technically be called beforehand as do not need any data generated by s7_analysis_GWAS-FST.r (as for th ehaplotype generation)
## genetic distance with the low load strains 
## pairwise comparison
## this time, use the whole genome to calculate the number of differences between the two strains and not only the CDS
## input: METADATA/listStrains.txt, METADATA/strains_SiberiaIndAm.csv, DATA/woHpGP_dnds.cvg70.aln
## output: HAPLOTYPES/nbDiff_straini.txt [number of differences between the straini and all the other strains from the dataset]
for i in {0..6874..700}; do
    python3 s7ter_nbDiffSeq.py $i &
done
cat HAPLOTYPES/nbDiff_* >> HAPLOTYPES/nbdiff_lowLoadStrains.txt
rm HAPLOTYPES/nbDiff_*

#######
### dN/dS & tree - separate differentiated and undifferentiated genes
#######

## generate the alignment files for the differentiated, undifferentiated and intermediate CDS
## based on the number of (un)differentiated SNPs they have
mkdir DIVERGENT_NO_CDS
mkdir DIVERGENT_NO_CDS/highGenesSeq
mkdir DIVERGENT_NO_CDS/highDivStrains
mkdir DIVERGENT_NO_CDS/lowDivStrains
mkdir DIVERGENT_NO_CDS/interDivStrains
## input: HAPLOTYPES/pos_pval10FST0.9.txt, HAPLOTYPES/pos_nonSigSites_pval10FST0.5.txt, HAPLOTYPES/pos_otherSNPs.txt, DATA/woHpGP_dnds.cvg70.aln, PREPROCESSING/posCDS.csv
## output: highCoreCDS.aln [differentiated genes sequence for all strains], lowCoreCDS.aln [undifferentiated genes sequence for all strains], interCoreCDS.aln [intermediate genes sequence for all strains], highGenesSeq [folder with one file per differentiated gene], highDivStrains [folder with one file per strain, each containing all differentiated genes], lowDivStrains [folder with one file per strain, each containing all undifferentiated genes], interDivStrains [folder with one file per strain, each containing all intermediate genes], geneHigh.txt [list of differentiated genes], , geneLow.txt [list of undifferentiated genes], geneInter.txt [list of intermediate genes], nbSigSitesHigh.txt [number of differentiated SNPs for each differentiated gene]
python3 s8_separationDiffUndiffGenes.py

nohup FastTree -nt highCoreCDS.aln > high.nwk 2> high.error &
nohup FastTree -nt lowCoreCDS.aln > low.nwk 2> low.error &
nohup FastTree -nt interCoreCDS.aln > inter.nwk 2> inter.error &


## dN/dS for the three sets of genes
## generate the list of comparisons to do
## in our case: outgroup = hpAfrica2
## save both the NG86 and YN0 values (but remove the output of PAML)
cd DIVERGENT_NO_CDS
mkdir res
echo "strain1 strain2   dN/dS (dN,dS)" > res/NG86
echo "strain1 strain2   S       N        t   kappa   omega     dN +- SE    dS +- SE" > res/YN00

## need to do it separately for differentiated (high), undifferentiated (low) and intermediate genes
## cannot directly run this part in the terminal
## need to do it for the three possible values of folderComparisons
folderComparisons="lowDivStrains"
folderComparisons="highDivStrains"
folderComparisons="interDivStrains"

## creation of the file that list the pairwise comparisons to be done
## input: METADATA/listStrains.txt, METADATA/strains_pop.csv, METADATA/strains_pangenome.txt [subset of strains that were used in some other analysis for which couldn't use the whole dataset]
## output: METADATA/pairwiseComparisons_outgroup.txt [between population comparisons; here outgroup = hpAfrica2], METADATA/pairwiseComparisons_HardyUbiquitous.txt [comparison between Hardy and Ubiquitous strains]
Rscript s9_ createFile_dNdScomparisons.r
fileComparisons="../METADATA/pairwiseComparisons_outgroup.txt"
wc -l $fileComparisons

## pairwise calculation of the dN/dS and dS
## define the comparisons to be done beforehands in fileComparisons
## not that to be faster, do a for loop over all the comparisons to be done (from 1 to the number of line of fileComparisons, obtained with wc -l); each are independant from the other
## input: fileComparison [list of pairwise comparisons to be done], folderComparison [folder where the strain sequences are saved], XX.fasta [sequence of the strain XX for a given set of genes], yn00.ctl [parameter file of YN00], path to yn00
## output: output of PAML + results of NG86 and YN00 estimation saved in two different files
for i in {1..213416..10000}; do
    bash s10_runPAML.sh $i 10000 $fileComparisons $folderComparisons > "$i".out 2> "$i"_error.out &
done

mkdir dNdS_between_low
mkdir dNdS_between_high
mkdir dNdS_between_inter

## as for folderComparison, need to change the name of the folder where the results will be saved (dNdS_between_XX/)
cat res/NG86* | uniq >> dNdS_between_inter/NG86
cat res/YN00* | uniq >> dNdS_between_inter/YN00

## remove the output from PAML (could still saved them if wanted to look at the different results in more details,
## but we directly extracted what was interesting for us)
rm -r res; rm *.out


#######
### comparison of the differentiated genes sequence with H. cetorum sequence
## H. cetorum: species with the higher number of hits for the different genes compared to the other Helicobacter species
#######

## BLAST the differnetiated gene sequence against the H. cetorum genome
## use one sequence for each ecospecies X focus population (+ Hacinonychis and primate) version
## ie Hardy - hspSiberia / Ubiquitous hspSiberia / Hardu hspIndigenousNAmerica ...

### differentiated genes trees with H. cetorum sequence
## !!!!! for the genes that returned a hit in H. cetorum !!!!
mkdir Hcetorum
cd Hcetorum
## get the list of differentiated genes
## ignore the hypothetical proteins without unique identifier
while IFS= read -r line; do
    product=$(echo $line | cut -d "=" -f4 | cut -d ";" -f 1)
    echo $product >> tmp
done < ../listGenes/geneHigh.txt
grep -v "^26695_" tmp  | cut -d '_' -f1 > listGene.txt

## need to check that the strain is in the dataset
## extract the info about the differentiated genes for all the strains (for which there is an annotation available)
mkdir listGene
while IFS= read -r line; do
    product=$(echo $line | cut -d "=" -f7)
    npr=${#product}
    if [[ $npr -gt 0 ]]; then
        echo $product
        for istr in ../../../HpGlobal/DATA/Annotation_Roberto/*; do
            strain=$(echo $istr | cut -d '/' -f7)
            inData=$(grep $strain ../../METADATA/listStrains.txt)
            inData=${#inData}
            if [[ $inData -gt 0 ]]; then
                grep "$product" $istr/$strain.tsv | grep "CDS" >> listGene/$strain.txt
            fi
        done
    fi
done < listGene.txt
rm tmp

cat listGene/* >> listAnn.txt
rm -r listGene

## listHitsHcetorum comes from the HpGP analysis
## I just added the last gene (length sequence is from the 26695 version)
## and removed the genes that were not differentiated without the HpGP strains

mkdir HcetorumSequence
cd HcetorumSequence
mkdir blast_results

## use the same seqGenes as the one used for the dataset with the HpGP strains 
## remove the HpGP-CHI strain (except for this one, all the other are in the dataset without HpGP)
## manually add the sequence of the two additional genes and remove that were not detected wo the HpGP
cd seqGenes
for strain in *; do
    strain=$(echo "$strain" | cut -d '.' -f1)
    echo $strain
    blastn -query "$strain".ffn -db ../blast_db/H_cetorum_MIT99-5656_enterobase -outfmt "6 qseqid bitscore sseq" >> ../blast_results/"$strain".out 2>> ../blast_results/"$strain".error
done

## for each hot found the the H.cetorum genome
## determine which version give the best hit
## ie is the H. cetorum version of the genome more like a Hardy or an Ubiquitous version?
cd ..
grep ">" seqGenes/* > listGenes.txt
while IFS= read -r line; do
    gene=$(echo $line | cut -d"," -f1)
    echo "$gene"
    for strains in blast_results/*.out; do
        strShort=$(echo "$strains" | cut -d"/" -f2 | cut -d"." -f1)
        index=$(grep "$gene" listGenes.txt | grep "$strShort" | cut -d":" -f2 | cut -d" " -f1 | cut -b2-)
        echo $index
        if [ ${#index} -gt 0 ]; then
            grep "$index"  "$strains" >> tmp
        fi
    done
    seqBest=$(sort -nk2 tmp | tail -n 1)
    strainBest=$(echo $seqBest | cut -d" " -f1)
    newStrain=$(grep "$strainBest" listGenes.txt | cut -d":" -f2 | cut -b2-)
    echo $seqBest | sed "s:$strainBest:$newStrain:" >> seqHitsHcetorum.txt
    rm tmp
done < ../listHitsHcetorum.csv

cd ..
## for each gene that returned a hit in H. cetorum
## extract its sequence for a subset of strain (and for which have the annotation)
## align them together using mafft
## and generate a ML tree using FastTree
## use the python script s11bis_extractGene.py inside (may need to adjust the path to the script depending on where it is)
## for s11bis_extractGene.py
## input: strain, gene name, type of sequence, folder to saveresu;ts, threshold of the length of the sequence above which the sequence is considered to a new version (and not a fragment of the gene) - based on the reference gene length, list of genes for this strain
## output: sequence for this strain (with different fragments concatenated)
## for s11_treeHitsHcetorum.sh
## input: listHitsHcetorum.csv [list of the genes that returned a hit when compared to H. cetorum]
## output: alnFile/"$geneNb"_nt.aln [alignment file for each gene], nwkFile/"$geneNb"_nt.nwk [tree file for each gene]
nohup bash s11_treeHitsHcetorum.sh > outTreeHits.out &

cat blast_results/*.out >> blast_results.out

#######
### focus on UreA and UreB
#######

## use the script s11bis_extractGene.py to extract UreA and UreB genes for all ythe strains
cd ..
mkdir UreAB
cd UreAB

## change the focus gene manually
## ie run the same code by change the values between VacA, UreA and UreB

focusGene="VacA"
geneTag="HP0887"
geneName="Vacuolating cytotoxin VacA"
## the sequence of the 26695 VacA is stored in the file 26695_VacA.fna

focusGene="UreA"
geneTag="HP0073"
geneName="Urease subunit alpha UreA"
## the sequence of the 26695 UreA is stored in the file 26695_UreA.fna

focusGene="UreB"
geneTag="HP0072"
geneName="Urease subunit beta UreB"
## the sequence of the 26695 UreA is stored in the file 26695_UreB.fna

cd $focusGene
mkdir Nucleotide
mkdir Nucleotide/seq_annotation

thrLength="3500" ## for VacA
thrLength="645" ## for UreA
thrLength="1540" ## for UreB
## extract the sequence for the reference
python3 s11bis_extractGene.py "26695" "$geneName" "ffn" "../../$focusGene/Nucleotide/seq_annotation/" "$thrLength"
## extract the sequences for all strains for which I have the annotation
for strain in *; do
    setStrain=$(echo "$strain" | cut -b1-4)
    if [ "$setStrain" = "HpGP" ]; then
        echo $strain
        python3 s11bis_extractGene.py "$strain" "$geneName" "ffn" "../../$focusGene/Nucleotide/seq_annotation/" "$thrLength"
    fi
done
cd ../../"$focusGene"/Nucleotide

## concatenate all the sequences into one file
cat seq_annotation/*.ffn >> "$focusGene".ffn

## separation of UreA and UreB into their two versions
## strains1.txt and strains2.txt was generated manually
while IFS= read -r line; do
    grep -A 1 $line UreA.ffn >> UreA_nt1.ffn
done < strains1.txt

while IFS= read -r line; do
    grep -A 1 $line UreA.ffn >> UreA_nt2.ffn
done < strains2.txt

## strainsB1.txt and strainsB2.txt was generated manually
while IFS= read -r line; do
    grep -A 1 $line UreB.ffn >> UreB_nt1.ffn
done < strainsB1.txt

while IFS= read -r line; do
    grep -A 1 $line UreB.ffn >> UreB_nt2.ffn
done < strainsB2.txt

## alignment of the different sequences
mafft --auto UreA.ffn > UreA.aln 2> outA.out &
mafft --auto UreA_nt1.ffn > UreA_nt1.aln 2> outA1.out &
mafft --auto UreA_nt2.ffn > UreA_nt2.aln 2> outA2.out &
mafft --auto UreB.ffn > UreB.aln 2> outB.out &
mafft --auto UreB_nt1.ffn > UreB_nt1.aln 2> outB1.out &
mafft --auto UreB_nt2.ffn > UreB_nt2.aln 2> outB2.out &

## generation of the phylogenetic trees
nohup FastTree -nt UreA.aln > UreA.nwk 2> treeA.error &
nohup FastTree -nt UreA_nt1.aln > UreA_nt1.nwk 2> treeA1.error &
nohup FastTree -nt UreA_nt2.aln > UreA_nt2.nwk 2> treeA2.error &
nohup FastTree -nt UreB.aln > UreB.nwk 2> treeB.error &
nohup FastTree -nt UreB_nt1.aln > UreB_nt1.nwk 2> treeB1.error &
nohup FastTree -nt UreB_nt2.aln > UreB_nt2.nwk 2> treeB2.error &


#######
### dN/dS between the different strains and the monkey strains
## check if lower dS between monkey strains and hpSahul
#######

# use the whole genome sequence without the STOP codons
# and the codons that are --- in more than 1% of the strains
cd ../..
mkdir PAML
cd PAML

# between the two Sahul monkey strains (high divergence ones)
# I could concatenate, for each strain the high / low / intermediate CDS to create a global CDS sequence...
strain1="HEL_BA9303AA_AS"
strain2="HEL_BA9707AA_AS"
folder="out_$strain1-$strain2"
mkdir $folder
cat ../PREPROCESSING/strains/$strain1.fasta ../PREPROCESSING/strains/$strain2.fasta > $folder/ali.fasta
cp /Users/elise/Desktop/HpGlobal/PAML/yn00.ctl $folder/
cd $folder 
/Users/elise/Desktop/HpGP1k-Genome/ANALYSIS/paml_hpgp/paml4.9j/bin/yn00 yn00.ctl > out
rm ali.fasta
cd ..

mkdir res
echo "strain1 strain2   dN/dS (dN,dS)" > res/NG86
echo "strain1 strain2   S       N        t   kappa   omega     dN +- SE    dS +- SE" > res/YN00
echo "strain1 strain2  nbCodons" > res/nbCodons

folderComparisons="../PREPROCESSING/strains"
fileComparisons="pairwiseComparisons.txt"
nbComp=$(wc -l $fileComparisons | cut -d " " -f4)

for i in {1..11606..10000}; do
    bash s10_runPAML.sh $i 10000 $fileComparisons $folderComparisons > "$i".out 2> "$i"_error.out &
done

mkdir monkeySahul
cat res/NG86* | uniq >> monkeySahul/NG86
cat res/YN00* | uniq >> monkeySahul/YN00
cat res/nbCodons* | uniq >> monkeySahul/nbCodons
rm -r res; rm *.out

#######
### parallel evolutionary history between hardy and ubiquitous ecotypes
#######

## first extract the strains from the populations with both ecotypes (Siberia + IndigenousNAmerica) + Africa2 (for a comparison to H. acinonychis) + hpSahul (for the primate strains)
mkdir TREE/evolutionaryHistory
cd TREE/evolutionaryHistory

## separate the differentiated genes into one file for the Hardy strains and one file for the Ubiquitous strains
## and one global file for the undifferentiated genes
while IFS= read -r line; do
    ecotype=$(echo $line | cut -d"," -f3)
    if [ $ecotype = "hardy" ]; then
        nameStrain=$(echo $line | cut -d"," -f1)
        cat ../../DIVERGENT_NO_CDS/alignment/highDivStrains/"$nameStrain".fasta >> hardyDiff.aln
        nameStrain=$(echo $line | cut -d"," -f2)
        cat ../../DIVERGENT_NO_CDS/alignment/highDivStrains/"$nameStrain".fasta >> hardyDiff.aln
    else
        nameStrain=$(echo $line | cut -d"," -f1)
        cat ../../DIVERGENT_NO_CDS/alignment/highDivStrains/"$nameStrain".fasta >> ubiDiff.aln
        nameStrain=$(echo $line | cut -d"," -f2)
        cat ../../DIVERGENT_NO_CDS/alignment/highDivStrains/"$nameStrain".fasta >> ubiDiff.aln
    fi
    nameStrain=$(echo $line | cut -d"," -f1)
    cat ../../DIVERGENT_NO_CDS/alignment/lowDivStrains/"$nameStrain".fasta >> undiff.aln
    nameStrain=$(echo $line | cut -d"," -f2)
    cat ../../DIVERGENT_NO_CDS/alignment/lowDivStrains/"$nameStrain".fasta >> undiff.aln
done < listStrainEvol.csv

## generate the "evolution tree", at the population level (instead of the strain level)
## input: /listStrainEvol.csv [list of strains to be used], hardyDiff.aln, ubiDiff.aln, undiff.aln
## output: hardyDiff.pdf, ubiDiff.pdf, undiff.pdf [plots]
Rscript s12_evolTreePlot.r


#######
### dN/dS and dS between the Hardy strains and the other strains
#######

# only for the (un)differentiated genes

cd ~/Desktop/HpGlobal_noHpGP/DIVERGENT_NO_CDS/PAML

cp -r ../alignment/highDivStrains .
cp -r ../alignment/lowDivStrains .

## as before, need to switch manually between the to folderComparisons
## the fileComparison was already generated using the script s9_createFile_dNdScomparisons.r
folderComparisons="highDivStrains"
folderComparisons="lowDivStrains"
fileComparisons="../../METADATA/pairwiseComparisons_HardyUbiquitous.txt"
nbComp=$(wc -l $fileComparisons | cut -d " " -f4)

mkdir res
echo "strain1 strain2   dN/dS (dN,dS)" > res/NG86
echo "strain1 strain2   S       N        t   kappa   omega     dN +- SE    dS +- SE" > res/YN00

for i in {1..20118..10000}; do
    bash s10_runPAML.sh $i 10000 $fileComparisons $folderComparisons > "$i".out 2> "$i"_error.out &
done

mkdir dNdS_high_HardyUbiquitous
mkdir dNdS_low_HardyUbiquitous

cat res/NG86* | uniq >> dNdS_high_HardyUbiquitous/NG86
cat res/YN00* | uniq >> dNdS_high_HardyUbiquitous/YN00

cat res/NG86* | uniq >> dNdS_low_HardyUbiquitous/NG86
cat res/YN00* | uniq >> dNdS_low_HardyUbiquitous/YN00

rm -r res; rm *.out

## the script to generate all the plots (in particular for the manuscript)
## better to run it manually than in the command line
Rscript s13_plots.r






