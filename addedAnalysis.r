# conda activate hpecospecies
# cd Desktop/hp

## in R
library(ape)

## use Tempest to root the trees
## for that, need to convert the trees from newick format to nexus format

tree <- system('ls reviews/newick_nonRootTrees', intern = TRUE)[c(2:10)]
tree <- c(tree, paste0("genes/", system('ls reviews/newick_nonRootTrees/genes', intern = TRUE)))
tree <- c(tree, paste0("whole/", system('ls reviews/newick_nonRootTrees/whole', intern = TRUE)))

tree <- unlist(strsplit(tree, split = '.nwk')) ## assume that all tree files end with .nwk

for(itree in tree) {
    tree <- read.tree(paste0("reviews/newick_nonRootTrees/", itree, ".nwk"))
    write.nexus(tree, file = paste0("reviews/nexus/", itree, ".nex"))
}

## exit R

## then need to manually run tempest
# ../../Software/TempEst_v1.5.3/bin/tempest 
## choose the file -> go into tree -> check best-fitting root
## function used: heuristic residual mean squared

## need to redo the plots
## 1) using the new rooted trees
## 2) need to change the colors to be ;ore colorblind-friendly
## see the script s13_plots.r (just need to replqce the colors and the names of the trees) from my github page
## script plots.r
