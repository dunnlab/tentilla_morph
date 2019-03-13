library(BAMMtools)
library(stringr)
library(readr)
setwd("~/tentilla_morph/BAMM")
tree <- read.tree("haploneme_elongation_tree.tre")
traits <- read.csv("Desmoneme.length..um..txt", sep="\t", header = F)
rownames(traits) = traits[,1]
traits <- traits[tree$tip.label,]
write_tsv(traits, "Desmoneme.length..um..txt", col_names = FALSE)
is.ultrametric(tree) # check for ultrametric
is.binary.tree(tree) # check if there are no polytomies
min(tree$edge.length) # check all branch lengths are > 0
setBAMMpriors(phy = tree, traits = "Desmoneme.length..um..txt")
setwd("~/tentilla_morph/")
#q()