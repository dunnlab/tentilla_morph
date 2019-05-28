library(BAMMtools)
library(stringr)
library(readr)
setwd("~/tentilla_morph/BAMM")
tree <- read.tree("desmoneme_elongation_tree.tre")
traits <- read.csv("desmoneme_elongation.txt", sep="\t", header = F)
rownames(traits) = traits[,1]
traits <- traits[tree$tip.label,]
write_tsv(traits, "desmoneme_elongation.txt", col_names = FALSE)
is.ultrametric(tree) # check for ultrametric
is.binary.tree(tree) # check if there are no polytomies
min(tree$edge.length) # check all branch lengths are > 0
setBAMMpriors(phy = tree, traits = "desmoneme_elongation.txt")
setwd("~/tentilla_morph/")
#q()