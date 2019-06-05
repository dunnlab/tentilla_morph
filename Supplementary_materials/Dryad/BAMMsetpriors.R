library(BAMMtools)
library(stringr)
library(readr)
#setwd("~/tentilla_morph/Supplementary_Materials/Dryad")
tree <- read.tree("haploneme_elongation_tree.tre")
traits <- read.csv("haploneme_elongation.txt", sep="\t", header = F)
rownames(traits) = traits[,1]
traits <- traits[tree$tip.label,]
write_tsv(traits, "haploneme_elongation.txt", col_names = FALSE)
is.ultrametric(tree) # check for ultrametric
is.binary.tree(tree) # check if there are no polytomies
min(tree$edge.length) # check all branch lengths are > 0
setBAMMpriors(phy = tree, traits = "haploneme_elongation.txt")
#q()