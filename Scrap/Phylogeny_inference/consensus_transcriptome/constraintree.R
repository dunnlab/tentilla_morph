library(tidyverse)
library(ape)
library(phytools)
library(seqinr)

setwd("~/Dropbox/tentilla_phylogeny2019/consensus_transcriptome")

rrna <- seqinr::read.fasta("../all18S_curated_sublist.fasta")
rrna_names <- names(rrna)

treelist = list()
for(tr in 1:length(list.files("trees"))){
  treelist[[tr]] = ape::read.tree(paste("./trees/", list.files("trees")[tr], sep=""))
}
root(treelist[[1]],25) -> treelist[[1]]
root(treelist[[2]],25) -> treelist[[2]]
root(treelist[[3]],25) -> treelist[[3]]
root(treelist[[4]],17) -> treelist[[4]]
root(treelist[[5]],17) -> treelist[[5]]
root(treelist[[6]],17) -> treelist[[6]]
str_consensus <- ape::consensus(treelist)
str_consensus$tip.label[which(str_consensus$tip.label == "Prayidae_D27D2")] <- "Craseoa_lathetica"
str_consensus$tip.label[which(str_consensus$tip.label == "Prayidae_D27SS7")] <- "Desmophyes_haematogaster"
str_consensus$tip.label[which(str_consensus$tip.label == "Physonect_sp_")] <- "Physonect_sp"
consensus_tl <- str_consensus$tip.label

missingtaxa <- rrna_names[which(!(rrna_names %in% consensus_tl))]
extrataxa <- consensus_tl[which(!(consensus_tl %in% rrna_names))]
pending_transcriptome_taxa[]

