library(tidyverse)
library(seqinr)
library(msa)
library(ape)
setwd("~/Dropbox/tentilla_phylogeny2019")

#Load FASTA files
LIST_16S = list()
for(f in 1:length(list.files("sequences16S"))){
  LIST_16S[[f]] = seqinr::read.fasta(paste("./sequences16S/", list.files("sequences16S")[f], sep=""))
}
LIST_18S = list()
for(f in 1:length(list.files("sequences18S"))){
  LIST_18S[[f]] = seqinr::read.fasta(paste("./sequences18S/", list.files("sequences18S")[f], sep=""))
}

#Transform list of FASTAs into dataframes
sequences16S = lapply(LIST_16S, getSequence) %>% lapply(function(x){lapply(x,function(y){paste(unlist(y), sep="", collapse="")})}) %>% unlist()
names16S = lapply(LIST_16S, names) %>% unlist()
annots16S = lapply(LIST_16S, getAnnot) %>% unlist()
ALL16S = data.frame(names16S,annots16S, rep("16S", length(names16S)), sequences16S)
names(ALL16S) = c("Taxon","Annotation","Gene","Sequence")

sequences18S = lapply(LIST_18S, getSequence) %>% lapply(function(x){lapply(x,function(y){paste(unlist(y), sep="", collapse="")})}) %>% unlist()
sequences18S = str_replace_all(sequences18S,"[^atgc]","")
names18S = lapply(LIST_18S, names) %>% unlist()
annots18S = lapply(LIST_18S, getAnnot) %>% unlist()
ALL18S = data.frame(names18S,annots18S,rep("18S", length(names18S)),sequences18S)
names(ALL18S) = c("Taxon","Annotation","Gene","Sequence")
ALL16S = ALL16S[order(ALL16S$Taxon),]
ALL18S = ALL18S[order(ALL18S$Taxon),]
write.fasta(as.list(ALL16S$Sequence), names=paste(ALL16S$Annotation, 1:nrow(ALL16S))  , file.out = "all16S.fasta")
write.fasta(as.list(ALL18S$Sequence), names=paste(ALL18S$Annotation, 1:nrow(ALL18S))  , file.out = "all18S.fasta")

#Align by gene
ALL16S_aligned = ALL16S
ALL18S_aligned = ALL18S
ALL16S$Sequence %>% as.vector() %>% msa(method="Muscle", type="dna") %>% msaConvert(., type="seqinr::alignment") -> aln16S
ALL16S_aligned$Sequence <- aln16S$seq
write.fasta(as.list(ALL16S_aligned$Sequence), names=paste(ALL16S_aligned$Annotation, 1:nrow(ALL16S_aligned))  , file.out = "aligned_all16S.fa")
write.csv(ALL16S_aligned, "aligned_all16S.csv", row.names=F)
ALL18S$Sequence %>% as.vector() %>% msa(method="Muscle", type="dna") %>% msaConvert(., type="seqinr::alignment") -> aln18S
ALL18S_aligned$Sequence <- aln18S$seq
write.fasta(as.list(ALL18S_aligned$Sequence), names=paste(ALL18S_aligned$Annotation, 1:nrow(ALL18S_aligned))  , file.out = "aligned_all18S.fa")
write.csv(ALL18S_aligned, "aligned_all18S.csv", row.names = F)

#Read aligned files
ALL16S_aligned = read.csv("aligned_all16S.csv")
ALL16S_aligned = data.frame(ALL16S_aligned, source16S)
ALL18S_aligned = read.csv("aligned_all18S.csv")
source16S = rep("?", nrow(ALL16S_aligned))
source18S = rep("?", nrow(ALL18S_aligned))
source16S[str_detect(ALL16S_aligned$Annotation, "NULL")]<-"DunnEtAl2005++2011"
source18S[str_detect(ALL18S_aligned$Annotation, "NULL")]<-"DunnEtAl2005++2011"

taxa16S = as.character(unique(ALL16S_aligned$Taxon))
taxa18S = as.character(unique(ALL18S_aligned$Taxon))

only16Staxa = taxa16S[which(!(taxa16S %in% taxa18S))]
sharedTaxa = taxa16S[which(taxa16S %in% taxa18S)]
only18Staxa = taxa18S[which(!(taxa18S %in% taxa16S))]
ribosomal_taxa = unique(c(taxa16S,taxa18S))

#Transcriptome tree coverage
transcriptome = read.tree('~/tentilla_morph/ultimatephylo2018.tre')
transcriptome$tip.label[transcriptome$tip.label=="Cordagalma_sp"] = "Cordagalma_ordinatum"
transcriptome$tip.label[transcriptome$tip.label=="Prayidae_D27SS7"] = "Desmophyes_haematogaster"
transcriptome$tip.label[transcriptome$tip.label=="Prayidae_D27D2"] = "Craseoa_lathetica"
transcriptome$tip.label[transcriptome$tip.label=="Physonect_sp_"] = "Physonect_sp"
transcriptomeTaxa = transcriptome$tip.label
need_to_pull16S = transcriptomeTaxa[which(!(transcriptomeTaxa %in% taxa16S))]
need_to_pull18S = transcriptomeTaxa[which(!(transcriptomeTaxa %in% taxa18S))]

total_possible_taxa = length(unique(c(ribosomal_taxa,transcriptomeTaxa)))

#Pull trascriptome sequences
LIST_16S = list()
for(f in 1:length(list.files("agalma16S"))){
  LIST_16S[[f]] = seqinr::read.fasta(paste("./agalma16S/", list.files("agalma16S")[f], sep=""))
}
LIST_18S = list()
for(f in 1:length(list.files("agalma18S"))){
  LIST_18S[[f]] = seqinr::read.fasta(paste("./agalma18S/", list.files("agalma18S")[f], sep=""))
}

#Transform list of FASTAs into dataframes
sequences16S = lapply(LIST_16S, getSequence) %>% lapply(function(x){lapply(x,function(y){paste(unlist(y), sep="", collapse="")})}) %>% unlist()
names16S = lapply(LIST_16S, getAnnot) %>% unlist() %>% str_replace_all(">","")
annots16S = lapply(LIST_16S, getAnnot) %>% unlist()
agalma16S = data.frame(names16S,annots16S, rep("16S", length(names16S)), sequences16S)
names(agalma16S) = c("Taxon","Annotation","Gene","Sequence")

sequences18S = lapply(LIST_18S, getSequence) %>% lapply(function(x){lapply(x,function(y){paste(unlist(y), sep="", collapse="")})}) %>% unlist()
sequences18S = str_replace_all(sequences18S,"[^atgc]","")
names18S = lapply(LIST_18S, getAnnot) %>% unlist() %>% str_replace_all(">","")
annots18S = lapply(LIST_18S, getAnnot) %>% unlist()
agalma18S = data.frame(names18S,annots18S,rep("18S", length(names18S)),sequences18S)
names(agalma18S) = c("Taxon","Annotation","Gene","Sequence")

# #Align by gene
agalma16S_aligned = agalma16S
agalma18S_aligned = agalma18S
agalma16S$Sequence %>% as.vector() %>% msa(method="Muscle", type="dna") %>% msaConvert(., type="seqinr::alignment") -> aln16S
agalma16S_aligned$Sequence <- aln16S$seq
write.fasta(as.list(agalma16S_aligned$Sequence), names=agalma16S_aligned$Taxon  , file.out = "aligned_agalma16S.fa")
write.csv(agalma16S_aligned, "aligned_agalma16S.csv", row.names=F)
agalma18S$Sequence %>% as.vector() %>% msa(method="Muscle", type="dna") %>% msaConvert(., type="seqinr::alignment") -> aln18S
agalma18S_aligned$Sequence <- aln18S$seq
write.fasta(as.list(agalma18S_aligned$Sequence), names=agalma18S_aligned$Taxon  , file.out = "aligned_agalma18S.fa")
write.csv(agalma18S_aligned, "aligned_agalma18S.csv", row.names = F)
