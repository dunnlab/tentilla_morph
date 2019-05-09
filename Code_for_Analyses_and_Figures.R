
## LOAD PACKAGES ##
#General
library(tidyverse)
library(stringr)
library(readr)
library(reshape2)
## Biological
library(Rphylip)
library(arbutus)
library(vegan)
library(ape)
library(phangorn)
library(phytools)
library(OUwie)
library(picante)
library(geiger)
library(phylobase)
library(fields)
library(phylosignal)
library(geomorph)
library(phylopath)
library(phylolm)
library(adegenet)
## Graphics
library(FactoMineR)
library(factoextra)
library(corrplot)
library(BAMMtools)
library(gridExtra)
library(colorRamps)

# Set paths to input data
setwd("~/tentilla_morph")

#Load raw data
read.csv("byslide.csv") -> numbers
numbers$Species = as.character(numbers$Species)
categorical <- read.csv("Homolog_Categorical.csv")[,-2]
rownames(categorical) = categorical$Species

#Correct species spellings
numbers$Species[which(numbers$Species == "Agalma okeni")] <- "Agalma okenii"
numbers$Species[which(numbers$Species == "Forskalia edwardsi")] <- "Forskalia edwardsii"
numbers$Species[which(numbers$Species == "Rhizophysa eysenhardti")] <- "Rhizophysa eysenhardtii"

#Merge Nanomias
numbers$Species[which(numbers$Species == "Nanomia cara")] <- "Nanomia sp"
numbers$Species[which(numbers$Species == "Nanomia bijuga")] <- "Nanomia sp"

#Polish fields
castnumbers = numbers
castnumbers[castnumbers=="needDIC"] <- NA
castnumbers[castnumbers=="needConfocal"] <- NA
castnumbers[castnumbers=="needMature"] <- NA
castnumbers[castnumbers=="needReslide"] <- NA
castnumbers[castnumbers==-1] <- NA
castnumbers[,c(-1,-2,-3)] <- sapply(castnumbers[,c(-1,-2,-3)],as.character)
castnumbers[,c(-1,-2,-3)] <- sapply(castnumbers[,c(-1,-2,-3)],as.numeric)

#Define functions to deal with NAs and decimals
mean.na <- function(a){mean(a, na.rm = TRUE)}
var.na <- function(a){var(a, na.rm = TRUE)}
round3 <- function(a){round(a, 3)}

#Compute morphometric ratios
coiledness = castnumbers$Cnidoband.free.length..um./castnumbers$Cnidoband.length..um.
heteroneme_elongation = castnumbers$Heteroneme.free.length..um./castnumbers$Heteroneme.width..um.
haploneme_elongation = castnumbers$Haploneme.free.length..um./castnumbers$Haploneme.width..um.
desmoneme_elongation = castnumbers$Desmoneme.length..um./castnumbers$Desmoneme.width..um.
rhopaloneme_elongation = castnumbers$Rhopaloneme.length..um./castnumbers$Rhopaloneme.width..um.
heteroneme_shaft_extension = castnumbers$Heteroneme.free.length..um./castnumbers$Heteroneme.shaft.free.length..um.
Heteroneme_to_CB = castnumbers$Heteroneme.free.length..um./(castnumbers$Cnidoband.free.length..um.+0.001)
total_heteroneme_volume = castnumbers$Heteroneme.volume..um3.*castnumbers$Heteroneme.number
total_haploneme_volume = (castnumbers$Cnidoband.free.length..um./castnumbers$Haploneme.width..um.)*((4*pi/3)*0.5*castnumbers$Haploneme.free.length..um.*((0.5*castnumbers$Haploneme.width..um.)^2))
cnidomic_index = log(total_heteroneme_volume + total_haploneme_volume)
cnidomic_index[which(is.na(cnidomic_index))] <- log(total_haploneme_volume[which(is.na(cnidomic_index))])
cnidomic_index[which(is.na(cnidomic_index))] <- log(castnumbers[which(is.na(cnidomic_index)), "Heteroneme.volume..um3."])
cnidomic_index[cnidomic_index==-Inf]<-NA

#Define surface and volume functions
surface_ellipsoid <- function(L, W){
  a <- L/2
  b <- W/2
  c = b
  S <- 4*pi*((((a*b)^1.6)+((a*c)^1.6)+((c*b)^1.6))/3)^(1/1.6)
  return(S)
}
volume_ellipsoid <- function(L, W){
  a <- L/2
  b <- W/2
  c = b
  V <- (4/3)*pi*a*b*c
  return(V)
}

SAV_haploneme = surface_ellipsoid(castnumbers$Haploneme.free.length..um., castnumbers$Haploneme.width..um.)/volume_ellipsoid(castnumbers$Haploneme.free.length..um., castnumbers$Haploneme.width..um.)
names(SAV_haploneme) = castnumbers$Species
SAV_heteroneme = surface_ellipsoid(castnumbers$Heteroneme.free.length..um., castnumbers$Heteroneme.width..um.)/volume_ellipsoid(castnumbers$Heteroneme.free.length..um., castnumbers$Heteroneme.width..um.)
names(SAV_heteroneme) = sharedmeans$Species

#Combine all morphometric variables
morphometrics = data.frame(castnumbers$slide_id, castnumbers$Species,coiledness, heteroneme_elongation, haploneme_elongation, desmoneme_elongation, rhopaloneme_elongation,heteroneme_shaft_extension, Heteroneme_to_CB, total_heteroneme_volume, total_haploneme_volume, cnidomic_index, SAV_haploneme)
names(morphometrics)[1:2] = c("slide_id"," Species")

#Combine morphometrics with regular characters
castnumbers = data.frame(castnumbers,morphometrics[,c(-1,-2)])

#Logtransform variables which are not normal
castlogs = castnumbers
non_normal = apply(castnumbers[,c(-1,-2,-3)], 2, shapiro.test) %>% lapply(function(x){x$p.value}) %>% unlist() 
non_normal = which(as.vector(non_normal) < 0.05)+3
castlogs[,non_normal] <- sapply(castnumbers[,non_normal], log)
castlogs[castlogs==-Inf]<-NA

#Means, varinces, standard errors
castmeans <- aggregate(. ~  Species, data = castnumbers[,c(-1,-3)], mean.na, na.action = na.pass)
castmean_logs <- aggregate(. ~  Species, data = castlogs[,c(-1,-3)], mean.na, na.action = na.pass)
castvariances <- aggregate(. ~  Species, data = castnumbers[,c(-1,-3)], FUN = var.na, na.action = na.pass) #[,-2]
castvariance_logs <- aggregate(. ~  Species, data = castlogs[,c(-1,-3)], FUN = var.na, na.action = na.pass)
castvariances[is.na(castvariances)] <- 0
castvariance_logs[is.na(castvariance_logs)] <- 0
castN <- aggregate(. ~  Species, data = castnumbers[,c(-1,-3)], length, na.action = na.pass) 
castN = data.frame(castN$Species, rowMeans(castN[,-1])) #[,-2]
names(castN) = c("Species", "N_specimens")
castSE <- cbind(castvariances$Species, sqrt(castvariances[,-1])/sqrt(castN$N_specimens))
castSE_logs <- cbind(castvariance_logs$Species, sqrt(castvariance_logs[,-1])/sqrt(castN[,-1]))
names(castSE)[1] <- "Species"
names(castSE_logs)[1] <- "Species"

#Make overview of data and heatmap
Ses = castSE[,-1]
Ses[Ses == 0] = NA
morphdata = cbind(castN, castmeans[,-1], Ses)
morphdata = morphdata[,c(1:3,33,4,34,5,35,6,36,7,37,8,38,9,39,10,40,11,41,12,42,13,43,14,44,15,45,16,46,17,47,18,48,19,49,20,50,21,51,22,52,23,53,24,54,25,55,26,56,27,57,28,58,29,59,30,60,31,61,32,62)]
names(morphdata) = str_replace_all(names(morphdata),".1","_SE")
write.csv(morphdata, "characterdata.csv")

heatdata = as.matrix(castmean_logs[,-1])
rownames(heatdata) = castmean_logs$Species
heatdata[is.nan(heatdata)]<- -1
hcolors = grDevices::terrain.colors(20)
hcolors<-c(rep("#000000FF",10), hcolors)
heatmap(heatdata, scale = "column", cexCol = 0.2, col=hcolors, keep.dendro = T)

#Load phylogenetic tree
consensus = read.nexus("TimeTree_CT_truerho.tre") %>% drop.tip(56:61)
consensus$tip.label = str_replace_all(consensus$tip.label,"_"," ")
#Switch Nanomia bijuga for Nanaomia sp
consensus$tip.label[which(consensus$tip.label == "Nanomia bijuga")] <- "Nanomia sp"

#Prune quant matrix to tree species
matrix = castnumbers[which(!is.na(castnumbers$Species)),]
matrix_logs = castlogs[which(!is.na(castlogs$Species)),]
sharedspp = matrix$Species[which(matrix$Species %in% consensus$tip.label)]
sharedmatrix = matrix[which(matrix$Species %in% sharedspp),]
sharedlogs = matrix_logs[which(matrix$Species %in% sharedspp),]
sharedmeans = castmeans[which(castmeans$Species %in% sharedspp),]
sharedmean_logs = castmean_logs[which(castmean_logs$Species %in% sharedspp),]
sharedvars = castvariances[which(castvariances$Species %in% sharedspp),]
sharedvar_logs = castvariance_logs[which(castvariance_logs$Species %in% sharedspp),]

#Prune categorical matrix to tree species
cat_rowNAs = apply(categorical[,-1], 1, function(x) sum(is.na(x)))
catmatrix = categorical[which(cat_rowNAs<1),] %>% .[,-1]
sharedspp_cat = rownames(catmatrix)[which(rownames(catmatrix) %in% consensus$tip.label)]
sharedcategorical = catmatrix[which(rownames(catmatrix) %in% sharedspp_cat),]

#Prune tree to quant matrix species
nodatatipnames = consensus$tip.label[which(!(consensus$tip.label %in% sharedmatrix$Species))]
nodatatips = c(1:length(consensus$tip.label))[which(consensus$tip.label %in% nodatatipnames)]
prunedtree = drop.tip(consensus, nodatatips)
prunedtree$tip.label
plot(prunedtree)
#ultram = chronos(prunedtree)
ultram = prunedtree
plot(ultram)

#Prune tree to categorical matrix species
nodatatipnames_cat = consensus$tip.label[which(!(consensus$tip.label %in% sharedspp_cat))]
nodatatips_cat = c(1:length(consensus$tip.label))[which(consensus$tip.label %in% nodatatipnames_cat)]
cat_tree = drop.tip(consensus, nodatatips_cat)
cat_tree$tip.label
plot(cat_tree)
#ultram_cat = chronos(cat_tree)
ultram_cat = cat_tree
plot(ultram_cat)
sharedcategorical = sharedcategorical[match(ultram_cat$tip.label, rownames(sharedcategorical)),]
cprunedmatrix = sharedcategorical[which(rownames(sharedcategorical)%in%rownames(sharedmatrix)),] %>% .[which(sapply(.,function(x) length(unique(x)))>1)]
sharedbinary = sharedcategorical
sharedbinary$Haploneme.type = as.character(sharedbinary$Haploneme.type)
sharedbinary[sharedbinary=="Isorhizas"] = 0
sharedbinary[sharedbinary=="Anisorhizas"] = 1
sharedbinary$Haploneme.type = as.numeric(sharedbinary$Haploneme.type)
sharedbinary$Heteroneme.type = as.character(sharedbinary$Heteroneme.type)
sharedbinary[sharedbinary=="Stenotele"] = 0
sharedbinary[sharedbinary=="Microbasic mastigophore"] = 1
sharedbinary[sharedbinary=="Eurytele"] = NA
sharedbinary$Heteroneme.type = as.numeric(sharedbinary$Heteroneme.type)

#Purely numerical character sets without slide or species
Q_sharedmatrix = sharedmatrix[,c(-1,-2,-3)]
rownames(sharedmeans) = sharedmeans$Species
Q_sharedmeans = sharedmeans[,-1]
rownames(sharedvars) = sharedvars$Species
Q_sharedvars = sharedvars[,-1]
Q_sharedmean_logs = sharedmean_logs[,-1]
rownames(Q_sharedmean_logs) = sharedmean_logs$Species

### COMPARATIVE ANALYSES ###

phylosignals = as.data.frame(matrix(ncol=3, nrow=ncol(sharedmean_logs[,-1])))
#phylosignals = as.data.frame(matrix(ncol=3, nrow=ncol(sharedmeans[,-1])))
names(phylosignals) = c("K", "P", "Ntaxa")
for(i in 2:ncol(sharedmean_logs)){
  CH_I=as.numeric(sharedmean_logs[,i])
  #CH_I=as.numeric(sharedmeans[,i])
  names(CH_I) = sharedmean_logs$Species
  #names(CH_I) = sharedmeans$Species
  SE_I = as.numeric(castSE[,i]) %>% log()
  SE_I[which(SE_I==-Inf)] = 0
  names(SE_I) = castSE$Species
  CH_I= CH_I[!is.na(CH_I)]
  SE_I= SE_I[!is.na(SE_I)]
  SE_I = SE_I[which(names(SE_I)%in%names(CH_I))]
  CH_I = CH_I[which(names(CH_I)%in%names(SE_I))]
  treeI = drop.tip(ultram,which(!(ultram$tip.label %in% names(CH_I))))
  class(treeI) = "phylo"
  rownames(phylosignals)[i-1] <- names(sharedmean_logs)[i]
  PSIG <- phylosig(treeI, CH_I, se = SE_I, test=T)
  phylosignals[i-1,1] <- PSIG$K
  phylosignals[i-1,2] <- PSIG$P
  phylosignals[i-1,3] <- length(CH_I)
  phylosignals$K = round(phylosignals$K,3)
}
#write.csv(phylosignals, "PhylosignalsWSE_log.csv")

#Model support
AICdf = as.data.frame(matrix(ncol=6,nrow=ncol(Q_sharedmean_logs)))
colnames(AICdf) = c("Variable", "white_noise", "starBM", "BM", "EB", "OU")
for(c in 1:ncol(Q_sharedmean_logs)){
  startree <- rescale(ultram, "lambda", 0)
  C = Q_sharedmean_logs[,c]
  names(C) = rownames(Q_sharedmean_logs)
  C = C[!is.na(C)]
  Ctree = drop.tip(ultram, which(!(ultram$tip.label %in% names(C))))
  startree = drop.tip(startree, which(!(startree$tip.label %in% names(C))))
  Cse = castSE[,c+1] %>% log() %>% abs()
  Cse[which(Cse == Inf)] <- 0
  names(Cse) = castSE$Species
  Cse = Cse[which(names(Cse) %in% names(C))]
  model_matrix = matrix("NA", nrow = 5, ncol = 3)
  colnames(model_matrix) = c("aicc","aicc_best","dAICc")
  row.names(model_matrix) = c("white", "starBM", "BM", "EB", "OU")
  for(j in 1:dim(model_matrix)[1]){
    if(j==2){
      temp_model = fitContinuous(startree, C, model="BM", SE = Cse)$opt
    }
    else{
      temp_model = fitContinuous(Ctree, C, model=row.names(model_matrix)[j], SE = Cse)$opt
    }
    model_matrix = apply(model_matrix,2, as.numeric)
    row.names(model_matrix) = c("white", "starBM", "BM", "EB", "OU")
    model_matrix[j, "aicc"] <- temp_model$aicc
  }
  model_matrix[,"aicc_best"] <- min(model_matrix[,"aicc"])
  model_matrix[,"dAICc"] <- model_matrix[, "aicc"] - model_matrix[j, "aicc_best"]
  print(names(Q_sharedmean_logs)[c])
  string_c <- c(names(Q_sharedmean_logs)[c], model_matrix[,3])
  names(string_c) = colnames(AICdf)
  AICdf[c,] <- string_c
}
AICdf[,2:6] = apply(AICdf[,2:6], 2, as.numeric) %>% apply(2, round3)
#write.csv(AICdf, "log_model_support.csv")

#Model adequacy
worthy_models = AICdf[which(AICdf$white_noise != 0 & AICdf$starBM != 0),]
MAD = as.data.frame(matrix(ncol = 6, nrow = nrow(worthy_models)))
names(MAD) = c("msig", "cvar", "svar", "sasr", "shgt", "dcfd")
rownames(MAD) = worthy_models$Variable
for(m in 1:nrow(worthy_models)){
  C = Q_sharedmean_logs[,which(names(Q_sharedmean_logs) == worthy_models$Variable[m])]
  names(C) = rownames(Q_sharedmean_logs)
  C = C[!is.na(C)]
  Ctree = drop.tip(ultram, which(!(ultram$tip.label %in% names(C))))
  C = C[match(Ctree$tip.label, names(C))]
  class(Ctree)="phylo"
  FC <- fitContinuous(Ctree,C,model=worthy_models$Best_model[m])
  UTC <- make_unit_tree(FC)
  picstat_data <- calculate_pic_stat(UTC)
  sim <- simulate_char_unit(UTC)
  picstat_sim <- calculate_pic_stat(sim)
  compare_pic_stat(picstat_data, picstat_sim) %>% .$p.values -> MAD[m,]
}
phy_models = data.frame(worthy_models,MAD)
phy_models[,7:12] = apply(phy_models[,7:12], 2,round3)
#write.csv(phy_models, "model_adequacy.csv")

#Nematocyst shape evolution - phylomorphospace figure
het_el = sharedmean_logs$heteroneme_elongation
names(het_el) = sharedmean_logs$Species
hap_el = sharedmean_logs$haploneme_elongation
names(hap_el) = sharedmean_logs$Species
het_el = het_el[which(!is.na(het_el))]
hap_el = hap_el[which(!is.na(hap_el))]
hap_el = hap_el[which(names(hap_el) %in% names(het_el))]
het_el = het_el[which(names(het_el) %in% names(hap_el))]
het_el = het_el[match(names(hap_el), names(het_el))]
elon_data = cbind(het_el, hap_el) %>% as.data.frame()
names(elon_data) = c("Heteroneme elongation (log um)", "Haploneme elongation (log um)")
elon_tree = drop.tip(ultram, which(!(ultram$tip.label %in% names(hap_el))))
phylomorphospace(elon_tree, elon_data, label="horizontal")

## SIMMAPS used to make categorical evolution figure ##
par(ask=F)
tentilla = sharedcategorical$Tentilla
names(tentilla) = rownames(sharedcategorical)
prox_het = sharedcategorical$Proximal.heteronemes
names(prox_het) = rownames(sharedcategorical)
desmo = sharedcategorical$Desmonemes
names(desmo) = rownames(sharedcategorical)
rhopalo = sharedcategorical$Rhopalonemes
names(rhopalo) = rownames(sharedcategorical)
dyn_cnido = sharedcategorical$Dynamic.cnidoband
names(dyn_cnido) = rownames(sharedcategorical)
elastic = sharedcategorical$Elastic.strand
names(elastic) = rownames(sharedcategorical)
distal_desmo = sharedcategorical$Distal.cnidoband.desmonemes
names(distal_desmo) = rownames(sharedcategorical)
coiled = sharedcategorical$Coiled.tentilla
names(coiled) = rownames(sharedcategorical)
heterotype = sharedcategorical$Heteroneme.type
heterotype=as.character(heterotype)
names(heterotype) = rownames(sharedcategorical)
haplotype = sharedcategorical$Haploneme.type
haplotype=as.character(haplotype)
names(haplotype) = rownames(sharedcategorical)

Simmap_list = list()
###SIMMAP Tentilla:
make.simmap(ultram_cat, tentilla, nsim = 100) -> tentilla_sim
Simmap_list[[1]] <- tentilla_sim
plotTree(ultram_cat, lwd = 4)
tentilla_sim %>% plotSimmap(lwd = 4, add = T)
colors = c("black", "red")
names(colors) = c("Present", "Absent")
nodelabels(pie=(describe.simmap(tentilla_sim, plot=F)$ace) ,piecol=colors,cex=0.35)
add.simmap.legend(colors = colors, x=0.6*par()$usr[1],y=0.3*par()$usr[4],prompt=FALSE)
densityMap(tentilla_sim)

###SIMMAP Proximal Heteronemes:
make.simmap(ultram_cat, prox_het, nsim = 100) -> prox_het_sim
Simmap_list[[2]] <- prox_het_sim
plotTree(ultram_cat, lwd = 4)
prox_het_sim %>% plotSimmap(lwd = 4, add = T)
colors = c("black", "red")
names(colors) = c("Present", "Absent")
nodelabels(pie=(describe.simmap(prox_het_sim, plot=F)$ace) ,piecol=colors,cex=0.35)
add.simmap.legend(colors = colors, x=0.6*par()$usr[1],y=0.3*par()$usr[4],prompt=FALSE)
densityMap(prox_het_sim)

###SIMMAP Desmonemes:
make.simmap(ultram_cat, desmo, nsim = 100) -> desmo_sim
Simmap_list[[3]] <- desmo_sim
plotTree(ultram_cat, lwd = 4)
desmo_sim %>% plotSimmap(lwd = 4, add = T)
colors = c("black", "red")
names(colors) = c("Absent", "Present")
nodelabels(pie=(describe.simmap(desmo_sim, plot=F)$ace) ,piecol=colors,cex=0.35)
add.simmap.legend(colors = colors, x=0.6*par()$usr[1],y=0.3*par()$usr[4],prompt=FALSE)
densityMap(desmo_sim)

###SIMMAP Rhopalonemes:
make.simmap(ultram_cat, rhopalo, nsim = 100) -> rhopalo_sim
Simmap_list[[4]] <- rhopalo_sim
plotTree(ultram_cat, lwd = 4)
rhopalo_sim %>% plotSimmap(lwd = 4, add = T)
colors = c("black", "red")
names(colors) = c("Absent", "Present")
nodelabels(pie=(describe.simmap(rhopalo_sim, plot=F)$ace) ,piecol=colors,cex=0.35)
add.simmap.legend(colors = colors, x=0.6*par()$usr[1],y=0.3*par()$usr[4],prompt=FALSE)
densityMap(rhopalo_sim)

###SIMMAP Dynamic Cnidoband:
make.simmap(ultram_cat, dyn_cnido, nsim = 100) -> dyn_cnido_sim
Simmap_list[[5]] <- dyn_cnido_sim
plotTree(ultram_cat, lwd = 4)
dyn_cnido_sim %>% plotSimmap(lwd = 4, add = T)
colors = c("black", "red")
names(colors) = c("Absent", "Present")
nodelabels(pie=(describe.simmap(dyn_cnido_sim, plot=F)$ace) ,piecol=colors,cex=0.35)
add.simmap.legend(colors = colors, x=0.6*par()$usr[1],y=0.3*par()$usr[4],prompt=FALSE)
densityMap(dyn_cnido_sim)

###SIMMAP Elastic Strand:
make.simmap(ultram_cat, elastic, nsim = 100) -> elastic_sim
Simmap_list[[6]] <- elastic_sim
plotTree(ultram_cat, lwd = 4)
elastic_sim %>% plotSimmap(lwd = 4, add = T)
colors = c("black", "red")
names(colors) = c("Absent", "Present")
nodelabels(pie=(describe.simmap(elastic_sim, plot=F)$ace) ,piecol=colors,cex=0.35)
add.simmap.legend(colors = colors, x=0.6*par()$usr[1],y=0.3*par()$usr[4],prompt=FALSE)
densityMap(elastic_sim)

###SIMMAP Distal CB Desmonemes:
make.simmap(ultram_cat, distal_desmo, nsim = 100) -> distal_desmo_sim
Simmap_list[[7]] <- distal_desmo_sim
plotTree(ultram_cat, lwd = 4)
distal_desmo_sim %>% plotSimmap(lwd = 4, add = T)
colors = c("black", "red")
names(colors) = c("Absent", "Present")
nodelabels(pie=(describe.simmap(distal_desmo_sim, plot=F)$ace) ,piecol=colors,cex=0.35)
add.simmap.legend(colors = colors, x=0.6*par()$usr[1],y=0.3*par()$usr[4],prompt=FALSE)
densityMap(distal_desmo_sim)

###SIMMAP Tentilla Coiledness:
make.simmap(ultram_cat, coiled, nsim = 100) -> coiled_sim
Simmap_list[[8]] <- coiled_sim
plotTree(ultram_cat, lwd = 4)
coiled_sim %>% plotSimmap(lwd = 4, add = T)
colors = c("black", "red")
names(colors) = c("Absent", "Present")
nodelabels(pie=(describe.simmap(coiled_sim, plot=F)$ace) ,piecol=colors,cex=0.35)
add.simmap.legend(colors = colors, x=0.6*par()$usr[1],y=0.3*par()$usr[4],prompt=FALSE)
densityMap(coiled_sim)

###SIMMAP Heteroneme type:
heterotype = heterotype[heterotype!=""]
HTtree = drop.tip(ultram_cat, which(!(ultram_cat$tip.label %in% names(heterotype))))
make.simmap(HTtree, heterotype, nsim = 100) -> heterotype_sim
Simmap_list[[9]] <- heterotype_sim
plotTree(HTtree, lwd = 4)
heterotype_sim %>% plotSimmap(lwd = 4, add = T)
colors = c("black", "red", "green")
names(colors) = c("Eurytele", "Microbasic mastigophore", "Stenotele")
nodelabels(pie=(describe.simmap(heterotype_sim, plot=F)$ace) ,piecol=colors,cex=0.35)
add.simmap.legend(colors = colors, x=0.6*par()$usr[1],y=0.3*par()$usr[4],prompt=FALSE)

###SIMMAP Haploneme type:
haplotype = haplotype[haplotype!=""]
HTtree = drop.tip(ultram_cat, which(!(ultram_cat$tip.label %in% names(haplotype))))
make.simmap(HTtree, haplotype, nsim = 100) -> haplotype_sim
Simmap_list[[10]] <- haplotype_sim
plotTree(HTtree, lwd = 4)
haplotype_sim %>% plotSimmap(lwd = 4, add = T)
colors = c("black", "red")
names(colors) = c("Isorhizas", "Anisorhizas")
nodelabels(pie=(describe.simmap(haplotype_sim, plot=F)$ace) ,piecol=colors,cex=0.35)
add.simmap.legend(colors = colors, x=0.6*par()$usr[1],y=0.3*par()$usr[4],prompt=FALSE)
densityMap(haplotype_sim)

## Character correlations ##
C = sharedmatrix[,c(-1,-3)] %>% .[which(!is.na(rowSums(.[,-1]))),]
Cspecies = C$Species %>% as.character()
C = as.matrix(C[,-1])
rownames(C) = Cspecies
Ctree = drop.tip(ultram,which(!(ultram$tip.label %in% rownames(C))))
class(Ctree) = "phylo"
PICi = Rcontrast(tree=Ctree, X=C, path="phylip-3.695/exe", cleanup=TRUE)

#Correlation visualizations R2
phy_corr = PICi$VarA.Correlations
intra_corr = PICi$VarE.Correlations
phy_corr[upper.tri(phy_corr)] <- NA
combicorr = phy_corr
combicorr[upper.tri(combicorr)] = intra_corr[upper.tri(intra_corr)]
rownames(combicorr) = colnames(C)
colnames(combicorr) = colnames(C)
corrplot(combicorr, diag=F, tl.cex = 0.4, tl.col="black") #phylo correlations vs intraspecific correlations
PICOLS = combicorr
PICOLS[upper.tri(PICOLS)]=cor(C)[upper.tri(cor(C))]
corrplot(PICOLS, diag=F, tl.cex = 0.4, tl.col="black") #phylo correlations vs regular correlations for figure

#Scatterplot phylo vs regular correlations
cbind(as.vector(phy_corr), as.vector(cor(C))) %>% as.data.frame()->phyreg
names(phyreg)<-c("Phylo", "Reg")
ggplot(phyreg, aes(x=Reg, y=Phylo, color=(Phylo+Reg)/2)) + geom_point() + geom_hline(yintercept = 0)  + geom_vline(xintercept = 0) + theme_bw()
abline(h=0)
abline(v=0)

## PCA ##

#Using simple characters
PCA(raw_matrix_notf) -> Pca_raw
Pca_raw %>% fviz_contrib(choice="var", axes=1, sort.val="desc")
Pca_raw %>% fviz_contrib(choice="var", axes=2, sort.val="desc")
Pca_raw %>% fviz_pca_biplot( col.var="contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE) #For figure A
phylomorphospace(tree=raw_tree, Pca_raw$ind$coord[,1:2], label = "horizontal", xlab = "PC1", ylab = "PC2") #For figure B
multiPhylosignal(Pca_raw$ind$coord, raw_tree)
physignal(Pca_raw$ind$coord, raw_tree)

#Using Morphometric characters
PCA(compound_matrix) -> Pca_compound
Pca_compound %>% fviz_contrib(choice="var", axes=1, sort.val="desc")
Pca_compound %>% fviz_contrib(choice="var", axes=2, sort.val="desc")
Pca_compound %>% fviz_pca_biplot( col.var="contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE)
phylomorphospace(tree=compound_tree, Pca_compound$ind$coord[,1:2], label = "horizontal", xlab = "PC1", ylab = "PC2")
physignal(Pca_compound$ind$coord, compound_tree)
multiPhylosignal(Pca_compound$ind$coord, compound_tree)

## BAMM ##

traitlist = list()
for(i in 2:ncol(sharedmean_logs)){
  traitlist[[i]] <- sharedmean_logs[,c(1,i)]
  names(traitlist)[i] = names(sharedmean_logs)[i]
}

traitlist <- lapply(traitlist, function(x){x<-x[which(!is.nan(x[,2])),]})

for(i in 2:length(traitlist)){
  treei = drop.tip(ultram, which(!(ultram$tip.label %in% traitlist[[i]]$Species)))
  treei$tip.label = str_replace_all(treei$tip.label, " ", "_")
  write.tree(treei, file = paste("BAMM/",names(traitlist[[i]])[2], "_tree.tre", sep=""))
}

for(i in 2:length(traitlist)){
  traitlist[[i]]$Species = str_replace_all(traitlist[[i]]$Species, " ", "_")
  write.table(traitlist[[i]],file = paste("BAMM/",names(traitlist[[i]])[2], ".txt", sep=""), sep="\t", col.names = F, row.names = F, quote = F)
}

#See BAMM/BAMMsetpriors.R
#See BAMM/BAMMplots.R
#See BAMM/BAMMconfiguration.txt

## Retrieve diet data ##
#Retrieve diet info from literature BINARY
GC = read.csv("Cmerged.csv", header = T, sep=',')[,c(2,4,5,6)] %>% .[which(grepl("diet",.$character) & .$state==1),]
GC$character = factor(GC$character, levels=unique(GC$character))
GC$species = as.character(GC$species)
#Fix typos#
GC$species[which(GC$species == "Nanomia bijuga")] <- "Nanomia sp"
GC$species[which(GC$species == "Rhizophysa eyesenhardti")] <- "Rhizophysa eysenhardtii"
GC$species[which(GC$species == "Agalma okeni")] <- "Agalma okenii"

#Prune to tree species#
GC = GC[which(GC$species %in% ultram$tip.label),]
GC = split(GC,GC$character)
nrowGC = purrr::map(GC,nrow) %>% as.numeric()
GC = GC[which(nrowGC>2)]
GC = purrr::map(GC,unique)
diet= matrix(ncol=length(GC),nrow=length(unique(sharedmeans$Species))) %>% as.data.frame()
names(diet) = names(GC)
rownames(diet) = unique(sharedmeans$Species)
for(E in GC){
  print(E$species)
  for(S in E$species){
    diet[which(rownames(diet) == S),as.character(unique(E$character))] = 1
  }
}
diet[is.na(diet)] <- 0
diet = diet[which(rowSums(diet)>0),which(colSums(diet)<nrow(diet))]

#Add personal observations of the authors
cladeB <- matrix(rep(c(0,1,0,0,0,0,0,0,0,0),3),nrow=10, ncol = 3) %>% t() %>% as.data.frame()
rownames(cladeB) = c("Erenna richardi", "Erenna sirena", "Stephanomia amphytridis")
krilleaters = matrix(rep(c(0,0,0,1,0,0,1,0,0,0),4),nrow=10, ncol = 4) %>% t() %>% as.data.frame()
rownames(krilleaters) = c("Praya dubia", "Resomia ornicephala", "Lychnagalma utricularia", "Bargmannia amoena")
Gelatinous_diet = rep(0,nrow(diet))
diet = cbind(diet, Gelatinous_diet)
names(diet)[ncol(diet)] = "Gelatinous diet"
gelateaters = matrix(rep(c(0,0,0,0,0,0,0,0,0,1),1),nrow=10, ncol = 1) %>% t() %>% as.data.frame()
names(gelateaters) = names(diet)
rownames(gelateaters) = c("Apolemia rubriversa")
names(gelateaters)=names(diet)
names(krilleaters)=names(diet)
names(cladeB)=names(diet)
diet = rbind(diet, cladeB, krilleaters, gelateaters)
#Decapod diet column actually encompasses decapods, krill, and mysids (large crustaceans/shrimp like animals)

# Retrive ROV annotation data
VARS <- read.csv("Choy2017.tsv", sep="\t")
VARS_curated = VARS[which(VARS$Pred_lowest_tax %in% ultram$tip.label | VARS$Pred_lowest_tax=="Nanomia bijuga"),]
VARS_cast = acast(VARS_curated, Pred_lowest_tax~Prey_order, fun.aggregate = length)

#Prune morphological matrix to species in diet
dprunedmatrix = sharedmeans[which(sharedmeans$Species%in%rownames(diet)),]
dprunedmatrix_logs = sharedmean_logs[which(sharedmean_logs$Species%in%rownames(diet)),]
#Prune tree to diet species
dprunedTree = drop.tip(ultram, which(!(ultram$tip.label %in% rownames(diet))))

## Retrieve prey selectivity ##
selectivity = read.csv("PurcellSiphonophoresDiet.csv", header=T, sep=",")[,c(1,31:39)]
selectivity$Species <- as.character(selectivity$Species)
selectivity$Species[which(selectivity$Species == "Nanomia bijuga")] <- "Nanomia sp"
selectivity$Species[which(selectivity$Species == "Rhizophysa eyesenhardti")] <- "Rhizophysa eysenhardtii"
selectivity = aggregate(. ~ Species, data=selectivity, FUN=mean)
rownames(selectivity) = selectivity$Species
selectivity=selectivity[,-1]
# Prune morphological data to selectivity tree
Sprunedmatrix = sharedmeans[which(sharedmeans$Species%in%rownames(selectivity)),]
Sprunedmatrix_logs = sharedmean_logs[which(sharedmean_logs$Species%in%rownames(selectivity)),]
selectivity = selectivity[which(rownames(selectivity) %in% sharedmeans$Species),]
# Prune tree to selectivity species
Sprunedtree = drop.tip(ultram, which(!(ultram$tip.label %in% rownames(selectivity))))

#Retrieve copepod prey length from literature
quantDiet = read.csv("Qmerged.csv", header = T, sep=',')[,c(2,4,6)]
preylength = quantDiet[which(grepl("prey",quantDiet$character)),] %>% .[,c(1,3)]
names(preylength) = c("Species","Copepod prey length (mm)")
PL_pruned_matrix = sharedmean_logs[which(sharedmean_logs$Species %in% preylength$Species),] %>% data.frame(preylength[which(preylength$Species %in% sharedmean_logs$Species),2])
PL_pruned_SEs = sharedvar_logs[which(sharedvar_logs$Species %in% preylength$Species),] %>% data.frame(preylength[which(preylength$Species %in% sharedmean_logs$Species),2])
names(PL_pruned_matrix)[ncol(PL_pruned_matrix)] <- "Copepod prey length (mm)"
names(PL_pruned_SEs)[ncol(PL_pruned_matrix)] <- "Copepod prey length (mm)"

#PGLS of characters vs copepod prey length
pGLSp_preylength = as.data.frame(matrix(nrow=ncol(PL_pruned_matrix[,-1]), ncol=1))
colnames(pGLSp_preylength) = "Copepod prey length"
row.names(pGLSp_preylength) = names(PL_pruned_matrix[,-1])
pGLS_PL_sign = pGLSp_preylength
pGLS_sign = pGLSp_preylength
for(i in 2:(ncol(PL_pruned_matrix)-1)){
  CH_I=as.numeric(PL_pruned_matrix[,i])
  names(CH_I) = PL_pruned_matrix$Species
  CH_I = CH_I[!is.na(CH_I)]
  CH_J=PL_pruned_matrix[,ncol(PL_pruned_matrix)]
  names(CH_J) = PL_pruned_matrix$Species
  CH_I = CH_I[!is.na(CH_I)]
  CH_I = CH_I[which(names(CH_I) %in% names(CH_J))]
  CH_J = CH_J[which(names(CH_J) %in% names(CH_I))]
  SE_I = PL_pruned_SEs[which(PL_pruned_SEs$Species %in% names(CH_I)),i]
  names(SE_I) = PL_pruned_SEs[which(PL_pruned_SEs$Species %in% names(CH_I)),1]
  SE_I = SE_I[match(names(CH_I), names(SE_I))]
  pgls_tree = drop.tip(ultram, which(!(ultram$tip.label %in% names(CH_I))))
  pg_data = as.data.frame(cbind(CH_I,CH_J)) %>% .[which(!is.na(rowSums(.))),]
  pgls.SEy(CH_I~CH_J, data=pg_data, tree = pgls_tree, se = SE_I) -> TTABLE
  pgls.SEy(CH_I~sample(CH_J, length(CH_J)), data=pg_data, tree = pgls_tree, se = SE_I) -> TTABLE_null
  TTABLE <- gls(CH_I ~ CH_J, correlation =  corBrownian(phy = pgls_tree), data = pg_data, method = "ML") %>% summary() %>% .$tTable %>% as.data.frame()
  TTABLE$`p-value`[2] -> pGLSp_preylength[i-1,1]
  if(pGLSp_preylength[i-1,1]<0.05){
    print(paste(rownames(pGLSp_preylength)[i-1],"Copepod prey length",pGLSp_preylength[(i-1),1],sep=" "))
    ValueIJ = TTABLE %>% .$Value %>% .[2]
    if(ValueIJ>0){pGLS_PL_sign[i-1,1] = "+"}
    else pGLS_PL_sign[i,1] = "-"
  }
}
write.csv(pGLSp_preylength, "PGLS_Pval_preylength.csv")
pGLS_PL_sign <- pGLS_PL_sign[which(rowSums(is.na(pGLS_PL_sign))!=ncol(pGLS_PL_sign)),which(colSums(is.na(pGLS_PL_sign))!=nrow(pGLS_PL_sign))]
write.csv(pGLS_PL_sign, "PGLS_SIGN_preylength.csv")

## Soft bodied / Hard bodied prey ##
#Soft vs Hard
soft_hard = cbind((diet$`Copepod diet`+diet$`Crustacean diet` + diet$`Decapod diet`+ diet$`Amphipod diet` + diet$`Ostracod diet`), (diet$`Fish diet`+diet$`Chaetognath diet`+diet$`Gelatinous diet`+diet$`Mollusc diet`+diet$`Polychaete diet`))
colnames(soft_hard) = c("Hard", "Soft")
rownames(soft_hard) = rownames(diet)
soft_hard <- as.data.frame(soft_hard)
soft_hard[soft_hard>0]<-1
softORhard = as.data.frame(soft_hard$Hard+soft_hard$Soft)
rownames(softORhard)=rownames(soft_hard)
names(softORhard) = "Type"
softORhard[softORhard==2]<-"Both"
softORhard[which(soft_hard$Hard==1 & soft_hard$Soft==0),]<-"Hard"
softORhard[which(soft_hard$Hard==0 & soft_hard$Soft==1),]<-"Soft"

##Test literature hypotheses ##

#Nematocyst type number vs soft/hard specialization test
nematocyst_number=as.numeric(categorical$Heteroneme.type!="")[which(!is.na(categorical$Desmonemes))]+as.numeric(categorical$Haploneme.type!="")[which(!is.na(categorical$Desmonemes))]+as.numeric(categorical$Rhopalonemes!="0")[which(!is.na(categorical$Desmonemes))]+as.numeric(categorical$Desmonemes!="0")[which(!is.na(categorical$Desmonemes))]
names(nematocyst_number)=categorical$Species[which(!is.na(categorical$Desmonemes))]
CH_J=softORhard[,1]
names(CH_J) = rownames(softORhard)
CH_J[CH_J!="Soft"]<-0
CH_J[CH_J=="Soft"]<-1
CH_I = nematocyst_number
CH_J = CH_J[which(names(CH_J)%in%names(CH_I))]
CH_I = CH_I[which(names(CH_I)%in%names(CH_J))]
treeIJ = drop.tip(ultram, which(!(ultram$tip.label %in% names(CH_J))))
CH_I = CH_I[match(treeIJ$tip.label, names(CH_I))]
CH_J = CH_J[match(treeIJ$tip.label, names(CH_J))]
datIJ = data.frame(preytype = CH_J, character = CH_I)
datIJ$preytype = as.numeric(as.character(datIJ$preytype))
logit = glm(preytype~character,data=datIJ, family="binomial")
print(logit)
phylogit = phyloglm(preytype~character,phy=treeIJ,data=datIJ,boot=10, btol=10)
print(phylogit)
c("Soft specilization", "Number of nematocyst types", phylogit$n, phylogit$aic, coef(summary(phylogit))[2,6], coef(summary(phylogit))[2,1], logit$aic, coef(summary(logit))[2,4], coef(summary(logit))[2,1])

#Desmoneme/Rhopaloneme vs crustacean
dprunedCat = binary[which(rownames(binary)%in% rownames(diet)),]
CH_I = dprunedCat$Rhopalonemes #same result with desmonemes, identical phylogenetic distribution
names(CH_I) = rownames(dprunedCat)
CH_J = diet$`Crustacean diet`
names(CH_J) = rownames(diet)
CH_J = CH_J[which(names(CH_J)%in%names(CH_I))]
CH_I = CH_I[which(names(CH_I)%in%names(CH_J))]
treeIJ = drop.tip(ultram, which(!(ultram$tip.label %in% names(CH_J))))
CH_I = CH_I[match(treeIJ$tip.label, names(CH_I))]
CH_J = CH_J[match(treeIJ$tip.label, names(CH_J))]
fitPagel(treeIJ, CH_I, CH_J)

#Differentiated cnidobands with heteronemes
CH_I = dprunedCat$Proximal.heteronemes
names(CH_I) = rownames(dprunedCat)
CH_J = diet$`Crustacean diet`
names(CH_J) = rownames(diet)
CH_J = CH_J[which(names(CH_J)%in%names(CH_I))]
CH_I = CH_I[which(names(CH_I)%in%names(CH_J))]
treeIJ = drop.tip(ultram, which(!(ultram$tip.label %in% names(CH_J))))
CH_I = CH_I[match(treeIJ$tip.label, names(CH_I))]
CH_J = CH_J[match(treeIJ$tip.label, names(CH_J))]
fitPagel(treeIJ, CH_I, CH_J) %>% .$P

## OUwie regime tree of feeding guilds ##

regimetree = dprunedTree
hypdiet = c("Small crustacean", "Small crustacean", "Small crustacean", "Small crustacean", "Small crustacean", "Large crustacean", "Mixed", "Mixed", "Mixed", "Large crustacean", "Large crustacean", "Mixed", "Small crustacean", "Large crustacean", "Fish", "Fish", "Fish", "Large crustacean", "Gelatinous", "Fish", "Fish", "Fish")
names(hypdiet) = regimetree$tip.label
DietAnc = ace(hypdiet, regimetree, type = "discrete")$lik.anc
HypDietAnc = 1:nrow(DietAnc)
for(row in 1:nrow(DietAnc)){
  SCR = scale(DietAnc[row,])
  print(SCR)
  HypDietAnc[row] <- rownames(SCR)[which(SCR[,1] == max(SCR[,1]))] %>% print()
}
regimetree$node.label = as.factor(HypDietAnc)
plotTree(regimetree)
nodelabels(text=regimetree$node.label,frame="none",adj=c(1.6,-0.45), cex=0.6);tiplabels(text=hypdiet, frame="none", cex=0.6, adj=c(0,2))
RTorder = regimetree$tip.label
diet[match(RTorder, rownames(diet)),] %>% as.matrix() %>% heatmap(Rowv=NA, Colv=NA, col=c("grey", "black"))

#PGLS of characters vs Purcell Selectivity
pGLSp_sel = as.data.frame(matrix(nrow=ncol(Sprunedmatrix_logs[,-1]), ncol=ncol(selectivity)))
pGLS__sel_sign = pGLSp_sel
colnames(pGLSp_sel) = names(selectivity)
row.names(pGLSp_sel) = names(Sprunedmatrix_logs[,-1])
pGLS_sel_sign = pGLSp_sel
for(i in 2:ncol(Sprunedmatrix_logs)){
  CH_I=as.numeric(Sprunedmatrix_logs[,i])
  names(CH_I) = Sprunedmatrix_logs$Species
  CH_I = CH_I[!is.na(CH_I)]
  for(j in 1:ncol(selectivity)){
    CH_J=as.numeric(selectivity[,j])
    names(CH_J) = rownames(selectivity)
    CH_J = CH_J[!is.na(CH_J)]
    CH_I = CH_I[which(names(CH_I) %in% names(CH_J))]
    CH_J = CH_J[which(names(CH_J) %in% names(CH_I))]
    SE_I = castSE_logs[which(castSE$Species %in% names(CH_I)),i]
    names(SE_I) = castSE_logs[which(castSE$Species %in% names(CH_I)),1]
    SE_I = SE_I[match(names(CH_I), names(SE_I))]
    pgls_tree = drop.tip(Sprunedtree, which(!(Sprunedtree$tip.label %in% names(CH_I))))
    pg_data = as.data.frame(cbind(CH_I,CH_J)) %>% .[which(!is.na(rowSums(.))),]
    if(length(unique(CH_J))>1){
      pgls.SEy(CH_I~CH_J, data=pg_data, tree = pgls_tree, se = SE_I) -> TTABLE
      pgls.SEy(CH_I~sample(CH_J, length(CH_J)), data=pg_data, tree = pgls_tree, se = SE_I) -> TTABLE_null
      LikRat <- 2*(TTABLE$logLik[1] - TTABLE_null$logLik[1])
      pchisq(LikRat, df=1, lower.tail = F) -> pGLSp_sel[i-1,j]
      if(pGLSp_sel[i-1,j]<0.05){
        print(paste(rownames(pGLSp_sel)[i-1],colnames(pGLSp_sel)[j],pGLSp_sel[(i-1),(j)],TTABLE$coefficients[2], sep=" "))
        ValueIJ = TTABLE %>% .$coefficients %>% .[2]
        if(ValueIJ>0){pGLS_sel_sign[i-1,j] = "+"}
        else pGLS_sel_sign[i,j] = "-"
      }
    }
    else TTABLE<-NA; TTABLE_null<-NA; LikRat<-NA; pGLSp_sel[i-1,j]<-NA
  }
}
write.csv(pGLSp_sel, "PGLS_Pval_Selectivity.csv")
pGLS_sel_sign <- pGLS_sel_sign[which(rowSums(is.na(pGLSp_sel))!=ncol(pGLS_sel_sign)),which(colSums(is.na(pGLS_sel_sign))!=nrow(pGLS_sel_sign))]
write.csv(pGLS_sel_sign, "PGLS_SIGN_Selectivity.csv")

## Logistic regressions ##
logits = data.frame(matrix(ncol = 9, nrow = 0))
names(logits) = c("Prey type", "Character", "Ntaxa", "phyloglmAIC", "phyloglmP","phyloglm_b","glmAIC", "glmP", "glm_b")
logit_diet = diet[,-7]
logits[1,] = rep(1,9)
for(prey in 1:ncol(logit_diet)){
  for(char in 2:ncol(dprunedmatrix_logs)){
    CH_J=logit_diet[,prey]
    names(CH_J) = rownames(logit_diet)
    CH_I = dprunedmatrix_logs[,char]
    names(CH_I) = dprunedmatrix_logs$Species
    CH_I = CH_I[!is.na(CH_I)]
    CH_J = CH_J[which(names(CH_J)%in%names(CH_I))]
    treeIJ = drop.tip(dprunedTree, which(!(dprunedTree$tip.label %in% names(CH_J))))
    CH_I = CH_I[match(treeIJ$tip.label, names(CH_I))]
    CH_J = CH_J[match(treeIJ$tip.label, names(CH_J))]
    datIJ = data.frame(preytype = CH_J, character = CH_I)
    if(length(unique(CH_I))>1 & length(unique(CH_J))>1){
      logit = glm(preytype~character,data=datIJ, family="binomial")
      print(logit)
      phylogit = phyloglm(preytype~character,phy=treeIJ,data=datIJ,boot=10, btol=30)
      print(phylogit)
      rowIJ = c(names(logit_diet)[prey], names(dprunedmatrix_logs)[char], phylogit$n, phylogit$aic, coef(summary(phylogit))[2,6], coef(summary(phylogit))[2,1], logit$aic, coef(summary(logit))[2,4], coef(summary(logit))[2,1])
    } else{
      rowIJ = c(names(logit_diet)[prey], names(dprunedmatrix_logs)[char], NA, NA, NA, NA, NA, NA, NA)
    }
    print(rowIJ)
    logits = rbind(logits, rowIJ)
  }
}
logits = logits[-1,]
logits[,4:9] = sapply(logits[,4:9], as.numeric)
logits[,4:9] = sapply(logits[,4:9], function(x){round(x,3)})
write.csv(logits, "diet_morph_logits.csv", row.names = F)

## OUWIE ##
OUwie_matrix = as.data.frame(matrix(ncol=5, nrow=nrow(dprunedmatrix_logs)))
names(OUwie_matrix) = c("Character", "Ntips", "dAICc_BM", "dAICc_OU1", "dAICc_OUm")
library(magrittr)
for(c in 2:ncol(dprunedmatrix_logs)){
  Cdata = data.frame(dprunedmatrix$Species,hypdiet,dprunedmatrix_logs[,c])
  names(Cdata) = c("species", "regime", "trait")
  Cdata$trait[which(Cdata$trait == "NaN")] <- NA
  Cdata = Cdata[!is.na(Cdata$trait),]
  Cdata$trait = as.numeric(Cdata$trait)
  Cdata$species = as.character(Cdata$species)
  Cdata$regime = as.character(Cdata$regime)
  Ctree = drop.tip(regimetree,which(!(regimetree$tip.label %in% Cdata$species)))
  BM <- OUwie(Ctree, Cdata, model="BM1")
  OU1 <- OUwie(Ctree, Cdata, model="OU1")
  OUm <- OUwie(Ctree, Cdata, model="OUMV")
  OUwie_matrix[c-1,1] = names(dprunedmatrix_logs)[c]
  OUwie_matrix[c-1,2] = nrow(Cdata)
  OUwie_matrix[c-1,3] = BM$AICc - min(c(BM$AICc, OU1$AICc, OUm$AICc))
  OUwie_matrix[c-1,4] = OU1$AICc - min(c(BM$AICc, OU1$AICc, OUm$AICc))
  OUwie_matrix[c-1,5] = OUm$AICc - min(c(BM$AICc, OU1$AICc, OUm$AICc))
}
OUwie_matrix[,3:5] <- apply(OUwie_matrix[,3:5], 2, round3)
write.csv(OUwie_matrix, "OU_diet_test.csv")

## DAPC ##
#General prediction of feeding guild (transforming NAs to zeroes)
ldamtrix = sharedmean_logs
ldamtrix[is.na(ldamtrix)]<-0
ldamtrix$Species = as.character(ldamtrix$Species)
ldamtrix = ldamtrix[which(ldamtrix$Species %in% names(hypdiet)),]
prunediets = hypdiet[which(names(hypdiet) %in% ldamtrix$Species)]
#prunediets["Hippopodius hippopus"]<-"Ostracod"
ldamtrix = cbind(ldamtrix, prunediets[match(ldamtrix$Species, names(prunediets))])
names(ldamtrix)[ncol(ldamtrix)]<-"Diet"
xval <- xvalDapc(ldamtrix[,c(-1,-ncol(ldamtrix))], ldamtrix$Diet, n.pca.max = 300, training.set = 1,
                 result = "groupMean", center = TRUE, scale = TRUE,
                 n.pca = NULL, n.rep = 30, xval.plot = TRUE)
dapc(ldamtrix[,c(-1,-ncol(ldamtrix))], grp=ldamtrix$Diet, n.pca = xval$DAPC$n.pca, n.da = xval$DAPC$n.da) -> DAPClogs
contrib <- loadingplot(DAPClogs$var.contr, axis=1, lab.jitter=1)
scatter(DAPClogs, clabel = 0.5)
assignplot(DAPClogs)
compoplot(DAPClogs, posi="bottomright", lab="",ncol=1, xlab="individuals", cex.leg = 0.1, col.pal = "funky")
predictionset <- castmean_logs[,which(names(castmean_logs) %in% names(ldamtrix))]
predictionset[is.na(predictionset)]<-0
rownames(predictionset)<-predictionset$Species
predictionset=predictionset[which(!(predictionset$Species %in% c("Thermopalia taraxaca", "Nectadamas richardi", "Halistemma cupulifera", "Cardianecta parchelion"))),] #Remove taxa with many unmeasured NAs that create zero-biases
predictionset=predictionset[which(!(predictionset$Species %in% ldamtrix$Species)),]
preDIET <- predict(DAPClogs, predictionset[,-1])
cbind(predictionset$Species, as.character(preDIET$assign)) %>% View()
preDIET$posterior %>% round(5) %>%  as.matrix() %>% heatmap(scale="row", cexCol=0.5, col=c("white","white",gray.colors(10)[10:1],"black"),Colv=NA)

#Particular diet items that have enough variability (Copepods, fish, large crustaceans)
prunediets = diet[which(rownames(diet) %in% ldamtrix$Species),]

#Copepods
ldamtrix = sharedmean_logs
ldamtrix[is.na(ldamtrix)]<-0
ldamtrix$Species = as.character(ldamtrix$Species)
ldamtrix = ldamtrix[which(ldamtrix$Species %in% rownames(diet)),]
ldamtrix = cbind(ldamtrix, prunediets[match(ldamtrix$Species, rownames(prunediets)), "Copepod diet"])
names(ldamtrix)[ncol(ldamtrix)]<-"CopepodDiet"
xval <- xvalDapc(ldamtrix[,c(-1,-ncol(ldamtrix))], ldamtrix$CopepodDiet, n.pca.max = 300, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = TRUE,
                 n.pca = NULL, n.rep = 30, xval.plot = TRUE)
dapc(ldamtrix[,c(-1,-ncol(ldamtrix))], grp=ldamtrix$CopepodDiet, n.pca = xval$DAPC$n.pca, n.da = xval$DAPC$n.da) -> DAPClogs
summary.dapc(DAPClogs)
contrib <- loadingplot(DAPClogs$var.contr, axis=1, lab.jitter=1)
scatter(DAPClogs, clabel = 0.5)
assignplot(DAPClogs)
compoplot(DAPClogs, posi="bottomright", lab="",ncol=1, xlab="individuals", cex.leg = 0.1, col.pal = "funky")
copselmatrix = cbind(Sprunedmatrix, selectivity$Selectivity.Copepods)
copselmatrix = copselmatrix[which(copselmatrix$Species %in% LDAtree$tip.label),]
copselmatrix = copselmatrix[match(LDAtree$tip.label, copselmatrix$Species),]
CopepodGLM <- glm(CopepodDiet~Tentacle.width..um.+Haploneme.row.number..um., data = ldamtrix)
summary(CopepodGLM)
1-(CopepodGLM$deviance/CopepodGLM$null.deviance) #76.2
CopepodSelGLM <- glm(selectivity$Selectivity.Copepods~Tentacle.width..um.+Haploneme.row.number..um.+total_haploneme_volume, data = copselmatrix)
summary(CopepodSelGLM)
1-(CopepodSelGLM$deviance/CopepodSelGLM$null.deviance) #76.6%
rownames(ldamtrix) = ldamtrix$Species
CopepodPGLM <- phyloglm(CopepodDiet~Tentacle.width..um.+Haploneme.row.number..um., data = ldamtrix, phy = LDAtree, btol=100)
summary(CopepodPGLM)
predictionset <- castmean_logs[,which(names(castmean_logs) %in% names(ldamtrix))]
predictionset[is.na(predictionset)]<-0
rownames(predictionset)<-predictionset$Species
predictionset=predictionset[which(!(predictionset$Species %in% c("Thermopalia taraxaca", "Nectadamas richardi", "Halistemma cupulifera", "Cardianecta parchelion"))),] #Remove taxa with many unmeasured NAs that create zero-biases
predictionset=predictionset[which(!(predictionset$Species %in% ldamtrix$Species)),]
preDIET <- predict(DAPClogs, predictionset[,-1])
cbind(predictionset$Species, as.character(preDIET$assign)) %>% View()
preDIET$posterior %>% round(5) %>%  as.matrix() %>% heatmap(scale="row", cexCol=0.5, col=c("white","white",gray.colors(10)[10:1],"black"),Colv=NA)

#Fish
ldamtrix = sharedmean_logs
ldamtrix[is.na(ldamtrix)]<-0
ldamtrix$Species = as.character(ldamtrix$Species)
ldamtrix = ldamtrix[which(ldamtrix$Species %in% rownames(diet)),]
ldamtrix = cbind(ldamtrix, prunediets[match(ldamtrix$Species, rownames(prunediets)), "Fish diet"])
names(ldamtrix)[ncol(ldamtrix)]<-"FishDiet"
xval <- xvalDapc(ldamtrix[,c(-1,-ncol(ldamtrix))], ldamtrix$FishDiet, n.pca.max = 300, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = TRUE,
                 n.pca = NULL, n.rep = 30, xval.plot = TRUE)
dapc(ldamtrix[,c(-1,-ncol(ldamtrix))], grp=ldamtrix$FishDiet, n.pca = xval$DAPC$n.pca, n.da = xval$DAPC$n.da) -> DAPClogs
summary.dapc(DAPClogs)
contrib <- loadingplot(DAPClogs$var.contr, axis=1, lab.jitter=1)
scatter(DAPClogs, clabel = 0.5)
assignplot(DAPClogs)
compoplot(DAPClogs, posi="bottomright", lab="",ncol=1, xlab="individuals", cex.leg = 0.1, col.pal = "funky")
FishGLM <- glm(FishDiet~Heteroneme.volume..um3.+total_haploneme_volume+Pedicle.width..um., data = ldamtrix)
summary(FishGLM)
1-(FishGLM$deviance/FishGLM$null.deviance)

#Fish selectivity GLM
fishselmatrix = cbind(Sprunedmatrix, selectivity$Selectivity.Fish)
fishselmatrix = fishselmatrix[which(fishselmatrix$Species %in% ldamtrix$Species),]
fishselmatrix = fishselmatrix[match(LDAtree$tip.label, fishselmatrix$Species),]
FishSelGLM <- glm(selectivity$Selectivity.Fish~Pedicle.width..um., data = fishselmatrix)
summary(FishSelGLM)
1-(FishSelGLM$deviance/FishSelGLM$null.deviance) #97.5%
predictionset <- castmean_logs[,which(names(castmean_logs) %in% names(ldamtrix))]
predictionset[is.na(predictionset)]<-0
rownames(predictionset)<-predictionset$Species
predictionset=predictionset[which(!(predictionset$Species %in% c("Thermopalia taraxaca", "Nectadamas richardi", "Halistemma cupulifera", "Cardianecta parchelion"))),] #Remove taxa with many unmeasured NAs that create zero-biases
predictionset=predictionset[which(!(predictionset$Species %in% ldamtrix$Species)),]
preDIET <- predict(DAPClogs, predictionset[,-1])
cbind(predictionset$Species, as.character(preDIET$assign)) %>% View()
preDIET$posterior %>% round(5) %>%  as.matrix() %>% heatmap(scale="row", cexCol=0.5, col=c("white","white",gray.colors(10)[10:1],"black"),Colv=NA)

#Large crustaceans
ldamtrix = sharedmean_logs
ldamtrix[is.na(ldamtrix)]<-0
ldamtrix$Species = as.character(ldamtrix$Species)
ldamtrix = ldamtrix[which(ldamtrix$Species %in% rownames(diet)),]
ldamtrix = cbind(ldamtrix, prunediets[match(ldamtrix$Species, rownames(prunediets)), "Decapod diet"])
names(ldamtrix)[ncol(ldamtrix)]<-"DecapodDiet"
xval <- xvalDapc(ldamtrix[,c(-1,-ncol(ldamtrix))], ldamtrix$DecapodDiet, n.pca.max = 300, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = TRUE,
                 n.pca = NULL, n.rep = 30, xval.plot = TRUE)
dapc(ldamtrix[,c(-1,-ncol(ldamtrix))], grp=ldamtrix$DecapodDiet, n.pca = xval$DAPC$n.pca, n.da = xval$DAPC$n.da) -> DAPClogs
summary.dapc(DAPClogs)
contrib <- loadingplot(DAPClogs$var.contr, axis=1, lab.jitter=1)
scatter(DAPClogs, clabel = 0.5)
assignplot(DAPClogs)
compoplot(DAPClogs, posi="bottomright", lab="",ncol=1, xlab="individuals", cex.leg = 0.1, col.pal = "funky")
DecapodGLM <- glm(DecapodDiet~Heteroneme.volume..um3.+coiledness+Cnidoband.length..um.+haploneme_elongation+heteroneme_elongation, data = ldamtrix)
summary(DecapodGLM)
1-(DecapodGLM$deviance/DecapodGLM$null.deviance)

#Large crustacean selectivity GLM
decselmatrix = cbind(Sprunedmatrix, selectivity$Selectivity.Decapod.larvae)
decselmatrix = decselmatrix[which(decselmatrix$Species %in% ldamtrix$Species),]
decselmatrix = decselmatrix[match(LDAtree$tip.label,decselmatrix$Species),]
DecapodSelGLM <- glm(selectivity$Selectivity.Decapod.larvae~Heteroneme.volume..um3.+coiledness+Cnidoband.length..um.+haploneme_elongation+heteroneme_elongation, data = decselmatrix)
summary(DecapodSelGLM)
1-(DecapodSelGLM$deviance/DecapodSelGLM$null.deviance) #89.1%
predictionset <- castmean_logs[,which(names(castmean_logs) %in% names(ldamtrix))]
predictionset[is.na(predictionset)]<-0
rownames(predictionset)<-predictionset$Species
predictionset=predictionset[which(!(predictionset$Species %in% c("Thermopalia taraxaca", "Nectadamas richardi", "Halistemma cupulifera", "Cardianecta parchelion"))),] #Remove taxa with many unmeasured NAs that create zero-biases
predictionset=predictionset[which(!(predictionset$Species %in% ldamtrix$Species)),]
preDIET <- predict(DAPClogs, predictionset[,-1])
cbind(predictionset$Species, as.character(preDIET$assign)) %>% View()
preDIET$posterior %>% round(5) %>%  as.matrix() %>% heatmap(scale="row", cexCol=0.5, col=c("white","white",gray.colors(10)[10:1],"black"),Colv=NA)

## SIMPLE ANALYSES OF KINEMATIC DATA ##

kine = read.csv("Kinematics.tsv", sep='\t', header=T) #%>% .[which(apply(., 1, function(x) sum(is.na(x)))<ncol(.)-2),-2]
kineWNA = kine[which(kine$Species %in% ultram$tip.label),]
rownames(kineWNA)=kineWNA$Specimen
rownames(kine)=kine$Specimen
#kineWNA = kineWNA[,-1]
#kine=kine[,-1]
kinetree = drop.tip(ultram, which(!(ultram$tip.label %in% kine$Species)))
kinetree$edge.length = 200*kinetree$edge.length
#kine_clean = kineWNA[,which(colSums(kineWNA)>1)]
#names(kine) = c("ADS", "MDS", "HeDS","HeMDS", "HSDS", "HFL", "HaDS")
kine_byspp = aggregate(. ~ kine$Species, data = kine[,c(-1,-2)], mean.na, na.action = na.pass)
names(kine_byspp)[1] <- "Species"
kinemorph <- castmeans[which(castmeans$Species %in% kine_byspp$Species),] %>% cbind(kine_byspp[which(kine_byspp$Species %in% castmeans$Species),])
plot(kinemorph$Cnidoband.free.length..um., kinemorph$Average.CB.discharge.speed..mm.s.)
calys = kinemorph$Species[c(5,9,12,15,16,17,18)]
euphys = kinemorph$Species[which(!(kinemorph$Species %in% calys))]
relspeed = data.frame(kinemorph$Species, c(kinemorph$Average.CB.discharge.speed..mm.s./kinemorph$Cnidoband.free.length..um.))
names(relspeed) = c("Species", "Speed/Length")
relspeed <- relspeed[which(!(is.na(relspeed$`Speed/Length`))),]
t.test(relspeed[which(relspeed$Species %in% calys),2], relspeed[which(relspeed$Species %in% euphys),2])
t.test(kine$Average.CB.discharge.speed..mm.s.[which(kine$Species %in% calys)], kine$Average.CB.discharge.speed..mm.s.[which(kine$Species %in% euphys)])
t.test(castnumbers$Cnidoband.free.length..um.[which(castnumbers$Species %in% calys)], castnumbers$Cnidoband.free.length..um.[which(castnumbers$Species %in% euphys)])
lm(kinemorph$Average.CB.discharge.speed..mm.s.~ kinemorph$Cnidoband.free.length..um.) %>% summary()
lm(kinemorph$Heteroneme.discharge.speed.MAX..mm.s.~ kinemorph$Heteroneme.volume..um3.) %>% summary()
t.test(kine$Heteroneme.discharge.speed.AVG..mm.s.[which(kine$Species %in% calys)], kine$Heteroneme.discharge.speed.AVG..mm.s.[which(kine$Species %in% euphys)])
t.test(kine$Haploneme.discharge.speed.AVG..mm.s.[which(kine$Species %in% calys)], kine$Haploneme.discharge.speed.AVG..mm.s.[which(kine$Species %in% euphys)])
t.test(kine$Haploneme.discharge.speed.AVG..mm.s., kine$Heteroneme.discharge.speed.AVG..mm.s.)
cor(kinemorph[,c(-1,-32)], use="pairwise.complete.obs") %>% .[c(1:30),c(31:41)] %>% corrplot(diag=F, tl.cex = 0.4, tl.col="black")
