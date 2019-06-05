# load libraries
library(BAMMtools)
library(geiger)
library(coda)
#setwd("~/tentilla_morph/Supplementary_Materials/Dryad/BAMM")
# load specific character subtree made in the R_Code ## BAMM ## block

tree <- read.tree("haploneme_elongation_tree.tre")

# ESS after 20% burnin
mcmcout <- read.csv("mcmc_out.txt")
burnstart <- floor(0.2*nrow(mcmcout))
postburn <- mcmcout[burnstart:nrow(mcmcout), ]
effectiveSize(postburn$N_shifts)
effectiveSize(postburn$logLik)

# load BAMM main output
edata <- getEventData(tree, eventdata = "event_data.txt", burnin=0.2, type="trait")
summary(edata)

# the 6 most likely scenario's within the 95% credible scenarios
# F value is the relative probability, F values from all scenarios add up to F=1
pdf("BAMM_haploneme_elongation_css.pdf", width=7, height=10)
css <- credibleShiftSet(edata, expectedNumberOfShifts = 1, threshold=5, set.limit=0.95)
BAMMtools:::plot.credibleshiftset(css, plotmax = 6, pal = "temperature", legend=T, use.plot.bammdata = T)
dev.off()

# net diversification graph
pdf("BAMM_haploneme_elongation_netdiv.pdf", width=7, height=5)
starttime <- max(branching.times(tree))
plotRateThroughTime(edata, intervalCol="darkgreen", avgCol="darkgreen", start.time=starttime, ylim=c(0,25))
dev.off()

pdf("haplooneme_elongation_ratesplot.pdf",width=16,height=16)
plot.bammdata(edata,lwd=2,legend=T,labels=T,cex=0.9,logcolor=T,breaksmethod="jenks",pal=c("#0D0887", "#CC4678", "#F0F921"))
dev.off()

#q()