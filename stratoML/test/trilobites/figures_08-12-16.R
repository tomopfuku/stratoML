#potential_figures_01-26-16.R

# Figure Guidelines
# Paleobiology: 1 column ( ~7 cm) or 2 column ( ~15 cm) figures only.
# Max height: ~21.5 cm
# .tiff
# Resolution of 600 dpi or higher. Line drawings should be submitted at 1200 dpi.

# source("C:\\dave\\research\\1 trilobites and cal3\\analysis 03-25-16\\figures_05-30-16.R")

library(paleotree)
library(beanplot)

wd<-paste0(getwd(),"/figures")
setwd(wd)

names(timeTrees)[names(timeTrees)=="basic"]<-"MND"
names(corrCoefListNoInf)[names(corrCoefListNoInf)=="basic"]<-"MND"
names(corrCoefListNoInfTT)[names(corrCoefListNoInfTT)=="basic"]<-"MND"
names(nHalfLife)[names(nHalfLife)=="basic"]<-"MND"

###########################################################

# sampling estimates

tiff(file="Fig1.sampEst.tif",height=4,
	width=15,units="cm",
	res=1200,compression="lzw+p")

layout(matrix(1:4,,4),widths=c(0.3,0.25,0.25,0.25))
oldPar<-par(no.readonly = T)
nbreaks=15
#
par(mar=c(5,4,1,1))
hist(-apply(timeStackPruned[[1]],1,diff),
	main="",breaks=nbreaks,
	xlab="Stratigraphic Durations\n for Species (Myr)")
#legend("topright","(A)",box.lty=0,cex=1.35)
text(x=2.5,y=7.4,
	labels="A",cex=1.7)
#
par(mar=c(5,2,1,1))
hist(sapply(likeResPruned,function(x) x$par[1]),
	main="",breaks=nbreaks,
	xlab="Extinction Rate\n (per LMyr)")
text(x=1.55,y=12,
	labels="B",cex=1.7)
#
par(mar=c(5,2,1,1))
hist(sapply(likeResPruned,function(x) x$par[2]),
	main="",ylab="",breaks=nbreaks,
	xlab="Sampling Rate\n (per LMyr)")
text(x=42,y=25,
	labels="C",cex=1.7)
#
par(mar=c(5,2,1,2))
hist(sampRatePruned,
	main="",ylab="",breaks=nbreaks,
	xlab="Minimum Sampling Rate\n (per LMyr)")
text(x=6.55,y=13,
	labels="D",cex=1.7)
layout(1)
par(oldPar)

dev.off()

##############################################################

# comparative dating

tiff(file="Fig2.nodedating.tif",height=10,
	width=15,units="cm",
	res=1200,compression="lzw+p")

oldPar<-par(no.readonly = T)
par(mar=c(6,4.5,0.5,0.5))
boxplot(rootAges,names=names(timeTrees),cex.lab=1.2,
	cex.axis=1.1,
      ylab="Age (Ma)",
	las=2,ylim=c(498.7,494.5))
boxplot(SigPteroAges,names=names(timeTrees),add=TRUE,
	las=2,border="gray50",cex.axis=1.1)
legend1<-expression(paste("MRCA of ",
		italic("Sigmocheilus"), 
		"\n & " ,
		italic("Pterocephalia")))
legend("bottomleft",
	pch=c(15,15),col=c("gray60",1),pt.cex=2,cex=.9,
	legend=c(legend1,
		"Root Divergence")
	)
par(oldPar)

dev.off()

#############################################################

# Tau Correlation of Ancestor-Descendant Distances
	# only no ZLB analyses

tiff(file="Fig3.ADdistRes.tif",height=12,
	width=15,units="cm",
	res=1200,compression="lzw+p")

layout(1:2,heights=c(0.4,0.6))
oldPar<-par(no.readonly = T)
par(mar=c(0,5,0.5,0.5),xaxt="n")
#
ylims<-c(-0.5,0.6)
#
# ASR w/Unit-Length Tree
boxplot(corrCoefListNoInf,ylab="",las=2,
         main="",cex.lab=1.3,ylim=ylims,cex.axis=1.15)
mtext(side=2,line=3.5,"Kendall's Tau",cex=1.3)
abline(h=0,lty=2,lwd=1.5)
text(x=9.9,y=-0.45,pos=2,
	labels="Ancestral Values via Unit-Length Tree",cex=1.15)
text(x=9.9,y=0.52,pos=2,
	labels="A",cex=1.7)
#
par(mar=c(6.2,5,0,0.5),xaxt="s")
#
# ASR w/TimeTrees
boxplot(corrCoefListNoInfTT,ylab="",las=2,
         main="",cex.lab=1.3,ylim=ylims,cex.axis=1.15)
mtext(side=2,line=3.5,"Kendall's Tau",cex=1.3)
abline(h=0,lty=2,lwd=1.5)
text(x=9.9,y=-0.45,pos=2,
	labels="Ancestral Values via Time-Scaled Trees",cex=1.15)
text(x=9.9,y=0.52,pos=2,
	labels="B",cex=1.7)
#
layout(1)
par(oldPar)

dev.off()

########################################

# mvMORPH analyses

tiff(file="Fig4.mvMorphRes.tif",height=9,
	width=15,units="cm",
	res=1200,compression="lzw+p")

oldPar<-par(no.readonly = T)
par(mar=c(6.2,6.5,0.5,0.5))
boxplot(nHalfLife,las=2,cex.axis=1.1,
         ylab="")
mtext(side=2,line=2.5,cex=1.2,
	text="Expected Number of Rate\n Doublings (or Halvings)\n over 5 Myr for PCs 1 - 4")
abline(h=0,lty=2,lwd=1.5)
par(oldPar)

dev.off()

###########################################################

# Ancestor-Descendant Pairs

tiff(file="Fig5.ADpairs.tif",height=16,
	width=15,units="cm",
	res=1000,compression="lzw+p")

oldPar<-par(no.readonly = T)
par(mar=c(5,10,0.5,0.5),lheight=0.8)
namesPairs<-paste0(rownames(tb)," ")
namesPairs<-gsub("_"," ",namesPairs)
barX<-barplot(t(tb[,1:2])/100,
	names.arg=namesPairs,
	horiz=T,
	las=1,xlim=c(0,1.1),
	font=3,xaxt="n",
	cex.names=0.70)
#add x axis
axis(1)
# a-priori markers
points(y=barX[g1[,2]],x=g1[,1]+0.05,
     pch=16,cex=2,col='gray20')
points(y=barX[g2[,2]],x=g2[,1]+0.05,
     pch=16,cex=2,col='gray75')
# x-axis label
mtext(side=1,line=3.3,
	text="Proportion of cal3 Runs Placed\n as Ancestor-Descendant Pair")
# legend
legend("bottomright",inset=c(0,0),bty='n',
       legend=c('Anagenetic AD Pair','Budding AD Pair',
		'a-priori: Anagenetic AD Pair','a-priori: Budding AD Pair'),
       pch=c(15,15,16,16),pt.cex=c(2,2,1.9,1.9),
	 col=c('gray20','gray75','gray20','gray75'),cex=0.75)
par(oldPar)

dev.off()