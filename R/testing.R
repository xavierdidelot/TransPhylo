testing=function() {
rm(list=ls())
set.seed(0)
neg=100/365
off.r=5
off.p=0.5
w.shape=10
w.scale=0.1
pi=0.1
simu <- simulateOutbreak(neg=neg,pi=pi,off.r=off.r,off.p=off.p,w.shape=w.shape,
                         w.scale=w.scale,dateStartOutbreak=2005,datePresent=2009)
ttree<-ttreeFromFullTree(simu)
plotTTree(ttree,w.shape,w.scale,showLabels=T)

ptree<-ptreeFromFullTree(simu)
record<-inferTTree(ptree,mcmcIterations=1000,w.shape=w.shape,w.scale=w.scale,
                   startNeg=neg,startPi=0.5,startOff.p=off.p,startOff.r=1,
                   updateOff.r=T,updatePi=T,updateNeg=T,updateOff.p=F,datePresent=2010)
cons=consTTree(record[seq(1,1000,10)])
plotTTree(cons,w.shape,w.scale,showLabels=T)
}