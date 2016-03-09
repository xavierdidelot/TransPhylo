rm(list=ls())
library('epiphylo')
set.seed(0)
neg=100/365
off.r=2
off.p=0.5
w.shape=5
w.scale=0.2
pi=0.5
#xs=seq(0,3,0.01);plot(xs,dgamma(xs,shape=w.shape,scale=w.scale),type='l')
ttree=matrix(NA,2,2)
while (nrow(ttree)!=8){
simu <- simulateOutbreak(neg=neg,pi=pi,off.r=off.r,off.p=off.p,w.shape=w.shape,
                         w.scale=w.scale,dateStartOutbreak=2005,datePresent=2010)
ttree<-ttreeFromFullTree(simu)
}
pdf('fig1.pdf',4,6)
#Plot transmission tree
par(fig=c(0.1,0.9,0.6,1),mar=c(0,0,0,0))
plotTTree(ttree,w.shape,w.scale,showLabels=F,maxTime=2010)
cols=rainbow(nrow(ttree)+1);cols=cols[-4]
for (i in 1:nrow(ttree)) text(2010.1,9-i,i,col=cols[i])
par(xpd=NA)
text(2004.5,8,'A',cex=1.5)

#Plot colored tree
par(fig=c(0.14,0.86,0.1,0.5),new=TRUE)
reorder=c(4,5,6,1,2,3,7,8)
plotBothTree(simu,showLabels=F,cols=cols[reorder],maxTime=2010)
ys=c(4,2,1,5,4,4.5,2,3)+0.3
for (i in 1:nrow(ttree)) text(ttree[i,1]+0.15,ys[i],reorder[i],col=cols[reorder[i]])
par(xpd=NA)
text(2004.5,5,'B',cex=1.5)
dev.off()
system('open fig1.pdf')