rm(list=ls())
source('probTTree.R')
source('makeTTree.R')
source('withinhost.R')
source('glueTrees.R')
source('hostFromFulltree.R')
source('ttreeFromFullTree.R')
source('probPTreeGivenTTree.R')
source('pTreeFromFullTree.R')
source('makeFullTreeFromPTree.R')
neg<-100/365#Within-host effective population size (Ne) times  generation duration (g)
R<-1#Basic reproduction number

#Create a transmission tree with ten individuals
n<-1
while (n!=10) {
  ttree<-makeTTree(R)[[1]]
  n<-nrow(ttree)
  if (is.null(ttree)) {n<-0}
}

#Create a within-host phylogenetic tree for each infected host
wtree<-vector('list',n)
for (i in (1:n)) {
  times<-c(ttree[i,2],ttree[which(ttree[,3]==i),1])-ttree[i,1]
  wtree[[i]]<-withinhost(times,neg)[[1]];
}

#Glue these trees together
truth<-glueTrees(ttree,wtree)
truth[,1]<-truth[,1]+2005#Epidemic started in 2005
makeFullTreeFromPTree(ptreeFromFullTree(truth))
