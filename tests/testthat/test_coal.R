context("Test coalescent functions")

test_that("Coalescent probability function gives expected result on small ultrametric example.", {
  neg=2.1
  a=TransPhylo:::probSubtree(cbind(c(2010+1e-10,2010+2e-10,2010,2009,2008,0),c(4,4,5,5,6,0)),neg)
  b=dexp(1,choose(3,2)/neg,T)+dexp(1,choose(2,2)/neg,T)+log(1/3)
  c=TransPhylo:::coalescent(c(2010+1e-10,2010+2e-10,2010),c(2009,2008),neg)
  d=TransPhylo:::probSubtreeR(cbind(c(2010+1e-10,2010+2e-10,2010,2009,2008,0),c(4,4,5,5,6,0)),neg)
  expect_equal(a,b)
  expect_equal(a,c)
  expect_equal(a,d)
})

test_that("Coalescent probability function gives expected result on small non-ultrametric example.", {
  neg=2.1
  a=TransPhylo:::probSubtree(cbind(c(2012,2011,2010,2009,2008,0),c(4,4,5,5,6,0)),neg)
  b=TransPhylo:::probSubtreeR(cbind(c(2012,2011,2010,2009,2008,0),c(4,4,5,5,6,0)),neg)
  c=TransPhylo:::coalescent(c(2012,2011,2010),c(2009,2008),neg)
  expect_equal(a,b)
  expect_equal(a,c)
})

test_that("Probabilities of a coalescent tree are the same when simulating and evaluating.",{
  set.seed(0)
  leaves=13:1
  tree=TransPhylo:::withinhost(leaves,1.1)
  tre=tree$nodes
  fathers=rep(0,nrow(tre))
  for (i in (length(leaves)+1):(nrow(tre)-1)) for (j in 2:3) fathers[tre[i,j]]=i 
  fathers[tre[nrow(tre),2]]=nrow(tre)
  ind=order(tre[,1],decreasing=T)
  revind=ind;revind[ind]=1:length(ind)
  subtree=cbind(tre[ind,1],c(revind[fathers[ind[1:(length(ind)-1)]]],0))
  p=TransPhylo:::probSubtree(subtree,1.1)
  p2=TransPhylo:::probSubtreeR(subtree,1.1)
  expect_equal(p,tree$prob)
  expect_equal(p2,tree$prob)
})

test_that("Probabilities of coalescence in an outbreak is same when simulation and evaluating.",{
  set.seed(1)
  n=5
  neg=1.1
  ttree=rbind(c(0,2000,0),c(2001,2002,1),c(2001.1,2002.1,1),c(2001.2,2002.2,1),c(2001.3,2002.3,1))
  wtree<-vector('list',n)
  p1=0
  for (i in (1:n)) {
    if (is.na(ttree[i,2])) {times<-c(ttree[which(ttree[,3]==i),1])-ttree[i,1]}
    else {times<-c(ttree[i,2],ttree[which(ttree[,3]==i),1])-ttree[i,1]}
    x=TransPhylo:::withinhost(times,neg)
    wtree[[i]]=x[[1]]
    p1=p1+x$prob
  }
  ft<-TransPhylo:::.glueTrees(ttree,wtree)
  ft=list(ctree=ft,nam=as.character(1:n))
  expect_equal(p1,probPTreeGivenTTree(ft$ctree,1.1))
})

test_that("Probabilities of coalescence in an outbreak is same when simulation and evaluating.",{
  set.seed(4)
  neg=1.1
  ttree=makeTTree(off.r = 1,off.p =0.5,pi=1,w.shape=1,w.scale=1,ws.shape=1,ws.scale=1,maxTime=Inf,nSampled = NA)
  ttree=ttree$ttree
  n=nrow(ttree)#with this seed we get a big outbreak with n=95
  wtree<-vector('list',n)
  p1=0
  for (i in (1:n)) {
    if (is.na(ttree[i,2])) {times<-c(ttree[which(ttree[,3]==i),1])-ttree[i,1]}
    else {times<-c(ttree[i,2],ttree[which(ttree[,3]==i),1])-ttree[i,1]}
    x=TransPhylo:::withinhost(times,neg)
    wtree[[i]]=x[[1]]
    p1=p1+x$prob
  }
  ft<-TransPhylo:::.glueTrees(ttree,wtree)
  ft=list(ctree=ft,nam=as.character(1:n))
  expect_equal(p1,probPTreeGivenTTree(ft$ctree,neg))
})

test_that("Probabilities of an outbreak is same when simulation and evaluating.",{
  set.seed(4)
  neg=1.1
  shape=1.1
  scale=1.2
  out=simulateOutbreak(off.r=1,neg=neg,pi = 1,w.shape=shape,w.scale=scale,dateStartOutbreak = 2000,dateT = Inf)
  expect_equal(out$probwithin,probPTreeGivenTTree(out$ctree,neg))
  expect_equal(out$probttree,probTTree(extractTTree(out)$ttree,rOff = 1,pOff = 0.5,pi = 1,shGen = shape,scGen = scale,shSam = shape,scSam = scale,dateT=Inf),tolerance=0.001,scale=1)
})


test_that("Probabilities probPTreeGivenTTree are the same in R and Rcpp.",{
  set.seed(4)
  neg=1.1
  shape=1.1
  scale=1.2
  out=simulateOutbreak(off.r=2,neg=neg,pi = 1,w.shape=shape,w.scale=scale,dateStartOutbreak = 2000,dateT = 2005)
  expect_equal(probPTreeGivenTTreeR(out$ctree,neg),probPTreeGivenTTree(out$ctree,neg))
})


test_that("Linear coalescent probability function gives expected result on small example.", {
  rate=2.2
  coaltime=1 #must be between 0 and 2
  a=TransPhylo:::probSubtreeLinear(cbind(c(3,2,coaltime,0),c(3,3,4,0)),rate)
  b=1/(rate*coaltime)*exp(-1/rate*(log(rate*2)-log(rate*coaltime)))
  b=log(b)
  expect_equal(a,b)
})

test_that("Linear coalescent probability function gives expected result on tree with rescaled time.", {
  rate=1.1#rate of growth
  pres=10
  coal=pres-1.2
  tree=cbind(c(pres+1e-10,pres,coal,0),c(3,3,4,0))
  b=dexp(pres-coal,choose(2,2),T)#probability of unscaled tree
  b=b+log(exp(rate*(pres-coal)))
  tree2=tree
  tree2[3,1]=pres-(1-exp(-(pres-coal)*rate))/rate
  tree2[4,1]=pres-1/rate
  a=TransPhylo:::probSubtreeLinear(tree2,rate)
  expect_equal(a,b)
})

test_that("Can calculate probability of a within-host linear tree.",{
  rate=1.1
  expect_silent(tree<-TransPhylo:::withinhostLinear(10:1,1/rate)[[1]])
  expect_silent(a<-TransPhylo:::probSubtreeLinear(tree,rate))
  expect_is(a,'numeric')
})