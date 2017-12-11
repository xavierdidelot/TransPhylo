context("Test coalescent functions")

test_that("Coalescent probability function gives expected result on small example.", {
  a=TransPhylo:::.probSubtree(cbind(c(2010+1e-10,2010+2e-10,2010,2009,2008,0),c(4,4,5,5,6,0)),2)
  b=dexp(1,3/2,T)+dexp(1,1/2,T)+log(1/3)
  expect_equal(a,b)
})

test_that("Probabilities of a coalescent tree are the same when simulating and evaluating.",{
 set.seed(0)
 leaves=2012:2000
 tree=TransPhylo:::.withinhost(leaves,1.1)
 tre=tree$nodes
 fathers=rep(0,nrow(tre))
 for (i in (length(leaves)+1):(nrow(tre)-1)) for (j in 2:3) fathers[tre[i,j]]=i 
 fathers[tre[nrow(tre),2]]=nrow(tre)
 MySort <- sort(tre[,1],index.return = TRUE,decreasing = T);ind<-MySort$ix
 revind=ind;revind[ind]=1:length(ind)
 subtree=cbind(tre[ind,1],c(revind[fathers[ind[1:(length(ind)-1)]]],0))
 p=TransPhylo:::.probSubtree(subtree,1.1)
 expect_equal(p,tree$prob)
})

test_that("Probabilities of coalescence in an outbreak is same when simulation and evaluating.",{
  set.seed(1)
  n=5
  neg=1.1
  ttree=rbind(c(0,2000,0),c(2001,2002,1),c(2001.1,2002.1,1),c(2001.2,2002.2,1),c(2001.3,2002.3,1))
  wtree<-vector('list',n)
  p1=0
  for (i in (1:n)) {
    if (is.na(ttree[i,2])) {times<-c(           ttree[which(ttree[,3]==i),1])-ttree[i,1]}
    else {times<-c(ttree[i,2],ttree[which(ttree[,3]==i),1])-ttree[i,1]}
    x=TransPhylo:::.withinhost(times,neg)
    wtree[[i]]=x[[1]]
    p1=p1+x$prob
  }
  ft<-TransPhylo:::.glueTrees(ttree,wtree)
  ft=list(ctree=ft,nam=as.character(1:n))
  expect_equal(p1,probPTreeGivenTTree(ft,1.1))
})
