context("Typical session")

test_that("Typical session can be run.", {
  set.seed(0)
  neg=100/365
  off.r=5
  w.shape=10
  w.scale=0.1
  pi=0.25
  dateT=2008
    expect_silent(simu <- simulateOutbreak(neg=neg,pi=pi,off.r=off.r,w.shape=w.shape,
                             w.scale=w.scale,dateStartOutbreak=2005,dateT=dateT))
   expect_is(simu,'list')
   expect_silent(plotCTree(simu))
   expect_silent(ttree<-extractTTree(simu))
   expect_is(ttree,'list')
   expect_silent(plotTTree(ttree,w.shape,w.scale))
   expect_silent(ptree<-extractPTree(simu))
   expect_is(ptree,'list')
   expect_silent(p<-phyloFromPTree(ptree))
   expect_is(p,'phylo')
   expect_silent(ptree<-ptreeFromPhylo(p,dateLastSample=max(simu$ctree[,1])))
   expect_is(ptree,'list')
   expect_is(capture_output(record<-inferTTree(ptree,mcmcIterations=100,w.shape=w.shape,w.scale=w.scale,dateT=dateT)),'character')
   expect_is(record,'list')
   lastIteration<-record[[length(record)]]
   expect_silent(plotCTree(lastIteration$ctree))
   expect_silent(plotTraces(record))
   expect_silent(mcmc<-convertToCoda(record))
   expect_is(mcmc,'mcmc')
   expect_silent(cons<-consTTree(record))
   expect_is(cons,'list')
   expect_silent(sel<-selectTTree(record))
   expect_is(sel,'numeric')
   expect_silent(plotTTree(cons,w.shape,w.scale))
   expect_silent(plotTTree2(cons,w.shape,w.scale))
   expect_silent(mat<-computeMatWIW(record))
   expect_is(mat,'matrix')
   expect_silent(mat<-computeMatTDist(record))
   expect_is(mat,'matrix')
   expect_silent(a<-getIncidentCases(record,show.plot = T))
   expect_is(a,'list')
   expect_silent(a<-getGenerationTimeDist(record,show.plot = T))
   expect_is(a,'matrix')
   expect_silent(a<-getSamplingTimeDist(record,show.plot = T))
   expect_is(a,'matrix')
   expect_silent(a<-getInfectionTimeDist(record,k='1',show.plot = T))
   expect_is(a,'numeric')
})


