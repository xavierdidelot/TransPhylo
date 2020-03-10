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
   expect_is(capture_output(res<-inferTTree(ptree,mcmcIterations=100,w.shape=w.shape,w.scale=w.scale,dateT=dateT)),'character')
   expect_is(res,'resTransPhylo')
   expect_is(capture_output(print(res)),'character')
   expect_is(capture_output(summary(res)),'character')
   lastIteration<-res[[length(res)]]
   expect_silent(plotCTree(lastIteration$ctree))
   expect_silent(plot(res))
   expect_silent(mcmc<-as.mcmc.resTransPhylo(res))
   expect_is(mcmc,'mcmc')
   expect_silent(cons<-consTTree(res))
   expect_is(cons,'list')
   expect_silent(sel<-selectTTree(res))
   expect_is(sel,'numeric')
   expect_silent(plotTTree(cons,w.shape,w.scale))
   expect_silent(plotTTree2(cons,w.shape,w.scale))
   expect_silent(mat<-computeMatWIW(res))
   expect_is(mat,'matrix')
   expect_silent(mat<-computeMatTDist(res))
   expect_is(mat,'matrix')
   expect_silent(a<-getIncidentCases(res,show.plot = T))
   expect_is(a,'list')
   expect_silent(a<-getGenerationTimeDist(res,show.plot = T))
   expect_is(a,'matrix')
   expect_silent(a<-getSamplingTimeDist(res,show.plot = T))
   expect_is(a,'matrix')
   expect_silent(a<-getInfectionTimeDist(res,k='1',show.plot = T))
   expect_is(a,'numeric')
})

test_that("Inference with multiple trees runs", {
  set.seed(0)
  neg=100/365
  off.r=5
  w.shape=10
  w.scale=0.1
  pi=0.25
  dateT=2008
  simu=simulateOutbreak(neg=neg,pi=pi,off.r=off.r,w.shape=w.shape,
  w.scale=w.scale,dateStartOutbreak=2005,dateT=dateT)
  ptree=extractPTree(simu)
  expect_is(capture.output(rec<-infer_multittree_share_param(list(ptree),mcmcIterations=100,w.shape=w.shape,w.scale=w.scale,dateT=dateT)),'character')
  expect_is(rec,'list')
})



