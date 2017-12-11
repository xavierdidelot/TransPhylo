context("Test function calculating probability of transmission tree")

test_that("Probability of a transmission tree is the same when simulating and calculating.", {
  set.seed(0)
  pi=1
  off.r=1
  off.p=0.5
  w.shape=1.1
  w.scale=1.1
  ttree=makeTTree(pi=pi,off.r=off.r,off.p=off.p,w.shape=w.shape,w.scale=w.scale,ws.shape=w.shape,ws.scale=w.scale,maxTime=Inf,nSampled = NA)
  p=TransPhylo:::probTTree(ttree = ttree$ttree,rOff = off.r,pOff = off.p,shGen = w.shape,shSam = w.shape,scGen = w.scale,scSam = w.scale,pi = 1,dateT = Inf)
 expect_equal(ttree$prob,p,tolerance=0.01)
})
