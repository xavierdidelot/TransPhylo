context("Test probTTree")

test_that("Probability of a transmission tree is the same when simulating and calculating.", {
  #Note this test does not work for pi>0 or maxTime<Inf because then there are unsampled cases not accounted for in the calculation of probability within the simulation procedure
  set.seed(0)
  pi=1
  off.r=1
  off.p=0.5
  w.shape=1.1
  w.scale=1.1
  ttree=makeTTree(pi=pi,off.r=off.r,off.p=off.p,w.shape=w.shape,w.scale=w.scale,ws.shape=w.shape,ws.scale=w.scale,maxTime=Inf,nSampled = NA)
  p=probTTree (ttree$ttree,off.r,off.p,pi,w.shape,w.scale,w.shape,w.scale,Inf)
  expect_equal(ttree$prob,p,tolerance=0.01)
})

test_that("Probability of a transmission tree is the same in R and C++.", {
  set.seed(0)
  pi=0.4
  off.r=2
  off.p=0.4
  w.shape=1.1
  w.scale=1.1
  ttree=makeTTree(pi=pi,off.r=off.r,off.p=off.p,w.shape=w.shape,w.scale=w.scale,ws.shape=w.shape,ws.scale=w.scale,maxTime=5,nSampled = NA)
  p =probTTree (ttree$ttree,off.r,off.p,pi,w.shape,w.scale,w.shape,w.scale,Inf)
  p2=TransPhylo:::probTTreeR(ttree$ttree,off.r,off.p,pi,w.shape,w.scale,w.shape,w.scale,Inf)
  expect_equal(p,p2,tolerance=0.01)
})
