context("Test probTTree")

test_that("Probability of a transmission tree is the same when simulating and calculating.", {
  #Note this test does not work for pi>0 or maxTime<Inf because then there are unsampled cases not accounted for in the calculation of probability within the simulation procedure
  set.seed(0)
  pi=1
  off.r=1
  off.p=0.5
  w.shape=1.1
  w.scale=1.2
  ttree=makeTTree(pi=pi,off.r=off.r,off.p=off.p,w.shape=w.shape,w.scale=w.scale,ws.shape=w.shape,ws.scale=w.scale,maxTime=Inf,nSampled = NA)
  p=probTTree (ttree$ttree,off.r,off.p,pi,w.shape,w.scale,w.shape,w.scale,Inf)
  p2=TransPhylo:::probTTreeR(ttree$ttree,off.r,off.p,pi,w.shape,w.scale,w.shape,w.scale,Inf)
  expect_equal(ttree$prob,p,tolerance=0.001,scale=1)
  expect_equal(ttree$prob,p2,tolerance=0.001,scale=1)
  
  set.seed(6)
  off.r=1.2
  off.p=0.5
  ttree=makeTTree(pi=pi,off.r=off.r,off.p=off.p,w.shape=w.shape,w.scale=w.scale,ws.shape=w.shape,ws.scale=w.scale,maxTime=Inf,nSampled = NA)
  p=probTTree (ttree$ttree,off.r,off.p,pi,w.shape,w.scale,w.shape,w.scale,Inf)
  p2=TransPhylo:::probTTreeR(ttree$ttree,off.r,off.p,pi,w.shape,w.scale,w.shape,w.scale,Inf)
  expect_equal(ttree$prob,p,tolerance=0.001,scale=1)
  expect_equal(ttree$prob,p2,tolerance=0.001,scale=1)
  
  set.seed(0)
  off.r=0.3333
  off.p=0.25
  ttree=makeTTree(pi=pi,off.r=off.r,off.p=off.p,w.shape=w.shape,w.scale=w.scale,ws.shape=w.shape,ws.scale=w.scale,maxTime=Inf,nSampled = NA)
  p=probTTree (ttree$ttree,off.r,off.p,pi,w.shape,w.scale,w.shape,w.scale,Inf)
  p2=TransPhylo:::probTTreeR(ttree$ttree,off.r,off.p,pi,w.shape,w.scale,w.shape,w.scale,Inf)
  expect_equal(ttree$prob,p,tolerance=0.001,scale=1)
  expect_equal(ttree$prob,p2,tolerance=0.001,scale=1)
})

test_that("Probability of a transmission tree is roughly the same if timeT=Inf and timeT=10000.", {
  set.seed(1)
  pi=1
  off.r=1
  off.p=0.5
  w.shape=1.1
  w.scale=1.2
  ttree=makeTTree(pi=pi,off.r=off.r,off.p=off.p,w.shape=w.shape,w.scale=w.scale,ws.shape=w.shape,ws.scale=w.scale,maxTime=Inf,nSampled = NA)
  p =TransPhylo:::probTTreeR(ttree$ttree,1.1,0.4,0.8,w.shape,w.scale,w.shape,w.scale,Inf)
  p2=TransPhylo:::probTTreeR(ttree$ttree,1.1,0.4,0.8,w.shape,w.scale,w.shape,w.scale,1000)
  expect_equal(p,p2,tolerance=0.001,scale=1)
  p =probTTree(ttree$ttree,1.1,0.4,0.8,w.shape,w.scale,w.shape,w.scale,Inf)
  p2=probTTree(ttree$ttree,1.1,0.4,0.8,w.shape,w.scale,w.shape,w.scale,100)
  expect_equal(p,p2,tolerance=0.001,scale=1)
})


test_that("Consistency with exact solution when r=1.", {
  set.seed(0)
  pi=0.4
  off.r=2
  off.p=0.4
  w.shape=1.1
  w.scale=1.2
  ttree=makeTTree(pi=pi,off.r=off.r,off.p=off.p,w.shape=w.shape,w.scale=w.scale,ws.shape=w.shape,ws.scale=w.scale,maxTime=1,nSampled = NA)
  p =probTTree (ttree$ttree,1,off.p,pi,w.shape,w.scale,w.shape,w.scale,1)
  p2=probTTree (ttree$ttree,1.0001,off.p,pi,w.shape,w.scale,w.shape,w.scale,1)
  expect_equal(p,p2,tolerance=0.001)
})

test_that("Probability of a transmission tree is the same in R and C++.", {
  set.seed(0)
  pi=0.4
  off.r=2
  off.p=0.4
  w.shape=1.1
  w.scale=1.2
  ttree=makeTTree(pi=pi,off.r=off.r,off.p=off.p,w.shape=w.shape,w.scale=w.scale,ws.shape=w.shape,ws.scale=w.scale,maxTime=5,nSampled = NA)
  p =probTTree (ttree$ttree,off.r,off.p,pi,w.shape,w.scale,w.shape,w.scale,Inf)
  p2=TransPhylo:::probTTreeR(ttree$ttree,off.r,off.p,pi,w.shape,w.scale,w.shape,w.scale,Inf)
  expect_equal(p,p2,tolerance=0.001)
  set.seed(0)
  ttree=makeTTree(pi=pi,off.r=off.r,off.p=off.p,w.shape=w.shape,w.scale=w.scale,ws.shape=w.shape,ws.scale=w.scale,maxTime=3,nSampled = NA)
  p =probTTree (ttree$ttree,off.r,off.p,pi,w.shape,w.scale,w.shape,w.scale,3)
  p2=TransPhylo:::probTTreeR(ttree$ttree,off.r,off.p,pi,w.shape,w.scale,w.shape,w.scale,3)
  expect_equal(p,p2,tolerance=0.01)
})
