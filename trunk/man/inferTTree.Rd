\name{inferTTree}
\alias{inferTTree}
\title{Infer transmission tree given a phylogenetic tree}
\usage{
inferTTree(ptree, dateLastSample = 2014, mcmcIterations = 1000,
  startNeg = 100/365, startR = 1, startPi = 0.5, updateNeg = TRUE,
  updateR = TRUE, updatePi = TRUE, testing = FALSE)
}
\arguments{
  \item{ptree}{Phylogenetic tree}

  \item{dateLastSample}{Date of the last sample}

  \item{mcmcIterations}{Number of MCMC iterations to run
  the algorithm for}

  \item{startNeg}{Starting value of within-host coalescent
  parameter Ne*g}

  \item{startR}{Starting value of basic reproduction number
  R}

  \item{startPi}{Starting value of sampling proportion pi}

  \item{updateNeg}{Whether of not to update the parameter
  Ne*g}

  \item{updateR}{Whether or not to update the parameter R}

  \item{updatePi}{Whether or not to update the parameter
  pi}

  \item{testing}{Whether or not to apply the test mode
  without likelihood}
}
\value{
posterior sample set of transmission trees
}
\description{
Infer transmission tree given a phylogenetic tree
}

