\name{simulateOutbreak}
\alias{simulateOutbreak}
\title{Simulate an outbreak}
\usage{
simulateOutbreak(R = 1, neg = 0.25, ninf = NA, pi = 0.5, w.shape = 2,
  w.scale = 1, dateLastSample = 2014)
}
\arguments{
  \item{R}{the basic reproduction number}

  \item{neg}{the within-host effective population size (Ne)
  times generation duration (g)}

  \item{ninf}{number of sampled infected individuals, or NA
  for any}

  \item{pi}{probability of sampling an infected individual}

  \item{w.shape}{Shape parameter of the Gamma probability
  density function representing the generation length w}

  \item{w.scape}{Scale parameter of the Gamma probability
  density function representing the generation length w}

  \item{dateLastSample}{Date given to the last sample}
}
\value{
Combined phylogenetic and transmission tree
}
\description{
Simulate an outbreak
}
\examples{
plotBothTree(simulateOutbreak())
}

