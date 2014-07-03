\name{simulateOutbreak}
\alias{simulateOutbreak}
\title{Simulate an outbreak}
\usage{
simulateOutbreak(R = 1, neg = 0.25, ninf = NA, pi = 0.5,
  dateLastSample = 2014)
}
\arguments{
  \item{R}{the basic reproduction number}

  \item{neg}{the within-host effective population size (Ne)
  times generation duration (g)}

  \item{ninf}{number of sampled infected individuals, or NA
  for any}

  \item{pi}{probability of sampling an infected individual}
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

