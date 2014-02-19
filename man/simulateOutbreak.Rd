\name{simulateOutbreak}
\alias{simulateOutbreak}
\title{Simulate an outbreak}
\usage{
simulateOutbreak(R = 1, neg = 0.25, ninf = 10)
}
\arguments{
  \item{R}{the basic reproduction number}

  \item{neg}{the within-host effective population size (Ne)
  times generation duration (g)}

  \item{ninf}{number of infected individuals}
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

