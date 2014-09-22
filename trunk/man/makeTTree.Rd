\name{makeTTree}
\alias{makeTTree}
\title{Simulate a transmission tree}
\usage{
makeTTree(R, pi, w.shape, w.scale)
}
\arguments{
  \item{R}{the basic reproduction number}

  \item{pi}{probability of sampling an infected individual}

  \item{w.shape}{Shape parameter of the Gamma probability
  density function representing the generation length w}

  \item{w.scape}{Scale parameter of the Gamma probability
  density function representing the generation length w}
}
\value{
A N*3 matrix in the following format with one row per
infected host, first column is time of infection, second
column is time of sampling, third column is infector
}
\description{
Simulate a transmission tree
}

