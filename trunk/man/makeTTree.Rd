\name{makeTTree}
\alias{makeTTree}
\title{Simulate a transmission tree
(assumes generation distribution w is Gamma(2,1))}
\usage{
makeTTree(R, pi)
}
\arguments{
  \item{R}{the basic reproduction number}

  \item{pi}{probability of sampling an infected individual}
}
\value{
A N*3 matrix in the following format with one row per
infected host, first column is time of infection, second
column is time of sampling, third column is infector
}
\description{
Simulate a transmission tree (assumes generation
distribution w is Gamma(2,1))
}

