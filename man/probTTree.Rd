\name{probTTree}
\alias{probTTree}
\title{Calculates the log-probability of a transmission tree}
\usage{
probTTree(ttree, R, pi, w.shape, w.scale)
}
\arguments{
  \item{ttree}{Transmission tree}

  \item{R}{Basic reproduction number}

  \item{pi}{probability of sampling an infected individual}

  \item{w.shape}{Shape parameter of the Gamma probability
  density function representing the generation length w}

  \item{w.scape}{Scale parameter of the Gamma probability
  density function representing the generation length w}
}
\value{
Probability of the transmission tree
}
\description{
Calculates the log-probability of a transmission tree
}

