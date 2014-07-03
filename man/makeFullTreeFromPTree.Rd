\name{makeFullTreeFromPTree}
\alias{makeFullTreeFromPTree}
\title{Deterministically create a transmission tree compatible with the provided phylogenetic tree}
\usage{
makeFullTreeFromPTree(tree, dateLastSample = 2014)
}
\arguments{
  \item{tree}{phylogenetic tree}

  \item{dateLastSample}{Date of the last sample}
}
\value{
A non-zero probability phylogenetic+transmission tree
}
\description{
Deterministically create a transmission tree compatible
with the provided phylogenetic tree
}

