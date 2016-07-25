# TransPhylo

# Introduction

This is the homepage of TransPhylo, a software package that can reconstruct infectious disease transmission using genomic data. The input is a dated phylogeny, where leaves correspond to pathogens sampled from the known infected hosts. The main output is a transmission tree which indicates who infected whom, including the potential existence of unsampled individuals who may have acted as missing transmission links. TransPhylo works by colouring the branches of the phylogeny using a separate colour for each host, sampled or not. Each section of the tree  coloured in a unique colour represents the pathogen evolution happening within a distinct host. Changes of colours on branches therefore correspond to transmission events from one host to another.

For a more formal description of TransPhylo, see the following preprint:

Didelot, Fraser, Gardy and Colijn (2016)
Genomic infectious disease epidemiology in partially sampled and ongoing outbreaks
http://biorxiv.org/content/early/2016/07/22/065334

# Installation

You can install TransPhylo in R using the following command:

`devtools::install_github('xavierdidelot/TransPhylo')`

# Tutorial

See R vignette for now. More coming soon.

# Getting help

If you need assistance using TransPhylo, you can get in touch by emailing `xavier.didelot@gmail.com`


# TransPhyloMatlab

If you are looking for the older Matlab version of TransPhylo described in MBE 2014 31:1869-1879, please note that this has now been moved to the repository TransPhyloMatlab, available at https://github.com/xavierdidelot/TransPhyloMatlab
