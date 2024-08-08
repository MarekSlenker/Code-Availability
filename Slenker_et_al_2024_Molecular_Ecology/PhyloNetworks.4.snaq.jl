#!/usr/bin/env julia


hmax = ____  # maximum number of hybridizations allowed
run = ____   # integer, number of run. We used to submit the calculations in 50 independent jobs. 
runs = _____ # number of independent starting points for the search

using PhyloNetworks
using CSV
using DataFrames
using PhyloPlots
using RCall
using Random


# !!!!!! use only one of the two following rows
HybridNetwork = readMultiTopology(nnet)[102]  # if the input is H0 astral tree
HybridNetwork = readTopology(nnet)            # if the input is network with hmax-1 reticulation events


raxmlCF = readTableCF("tableCF.astral.speciesNames.csv")


name = string("net", hmax, ".RAxML.", run, ".", rand(1234:9999), ".snaq")

netRes = snaq!(HybridNetwork,  raxmlCF, hmax=hmax, filename=name, runs=runs)

exit()

