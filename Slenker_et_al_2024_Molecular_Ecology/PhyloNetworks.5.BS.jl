#!/usr/bin/env julia


finalTree = ________  # final tree with the optimal number of reticulation events
hmax = ____           # maximum number of hybridizations allowed
rruns = _____ # number of independent starting points for the search
run = ____   # integer, number of run. We used to submit the BS calculations in 500 independent jobs. 

using PhyloNetworks
using CSV
using DataFrames
using PhyloPlots
using RCall
using Random




net = readTopology(finalTree)

bootTrees = readBootstrapTrees("astral.BSlistfiles")

ssed = rand(1234:9999)

name = string("net", hhmax, ".bootsnaq.run", run, ".sed", ssed, ".snaq")

bootnet = bootsnaq(net, bootTrees, hmax=hmax, nrep=1, runs=rruns,
                   filename=name, seed=ssed)


exit()
