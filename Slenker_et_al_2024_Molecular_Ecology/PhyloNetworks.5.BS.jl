#!/usr/bin/env julia


finalTree = ________  # final tree with the optimal number of reticulation events
hmax = ____           # maximum number of hybridizations allowed

using PhyloNetworks
using CSV
using DataFrames
using PhyloPlots
using RCall
using Random




net = readTopology(finalTree)

bootTrees = readBootstrapTrees("astral.BSlistfiles")

ssed = rand(1234:9999)

name = string("net.1.final.bootsnaq.", ssed, ".snaq")

bootnet = bootsnaq(net, bootTrees, hmax=hmax, nrep=1, runs=15,
                   filename=name, seed=ssed)


exit()
