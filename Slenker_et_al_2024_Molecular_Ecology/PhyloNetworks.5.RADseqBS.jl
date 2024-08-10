#!/usr/bin/env julia


using PhyloNetworks
using CSV
using DataFrames
using PhyloPlots
using RCall
using Random



finalTree = ________  # final tree with the optimal number of reticulation events
ccf = _______         # Concordance Factors
hhmax = ____           # maximum number of hybridizations allowed
rruns = _____ # number of independent starting points for the search
run = ____   # integer, number of run. We used to submit the BS calculations in 500 independent jobs. 



net = readTopology(finalTree)
buckyDat = CSV.read(ccf, DataFrame) # names like: CF12.34, CF12.34_lo etc.
rename!(x -> Symbol(replace(String(x), "." => "_")), buckyDat) # bootsnaq requires these colunm names

ssed = rand(1234:9999)

name = string("net", hhmax, ".bootsnaq.run", run, ".sed", ssed, ".snaq")

bootnet = bootsnaq(net, buckyDat, hmax=hhmax, nrep=1, runs=rruns,
                   filename=name, seed=ssed)



exit()
