#!/usr/bin/env julia


using PhyloNetworks
using CSV
using DataFrames
using PhyloPlots
using RCall



besttrees = readMultiTopology("besttrees.tre");
q,t = countquartetsintrees(besttrees);
individualsCF = writeTableCF(q,t)
CSV.write("tableCF.astral.individuals.csv", individualsCF); # to save the data frame to a file

mapAllelesCFtable("allele-species-map.csv", "tableCF.astral.individuals.csv"; filename = "tableCF.astral.speciesNames.csv")


exit()
