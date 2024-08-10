#!/usr/bin/env julia



using PhyloNetworks
using CSV
using DataFrames
using PhyloPlots
using RCall



bootnet = readMultiTopology("bootTrees");
net = readTopology("net3.RAxML.FINAL.snaq.out")


# default figure
R"pdf"("plot1.pdf", width=15, height=10);R"par(mar=c(0,0,0,1))";plot(net, :R, showNodeNumber=true);R"dev.off()"; 


# root and rotate
rootatnode!(net, 38)

rotate!(net, 47)
rotate!(net, 43)
rotate!(net, 42)
rotate!(net, 51)
rotate!(net, 52)

R"pdf"("plot1b.pdf", width=15, height=10);R"par(mar=c(0,0,0,1))";plot(net, :R, showNodeNumber=true, showGamma=true);R"dev.off()"; 



BSe_tree, tree1 = treeEdgesBootstrap(bootnet,net)
show(BSe_tree, allrows=true, allcols=true)

BSe_tree[BSe_tree[!,"proportion"] .< 100, :]


R"pdf"("plot2.pdf", width=9, height=6);R"par(mar=c(0,0,0,1))";plot(net, :R, edgeLabel=BSe_tree);R"dev.off()"; 
R"pdf"("plot2b.pdf", width=9, height=6);R"par(mar=c(0,0,0,1))";plot(net, :R; showNodeNumber = true, showEdgeNumber = true);R"dev.off()"; 


BSn, BSe, BSc, BSgam, BSedgenum = hybridBootstrapSupport(bootnet, net);



R"pdf"("plot3.pdf", width=9, height=6);R"par(mar=c(0,0,0,1))";plot(net, :R, edgeLabel=BSe[!,["edge", "BS_hybrid_edge"]]);R"dev.off()"; 


BSe[!,["edge", "BS_hybrid_edge"]]
CSV.write("BSe.txt", BSe)


tmp = filter(row -> !ismissing(row[:edge]), BSe) # filter rows
select!(tmp, [:edge,:BS_hybrid_edge])            # select 2 columns only
rename!(tmp, :BS_hybrid_edge => :proportion)     # rename those columns, to match names in BSe_tree
rename!(tmp, :edge => :edgeNumber)
tmp = vcat(BSe_tree, tmp)


R"pdf"("plot_final.pdf", width=9, height=6);R"par(mar=c(0,0,0,1))";plot(net, :R, edgeLabel=tmp, nodeLabel=BSn[:, [:hybridnode,:BS_hybrid_samesisters]], showGamma=true);R"dev.off()"; 


writeTopology(net, "net3.RAxML.FINAL.snaq.BS.out")


exit()

