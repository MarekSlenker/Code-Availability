#!/bin/bash


export ASTRALNAMEMAPFILE # namemapfile maps samples and their species


# collect all BSlistfiles
cat astral.Assembly_*/BSlistfiles > BSlistfiles

# collect all besttrees
cat raxml.Assembly_*/besttrees.tre > besttrees.tre


module add jdk-8

java -jar ./PhyloNetworks/Astral_binary/astral.5.7.8.jar \
-i besttrees.tre -b BSlistfiles \
-r 100 \
-o astral.tre \
-a "$ASTRALNAMEMAPFILE"


exit
