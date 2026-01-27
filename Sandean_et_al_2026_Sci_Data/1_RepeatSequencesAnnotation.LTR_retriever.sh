#!/bin/bash




module add blast-plus/2.12.0-gcc-10.2.1-2phsggo
module add repeatmasker/4.1.2-p1-gcc-10.2.1-uymtbpy
module add cdhit/4.8.1

module add mambaforge
mamba activate TEsorter


DATADIR="<PATH>/LTR_FINDER"
RESDIR="<PATH>/LTR_retriever"

cd $SCRATCHDIR
cp $DATADIR/* .


<PATH>/LTR_retriever -genome C087_203_hap1.fa -inharvest C087_203_hap1.fa.finder.combine.scn -threads 16 -noanno | tee LTR_retriever.log


# C087_203_hap1.fa.LTRlib.fa is the output





