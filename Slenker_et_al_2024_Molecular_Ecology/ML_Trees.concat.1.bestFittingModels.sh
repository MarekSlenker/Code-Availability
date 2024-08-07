#!/bin/bash


export ALIGNMENT       # fasta file with concatenated sequences
export PARTITIONS  # partition of the fasta file
export PBS_NUM_PPN # number of cores

module add iqtree/1.6.12


# Going to working directory $SCRATCHDIR
cd "$SCRATCHDIR"


cp $ALIGNMENT .
cp $PARTITIONS .


iqtree -s $ALIGNMENT -spp $PARTITIONS -m TESTMERGEONLY -st DNA -nt $PBS_NUM_PPN -pre iqtreeoutput -cptime 3600 -rclusterf 10 -safe >>iqtree.output.log 2>&1


exit
