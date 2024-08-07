#!/bin/bash


export ALIGNMENT # fasta file with concatenated sequences
export MODEL="<ALIGNMENT>.raxml.bestModel"  # bestModel file produced by script ML_Trees.concat.2.bestTree.sh
export PBS_NUM_PPN # number of cores
export SEED=$(($RANDOM*$RANDOM))




module add raxml-ng-8



# just for illustration. Do the following in separate jobs !!!!
for((TREE=0;TREE<500;TREE++)); do

    raxml-ng --bootstrap --msa "$ALIGNMENT" --model "$MODEL" --bs-trees 1 --blopt nr_safe --seed "$SEED" --threads $PBS_NUM_PPN --prefix "${ALIGNMENT%.*}"_"$TREE" >> "${ALIGNMENT%.*}"_"$TREE".log

done 
