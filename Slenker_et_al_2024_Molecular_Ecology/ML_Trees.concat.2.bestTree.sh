#!/bin/bash


export ALIGNMENT # fasta file with concatenated sequences
export MODEL="partitions.best_scheme"
export PBS_NUM_PPN # number of cores


# Going to working directory $SCRATCHDIR
cd "$SCRATCHDIR"

cp $ALIGNMENT .
cp $MODEL .


module add raxml-ng-8

raxml-ng --msa "$ALIGNMENT" --model "$MODEL" --tree pars{10},rand{10} --blopt nr_safe --threads $PBS_NUM_PPN --redo --prefix ${ALIGNMENT%.*} >> "${ALIGNMENT%.*}".log



exit

