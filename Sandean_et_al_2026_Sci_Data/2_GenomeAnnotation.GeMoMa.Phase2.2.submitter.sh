#!/bin/bash




for PROT in $(find ./ -type d | tail -n+2 | sort); do
    echo $PROT
    export PROT

    qsub -l walltime=440:0:0 -l select=1:ncpus=1:mem=24gb:scratch_local=100gb -V 2_GenomeAnnotation.GeMoMa.Phase2.2.sh

done