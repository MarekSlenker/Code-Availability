#!/bin/bash





RELIB="merged-library.centroids.fa.classified.corrected.fasta"
SEQ="C087_203_hap2.fa"



singularity run \
-B <PATH>/dfam-tetools/Libraries:/opt/RepeatMasker/Libraries \
<PATH>/dfam-tetools/dfam-tetools-latest.sif \
RepeatMasker \
-lib $RELIB \
-nolow \
-s \
-parallel 16 \
-gccalc \
-xsmall \
-gff \
-html \
-noisy \
$(basename $SEQ)


