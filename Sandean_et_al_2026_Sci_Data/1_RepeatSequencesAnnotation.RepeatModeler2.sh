#!/bin/bash




SEQ="C087_203_hap1.fa"


singularity run \
-B <PATH>/dfam-tetools/Libraries:/opt/RepeatMasker/Libraries \
<PATH>/dfam-tetools/dfam-tetools-latest.sif \
BuildDatabase -name "${SEQ%.*}" $SEQ      2>&1 | tee 00_BuildDatabase.log



singularity run \
-B <PATH>/dfam-tetools/Libraries:/opt/RepeatMasker/Libraries \
<PATH>/dfam-tetools/dfam-tetools-latest.sif \
RepeatModeler -database "${SEQ%.*}" -threads 16 -LTRStruct             2>&1 | tee 00_repeatmodeler.log



# C087_203_hap1-families.fa  is the output





