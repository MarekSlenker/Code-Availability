#!/bin/bash




SEQ="C087_203_hap1.fa"

# <PATH>/dfam-tetools/Libraries is the location of Dfam (v3.9) and RepBase (v20181026) databases. Databases were configured with /opt/RepeatMasker/tetoolsDfamUpdate.pl

singularity run \
-B <PATH>/dfam-tetools/Libraries:/opt/RepeatMasker/Libraries \
<PATH>/dfam-tetools/dfam-tetools-latest.sif \
BuildDatabase -name "${SEQ%.*}" $SEQ      2>&1 | tee 00_BuildDatabase.log



singularity run \
-B <PATH>/dfam-tetools/Libraries:/opt/RepeatMasker/Libraries \
<PATH>/dfam-tetools/dfam-tetools-latest.sif \
RepeatModeler -database "${SEQ%.*}" -threads 16 -LTRStruct             2>&1 | tee 00_repeatmodeler.log



# C087_203_hap1-families.fa  is the output





