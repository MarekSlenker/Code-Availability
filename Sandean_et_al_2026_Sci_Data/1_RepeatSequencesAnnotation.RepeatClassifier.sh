#!/bin/bash




RELIB="merged-library.fa"



# Classification using RepeatClassifier

singularity run \
-B /auto/pruhonice1-ibot/nfs4/home/mslenker/Projects/Anotacie2/bin/dfam-tetools/Libraries:/opt/RepeatMasker/Libraries \
/auto/pruhonice1-ibot/nfs4/home/mslenker/Projects/Anotacie2/bin/dfam-tetools/dfam-tetools-latest.sif \
RepeatClassifier \
-consensi $RELIB





