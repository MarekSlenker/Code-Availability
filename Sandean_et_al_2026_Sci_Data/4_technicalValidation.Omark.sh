#!/bin/bash



module add mambaforge
mamba activate omark





<PATH>/bin/omamer-2.1.0/bin/omamer search \
--db  <PATH>/bin/omamer-2.1.0/LUCA.h5 \
--query ../C087_203_hap1.PASA.aa \
--out C087_203_hap1.PASA.aa.omamer


omark \
-f C087_203_hap1.PASA.aa.omamer \
-d <PATH>/bin/omamer-2.1.0/LUCA.h5 \
-i ../isoforms.txt \
-o .


