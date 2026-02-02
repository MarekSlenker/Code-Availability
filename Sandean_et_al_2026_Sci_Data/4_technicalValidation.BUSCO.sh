#!/bin/bash



singularity run \
<PATH>/bin/busco/busco.img busco \
-i ../C087_203_hap1.PASA.longest_isoform.aa -l brassicales_odb10 -o C087_203_hap1.final.brassicales_odb10 -m proteins --cpu 8

