#!/bin/bash



# <PATH>/dfam-tetools/Libraries is the location of Dfam (v3.9) and RepBase (v20181026) databases. Databases were configured with /opt/RepeatMasker/tetoolsDfamUpdate.pl

singularity run \
-B <PATH>/dfam-tetools/Libraries:/opt/RepeatMasker/Libraries \
<PATH>/dfam-tetools/dfam-tetools-latest.sif

./famdb.py \
-i Libraries/famdb families \
--format fasta_name \
--ancestors --descendants Brassicaceae --include-class-in-name >~/Brassicaceae.famdbLib.fasta

