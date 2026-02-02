#!/bin/bash


module add parallel-20200322
module add samtools-1.9

cp <PATH>/Cardamine-amara-transcriptom/*trm* .
cp <PATH>/bioproject_PRJDB4989_Zurich/*trm* .
cp <PATH>/bioproject_PRJDB9426_underboden/*trm* .
cp <PATH>/bioproject_PRJNA575831_mandakova/*trm* .



parallel -j $PBS_NUM_PPN bunzip2 ::: *.bz2


R1=""
for f in *trm.R1*; do R1=$R1,$f; done
RR1=$(echo $R1 |sed 's/^,//' )


R2=""
for f in *trm.R2*; do R2=$R2,$f; done
RR2=$(echo $R2 |sed 's/^,//' )



singularity run \
<PATH>/bin/BRAKER3/braker3.sif /opt/ETP/tools/hisat2 \
-x C087_203_hap1.fa.masked \
-1 $RR1 \
-2 $RR2 \
--threads $PBS_NUM_PPN --dta \
-S Cardamine_RNASeq.sam


samtools view -bS Cardamine_RNASeq.sam > Cardamine_RNASeq.bam
samtools sort Cardamine_RNASeq.bam -o Cardamine_RNASeq.sorted.bam
samtools index Cardamine_RNASeq.sorted.bam


