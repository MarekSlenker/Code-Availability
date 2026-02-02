#!/bin/bash




singularity run \
<PATH>/bin/trinityrnaseq-v2.15.2/trinityrnaseq.v2.15.2.simg Trinity \
--genome_guided_bam Cardamine_RNASeq.sorted.bam \
--genome_guided_max_intron 10000 \
--CPU 16 \
--max_memory 450G --output trinity.C087_203_hap1 >>stdout.log 2>stderr.log




