#!/bin/bash



module add jdk/8
module add mambaforge
mamba activate gemoma


java -jar <PATH>/bin/gemoma2.3/GeMoMa-1.9.jar CLI GeMoMa \
    s=${PROT}/search.txt \
    t=$GENOM \
    c=${PROT}/cds-parts.fasta \
    a=${PROT}/assignment.tabular \
    sort=false \
    Score="Trust" \
    outdir=${PROT} \
    i=introns.gff \
    coverage=UNSTRANDED \
    coverage_unstranded=./coverage.bedgraph

