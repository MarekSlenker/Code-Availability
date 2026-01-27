#!/bin/bash



SEQ="C087_203_hap1.fa"



export PATH="$PATH":<PATH>/LTR_FINDER_parallel-1.3/bin/LTR_FINDER.x86_64-1.0.7


perl <PATH>/LTR_FINDER_parallel-1.3/LTR_FINDER_parallel \
-seq $SEQ -size 5000000 -time 500 -threads 8 -harvest_out


# results are forwarded to LTR_retriever





