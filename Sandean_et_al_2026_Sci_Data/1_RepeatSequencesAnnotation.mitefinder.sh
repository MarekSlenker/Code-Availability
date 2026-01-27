#!/bin/bash



SEQ="C087_203_hap1.fa"



<PATH>/miteFinder/bin/miteFinder_linux_x64 \
-input "$SEQ" \
-threshold 0.5 \
-output "$SEQ"_mitefinder.fasta \
-pattern_scoring <PATH>/miteFinder/profile



# C087_203_hap1_mitefinder.fasta  is the output









