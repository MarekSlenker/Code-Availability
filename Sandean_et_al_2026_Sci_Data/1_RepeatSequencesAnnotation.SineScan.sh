#!/bin/bash



SEQ="C087_203_hap1.fa"



perl SINE_Scan-v1.1.1/SINE_Scan_process.pl \
            -g $SEQ \
            -d sinescan_output \
            -k $PBS_NUM_PPN \
            -s 123 \
            -o sinescan_output


# C087_203_hap1.sine.fa  is the output






