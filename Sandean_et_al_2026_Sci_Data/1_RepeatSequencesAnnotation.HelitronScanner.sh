#!/bin/bash



SEQ="C087_203_hap1.fa"

module add jdk/8 




## Step 1: Execute the HelitronScanner scanHead command
java -jar <PATH>/HelitronScanner_V1.0/HelitronScanner/HelitronScanner.jar \
scanHead \
-lf <PATH>/HelitronScanner_V1.0/TrainingSet/head.lcvs \
-g $SEQ \
-bs 0 \
-o "$SEQ"_helitron.head \
-tl $PBS_NUM_PPN


## Step 2: Execute the HelitronScanner scanTail command
java -jar <PATH>/HelitronScanner_V1.0/HelitronScanner/HelitronScanner.jar \
scanTail \
-lf <PATH>/HelitronScanner_V1.0/TrainingSet/tail.lcvs \
-g $SEQ \
-bs 0 \
-o "$SEQ"_helitron.tail \
-tl $PBS_NUM_PPN


## Step 3: Execute the HelitronScanner pairEnds command
java -jar <PATH>/HelitronScanner_V1.0/HelitronScanner/HelitronScanner.jar \
pairends \
-hs "$SEQ"_helitron.head \
-ts "$SEQ"_helitron.tail \
-o "$SEQ"_helitron.paired


## Step 4: Execute draw command - Create the fasta sequences from paired for each helitrons
java -jar <PATH>/HelitronScanner_V1.0/HelitronScanner/HelitronScanner.jar \
draw \
-p "$SEQ"_helitron.paired \
-g $SEQ \
-o "$SEQ"_helitron \
-pure_helitron


# C087_203_hap1.fa_helitron.hel.fa  is the output





