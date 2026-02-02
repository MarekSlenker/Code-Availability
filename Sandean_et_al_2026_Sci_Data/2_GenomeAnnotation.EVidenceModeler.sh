#!/bin/bash






export sample_id="hap1"
export genome="C087_203_hap1"


# braker
genePredictions1="<PATH>/Braker3/augustus.hints.gff3"
genePredictions2="<PATH>/Braker3/GeneMark-ETP/genemark.EVM.gff3"

# gemoma
export proteinAlignments="<PATH>/GeMoMa/hap1.filtered_predictions.EVM.gff3"

# Transcripts
TD2="<PATH>/TD2/transcripts.fasta.TD2.genome.gff3"



export weights="<PATH>/weights.txt"
export RESDIR="<PATH>/EVM.rawBraker3+GeMoMa+TD2"



cat $genePredictions1 $genePredictions2 >genePredictions
cat $TD2 > transcripts


singularity exec <PATH>/EVidenceModeler/EVidenceModeler.v2.1.0.simg EVidenceModeler \
--sample_id $sample_id \
--genome $(basename $genome) \
--gene_predictions genePredictions \
--protein_alignments $(basename $proteinAlignments) \
--transcript_alignments transcripts \
--weights $(basename $weights) \
--segmentSize 1000000 \
--overlapSize 10000 
--CPU $PBS_NUM_PPN






