#!/bin/bash



module add jdk/8
module add mambaforge
mamba activate gemoma


tblastn -query ${PROT}/cds-parts.fasta -db blastdb -evalue 100.0 -out ${PROT}/search.txt -outfmt "6 std sallseqid score nident positive gaps ppos qframe sframe qseq sseq qlen slen salltitles" -db_gencode 1 -matrix BLOSUM62 -seg no -word_size 3 -comp_based_stats F -gapopen 11 -gapextend 1 -num_threads $PBS_NUM_PPN

