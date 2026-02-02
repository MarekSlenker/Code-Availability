<PATH>#!/bin/bash





module add mambaforge     # load the module
mamba activate /storage/pruhonice1-ibot/home/mslenker/.conda/envs/TD2_env




# # # # # # # # # # # # # # # # # # # # # # # # # # #
# Step 0
# # # # # # # # # # # # # # # # # # # # # # # # # # #


cp <PATH>/StringTie2/STRG.out.gtf .

# Construct the transcript fasta file using the genome and the transcripts.gtf file
<PATH>/bin/TD2/util/gtf_genome_to_cdna_fasta.pl \
STRG.out.gtf \
<PATH>/C087_203/C087_203_hap1/InputData/C087_203_hap1.fa.masked.fa \
> transcripts.fasta 

# Convert the transcript structure GTF file to an alignment-GFF3 formatted file
<PATH>/bin/TD2/util/gtf_to_alignment_gff3.pl STRG.out.gtf > STRG.out.gff3




# # # # # # # # # # # # # # # # # # # # # # # # # # #
# Step 1: extract the long open reading frames
# # # # # # # # # # # # # # # # # # # # # # # # # # #

TD2.LongOrfs -t transcripts.fasta --threads 16





# # # # # # # # # # # # # # # # # # # # # # # # # # #
# Step 2: identify ORFs with homology to known proteins via MMSeqs2, blastp or HMMER3 searches.
# # # # # # # # # # # # # # # # # # # # # # # # # # #

<PATH>/bin/mmseqs2/mmseqs/bin/mmseqs \
    databases UniProtKB/Swiss-Prot swissprot tmp

<PATH>/bin/mmseqs2/mmseqs/bin/mmseqs \
    easy-search transcripts/longest_orfs.pep swissprot alnRes.m8 tmp -s 7.0




# # # # # # # # # # # # # # # # # # # # # # # # # # #
# Step 3: predict the likely coding regions
# # # # # # # # # # # # # # # # # # # # # # # # # # #

TD2.Predict -t transcripts.fasta --retain-mmseqs-hits alnRes.m8





# # # # # # # # # # # # # # # # # # # # # # # # # # #
# Step 4: finally, generate a genome-based coding region annotation file:
# # # # # # # # # # # # # # # # # # # # # # # # # # #

<PATH>/bin/TD2/util/cdna_alignment_orf_to_genome_orf.pl \
    transcripts.fasta.TD2.gff3 \
    STRG.out.gff3 \
    transcripts.fasta \
    > transcripts.fasta.TD2.genome.gff3



