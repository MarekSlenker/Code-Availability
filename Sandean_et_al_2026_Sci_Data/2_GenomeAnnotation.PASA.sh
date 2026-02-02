#!/bin/bash




cp <PATH>/C087_203_hap1.EVM.gff3 .   # EVM output, first round. PASA will do 4 rounds
cp <PATH>/my_pasaDB.CHROM.sqlite.gene_structures_post_PASA_updates.*.gff3 # output of previous PASA run

cp <PATH>/C087_203_hap1.fa.masked.fa .
cp <PATH>/trinity/trinity.hap1-GG.fasta .

cp <PATH>/alignAssembly.config .
cp <PATH>/annotCompare.config .

mkdir DB  
sed -i "s@PATH@"$PWD"@"  alignAssembly.config
sed -i "s@PATH@"$PWD"@"  annotCompare.config

# 0) generate the transcript sequences --- by trinity

# executed in a separate script:    3_functionAnnotation.trinity.sh


# 1) cleaning the transcript sequences --- by trinity

singularity exec -B $PWD \
<PATH>/bin/PASA/pasapipeline_latest.sif /usr/local/src/PASApipeline/bin/seqclean trinity.hap1-GG.fasta -c 8 



# 2) Transcript alignments followed by alignment assembly

export  CPU=8
singularity exec -B $PWD \
<PATH>/bin/PASA/pasapipeline_latest.sif /usr/local/src/PASApipeline/Launch_PASA_pipeline.pl \
-c alignAssembly.config -R -C -g C087_203_hap1.fa.masked.fa -t trinity.hap1-GG.fasta.clean -T -u trinity.hap1-GG.fasta --ALIGNERS gmap,blat --CPU 8 >>Launch_PASA_pipeline.stdout 2>Launch_PASA_pipeline.stderr



# !!!!!!
# !!!!!!
# !!!!!!
# !!!!!!
# It usually requires at least two cycles of annotation loading, annotation comparison, and annotation updates
# in order to maximize the incorporation of transcript alignments into gene structures.  Updates made to gene structures 
# in the first round often lead to the capacity to incorporate additional transcript alignments that did not fit well in the context of the earlier gene structures.
# You can use the PASA-updated annotations in the GFF3 file created at the end of the annotation comparison step as input for a subsequent annotation comparison round.
# !!!!!!
# !!!!!!
# !!!!!!
# !!!!!!



# 3) Annotation Comparisons and Annotation Updates


# Loading preexisting protein-coding gene annotations

# first round, read EVM's gff3 file  (-P)
singularity exec -B $PWD \
<PATH>/bin/PASA/pasapipeline_latest.sif /usr/local/src/PASApipeline/scripts/Load_Current_Gene_Annotations.dbi \
-c alignAssembly.config -g C087_203_hap1.fa.masked.fa \
-P hap*.gff3  >>Load_Current_Gene_Annotations.stdout 2>Load_Current_Gene_Annotations.stderr


# second and following rounds, read results from the previous round (-P)
singularity exec -B $PWD \
<PATH>/bin/PASA/pasapipeline_latest.sif /usr/local/src/PASApipeline/scripts/Load_Current_Gene_Annotations.dbi \
-c alignAssembly.config -g C087_203_hap1.fa.masked.fa \
-P my_pasaDB.CHROM.sqlite.gene_structures_post_PASA_updates.*.gff3  >>Load_Current_Gene_Annotations.stdout 2>Load_Current_Gene_Annotations.stderr


# Performing an annotation comparison and generating an updated gene set


export  CPU=1
singularity exec -B $PWD \
<PATH>/bin/PASA/pasapipeline_latest.sif /usr/local/src/PASApipeline/Launch_PASA_pipeline.pl \
        -c annotCompare.config \
        -A \
        -g C087_203_hap1.fa.masked.fa \
        -t trinity.hap1-GG.fasta.clean \
        --CPU 8  >>Launch_PASA_pipeline.2.stdout 2>Launch_PASA_pipeline.2.stderr




# PASA will output a new GFF3 file that contains the PASA-updated version of the genome annotation, 
# including those gene models successfully updated by PASA, and those that remained untouched.  
# This file will be named '${mysql_db}.gene_structures_post_PASA_updates.$pid.gff3', 
# where $pid is the process ID for this annotation comparison computation.





exit













# 0) Initialize the PASA SQLite Schema

singularity exec -B $PWD \
<PATH>/bin/PASA/pasapipeline_latest.sif /usr/local/src/PASApipeline/Launch_PASA_pipeline.pl \
-c alignAssembly.config -R -C -g C087_203_hap1.fa.masked.fa -t trinity.hap1-GG.fasta.clean --ALIGNERS gmap,blat



# 1)  be sure to check them for PASA compatibility like so:
singularity exec -B $PWD \
<PATH>/bin/PASA/pasapipeline_latest.sif /usr/local/src/PASApipeline/misc_utilities/pasa_gff3_validator.pl \
C087_203_hap1.fa.masked.fa # orig_annotations


singularity exec -B $PWD \
<PATH>/bin/PASA/pasapipeline_latest.sif /usr/local/src/PASApipeline/scripts/Load_Current_Gene_Annotations.dbi \
     -c alignAssembly.config -g C087_203_hap1.fa.masked.fa \
     -P C087_203_hap1.fa.masked.fa        # orig_annotations



# 2) Now that the original annotations are loaded, we can perform a comparison of the PASA alignment assemblies to these preexisting gene annotations


# 2.1 cleaning the transcript sequences --- by trinity

singularity exec -B $PWD \
<PATH>/bin/PASA/pasapipeline_latest.sif /usr/local/src/PASApipeline/bin/seqclean  trinity.hap1-GG.fasta




# 2.2

singularity exec -B $PWD \
<PATH>/bin/PASA/pasapipeline_latest.sif /usr/local/src/PASApipeline/Launch_PASA_pipeline.pl \
        -c annotCompare.config -A \
        -g C087_203_hap1.fa.masked.fa \
        -t trinity.hap1-GG.fasta.clean






