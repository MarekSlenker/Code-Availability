#!/bin/bash




mv OrthoDB12.Viridiplantae.fa proteins.fasta
cat uniprot_SwissProt.fasta >> proteins.fasta

# sequences from seven related species, as mentioned in the Sandean et al., 2026. Sci. Data
cat *protein.faa >> proteins.fasta



wd=Braker3
singularity run \
<PATH>/bin/BRAKER3/braker3.sif \
braker.pl \
--AUGUSTUS_CONFIG_PATH=<PATH>/bin/BRAKER3/Augustus/config \
--genome=C087_203_hap1.fa.masked.fa \
--bam=Cardamine_RNASeq.sorted.bam \
--prot_seq=proteins.fasta \
--workingdir=${wd} --GENEMARK_PATH=/opt/ETP/bin --threads $PBS_NUM_PPN --busco_lineage brassicales_odb10 \
--gff3 --verbosity=4 --nocleanup

