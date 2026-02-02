

<PATH>/bin/diamond/diamond makedb \
--in uniprot_trembl.fasta \
--threads 80 \
--tmpdir <PATH>/AnnotationDTBSs/RefSeq/tmp \
--db trembl


cd <PATH>/<PATH>/function_annotation/refseq
<PATH>/bin/diamond/diamond blastp \
-d <PATH>/AnnotationDTBSs/RefSeq/refseq \
-q ../../C087_203_hap1.PASA.longest_isoform.aa \
--evalue 0.00001 \
--outfmt 6 qseqid qlen qstart qend sseqid slen sstart send pident ppos length score evalue stitle \
--salltitles --max-target-seqs 1 \
-o refseq.matches.tsv


cd <PATH>/<PATH>/function_annotation/swissprot
<PATH>/bin/diamond/diamond blastp \
-d <PATH>/AnnotationDTBSs/Swiss-Prot/swissprot \
-q ../../C087_203_hap1.PASA.longest_isoform.aa \
--evalue 0.00001 \
--outfmt 6 qseqid qlen qstart qend sseqid slen sstart send pident ppos length score evalue stitle \
--max-target-seqs 1 \
-o swissprot.matches.tsv



cd <PATH>/<PATH>/function_annotation/TrEMBL
<PATH>/bin/diamond/diamond blastp \
-d <PATH>/AnnotationDTBSs/TrEMBL/trembl \
-q ../../C087_203_hap1.PASA.longest_isoform.aa \
--evalue 0.00001 \
--outfmt 6 qseqid qlen qstart qend sseqid slen sstart send pident ppos length score evalue stitle \
--max-target-seqs 1 \
-o trembl.matches.tsv






