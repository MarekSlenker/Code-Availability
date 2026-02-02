

module add python/3.7.7-intel-19.0.4-mgiwa7z



# create_dbs.py -m diamond --dbname bacteria --taxa 


emapper.py --cpu 24 \
--data_dir <PATH>/bin/eggnog-mapper/data \
-o C087_203_hap1 \
--output_dir <PATH>/function_annotation/eggNOG \
--temp_dir <PATH>/function_annotation/eggNOG \
--override -m diamond --dmnd_ignore_warnings \
-i ../../C087_203_hap1.post_PASA.longest_isoform.aa \
--evalue 0.00001 --score 60 --pident 40 --query_cover 20 --subject_cover 20 --itype proteins \
--tax_scope auto --target_orthologs all --go_evidence non-electronic --pfam_realign none \
--report_orthologs \
--decorate_gff yes --excel \
> <PATH>/function_annotation/eggNOG/emapper.out \
2> <PATH>/function_annotation/eggNOG/emapper.err 






