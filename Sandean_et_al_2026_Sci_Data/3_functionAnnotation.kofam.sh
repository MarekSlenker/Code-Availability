#!/bin/bash






module add parallel-20160622
cp <PATH>/C087_203_hap1.selected_isoforms.aa .


<PATH>/bin/kofam_scan-1.3.0/exec_annotation \
-o C087_203_hap1.selected_isoforms.aa.kofam \
--e-value 0.00001 --cpu 40 \
--format detail-tsv \
../../C087_203_hap1.selected_isoforms.aa



