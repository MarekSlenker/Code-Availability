#!/bin/bash

# Input trees are bestTrees made by RAxML. We used consensus sequences of 2x samples, and phased sequences for polyploids. 
# however, dataset contains missing samples, so we have to check presence of each sample in tree


for f in *bestTree; do
touch namemapfile
toPrint=""
if grep -q carniolicum "$f"; then toPrint=$toPrint"carniolicum2x_carniolicum:"; fi
if grep -q  Podjelak2 "$f"; then toPrint=$toPrint"Podjelak2_carniolicum,"; fi
if grep -q  Svibno7 "$f"; then toPrint=$toPrint"Svibno7_carniolicum,"; fi
echo $toPrint >> namemapfile
toPrint=""
if grep -q crassistylum "$f"; then toPrint=$toPrint"crassistylum2x_crassistylum:"; fi
if grep -q  652RDU2 "$f"; then toPrint=$toPrint"652RDU2_crassistylum,"; fi
echo $toPrint >> namemapfile
toPrint=""
if grep -q cuspidatum "$f"; then toPrint=$toPrint"cuspidatum2x_cuspidatum:"; fi
if grep -q  KONJ14 "$f"; then toPrint=$toPrint"KONJ14_cuspidatum,"; fi
if grep -q  SVI9s "$f"; then toPrint=$toPrint"SVI9s_cuspidatum,"; fi
echo $toPrint >> namemapfile
toPrint=""

if grep -q kuemmerlei "$f"; then toPrint=$toPrint"kuemmerlei2x_kuemmerlei:"; fi
if grep -q  PRST2_kuemmerlei "$f"; then toPrint=$toPrint"PRST2_kuemmerlei,"; fi
if grep -q  SHTRs3_kuemmerlei "$f"; then toPrint=$toPrint"SHTRs3_kuemmerlei,"; fi
echo $toPrint >> namemapfile
toPrint=""

if grep -q linariifolium "$f"; then toPrint=$toPrint"linariifolium2x_linariifolium:"; fi
if grep -q  DUCA8_linariifolium "$f"; then toPrint=$toPrint"DUCA8_linariifolium,"; fi
if grep -q  Kavlar1lin_linariifolium "$f"; then toPrint=$toPrint"Kavlar1lin_linariifolium,"; fi
echo $toPrint >> namemapfile
toPrint=""

if grep -q malcolmia "$f"; then toPrint=$toPrint"malcolmia2x_malcolmia:"; fi
if grep -q  MalcolmiaI1_malcolmia "$f"; then toPrint=$toPrint"MalcolmiaI1_malcolmia,"; fi
if grep -q  MalcolmiaII1_malcolmia "$f"; then toPrint=$toPrint"MalcolmiaII1_malcolmia,"; fi
echo $toPrint >> namemapfile
toPrint=""

if grep -q microstylum "$f"; then toPrint=$toPrint"microstylum2x_microstylum:"; fi
if grep -q  Alona1_microstylum "$f"; then toPrint=$toPrint"Alona1_microstylum,"; fi
echo $toPrint >> namemapfile
toPrint=""

if grep -q pectinatum "$f"; then toPrint=$toPrint"pectinatum2x_pectinatum:"; fi
if grep -q  PTER1_pectinatum "$f"; then toPrint=$toPrint"PTER1_pectinatum,"; fi
if grep -q  TGTS1_pectinatum "$f"; then toPrint=$toPrint"TGTS1_pectinatum,"; fi
echo $toPrint >> namemapfile
toPrint=""

if grep -q saxosum "$f"; then toPrint=$toPrint"saxosum2x_saxosum:"; fi
if grep -q  BUILA1s_saxosum "$f"; then toPrint=$toPrint"BUILA1s_saxosum,"; fi
if grep -q  CAKR14_saxosum "$f"; then toPrint=$toPrint"CAKR14_saxosum,"; fi
echo $toPrint >> namemapfile
toPrint=""

if grep -q sylvestre "$f"; then toPrint=$toPrint"sylvestre2x_sylvestre:"; fi
if grep -q  Dovlici1_sylvestre "$f"; then toPrint=$toPrint"Dovlici1_sylvestre,"; fi
if grep -q  Sabotin11_sylvestre "$f"; then toPrint=$toPrint"Sabotin11_sylvestre,"; fi
echo $toPrint >> namemapfile
toPrint=""

if grep -q vitekii "$f"; then toPrint=$toPrint"vitekii2x_vitekii:"; fi
if grep -q  CRKV3_vitekii "$f"; then toPrint=$toPrint"CRKV3_vitekii,"; fi
if grep -q  Neves1_vitekii "$f"; then toPrint=$toPrint"Neves1_vitekii,"; fi
echo $toPrint >> namemapfile
toPrint=""

if grep -q witmannii "$f"; then toPrint=$toPrint"witmannii2x_witmannii:"; fi
if grep -q  GEMHs1_witmannii "$f"; then toPrint=$toPrint"GEMHs1_witmannii,"; fi
if grep -q  LESc_witmannii "$f"; then toPrint=$toPrint"LESc_witmannii,"; fi
if grep -q  NASZ8_witmannii "$f"; then toPrint=$toPrint"NASZ8_witmannii,"; fi
if grep -q  SUGA1sm_witmannii "$f"; then toPrint=$toPrint"SUGA1sm_witmannii,"; fi
if grep -q  TUN2s_witmannii "$f"; then toPrint=$toPrint"TUN2s_witmannii,"; fi
if grep -q  ZARN6s_witmannii "$f"; then toPrint=$toPrint"ZARN6s_witmannii,"; fi
echo $toPrint >> namemapfile
toPrint=""

if grep -q CIUC6R "$f"; then
echo "CIUC6R.h1_odoratum22chrom:CIUC6R.h1_odoratum22chrom" >> namemapfile
echo "CIUC6R.h2_odoratum22chrom:CIUC6R.h2_odoratum22chrom" >> namemapfile
echo "CIUC6R.h3_odoratum22chrom:CIUC6R.h3_odoratum22chrom" >> namemapfile
echo "CIUC6R.h4_odoratum22chrom:CIUC6R.h4_odoratum22chrom" >> namemapfile
fi

if grep -q Maliscak9 "$f"; then
echo "Maliscak9.h1_odoratum22chrom:Maliscak9.h1_odoratum22chrom" >> namemapfile
echo "Maliscak9.h2_odoratum22chrom:Maliscak9.h2_odoratum22chrom" >> namemapfile
echo "Maliscak9.h3_odoratum22chrom:Maliscak9.h3_odoratum22chrom" >> namemapfile
echo "Maliscak9.h4_odoratum22chrom:Maliscak9.h4_odoratum22chrom" >> namemapfile
fi


if grep -q BREITs1 "$f"; then
echo "BREITs1.h1_odoratum4x:BREITs1.h1_odoratum4x" >> namemapfile
echo "BREITs1.h2_odoratum4x:BREITs1.h2_odoratum4x" >> namemapfile
echo "BREITs1.h3_odoratum4x:BREITs1.h3_odoratum4x" >> namemapfile
echo "BREITs1.h4_odoratum4x:BREITs1.h4_odoratum4x" >> namemapfile
fi

if grep -q Kavlar2 "$f"; then
echo "Kavlar2.h1_odoratum4x:Kavlar2.h1_odoratum4x" >> namemapfile
echo "Kavlar2.h2_odoratum4x:Kavlar2.h2_odoratum4x" >> namemapfile
echo "Kavlar2.h3_odoratum4x:Kavlar2.h3_odoratum4x" >> namemapfile
echo "Kavlar2.h4_odoratum4x:Kavlar2.h4_odoratum4x" >> namemapfile
fi

if grep -q GLAVs1 "$f"; then
echo "GLAVs1.h1_odoratum6x:GLAVs1.h1_odoratum6x" >> namemapfile
echo "GLAVs1.h2_odoratum6x:GLAVs1.h2_odoratum6x" >> namemapfile
echo "GLAVs1.h3_odoratum6x:GLAVs1.h3_odoratum6x" >> namemapfile
echo "GLAVs1.h4_odoratum6x:GLAVs1.h4_odoratum6x" >> namemapfile
echo "GLAVs1.h5_odoratum6x:GLAVs1.h5_odoratum6x" >> namemapfile
echo "GLAVs1.h6_odoratum6x:GLAVs1.h6_odoratum6x" >> namemapfile
fi

if grep -q VGrad1 "$f"; then
echo "VGrad1.h1_odoratum6x:VGrad1.h1_odoratum6x" >> namemapfile
echo "VGrad1.h2_odoratum6x:VGrad1.h2_odoratum6x" >> namemapfile
echo "VGrad1.h3_odoratum6x:VGrad1.h3_odoratum6x" >> namemapfile
echo "VGrad1.h4_odoratum6x:VGrad1.h4_odoratum6x" >> namemapfile
echo "VGrad1.h5_odoratum6x:VGrad1.h5_odoratum6x" >> namemapfile
echo "VGrad1.h6_odoratum6x:VGrad1.h6_odoratum6x" >> namemapfile
fi

if grep -q RET1s "$f"; then
echo "RET1s.h1_odoratumRetezat:RET1s.h1_odoratumRetezat" >> namemapfile
echo "RET1s.h2_odoratumRetezat:RET1s.h2_odoratumRetezat" >> namemapfile
echo "RET1s.h3_odoratumRetezat:RET1s.h3_odoratumRetezat" >> namemapfile
echo "RET1s.h4_odoratumRetezat:RET1s.h4_odoratumRetezat" >> namemapfile
fi

if grep -q SCO3 "$f"; then
echo "SCO3.h1_odoratumRetezat:SCO3.h1_odoratumRetezat" >> namemapfile
echo "SCO3.h2_odoratumRetezat:SCO3.h2_odoratumRetezat" >> namemapfile
echo "SCO3.h3_odoratumRetezat:SCO3.h3_odoratumRetezat" >> namemapfile
echo "SCO3.h4_odoratumRetezat:SCO3.h4_odoratumRetezat" >> namemapfile
fi


java -jar ~/bin/astral.5.7.8/astral.5.7.8.jar -i $f -o "$f".astralTree --namemapfile namemapfile -t 0
rm namemapfile; 
done
sed -i 's/carniolicum2x_carniolicum/carniolicum2x_carniolicum:001/' *ee
sed -i 's/crassistylum2x_crassistylum/crassistylum2x_crassistylum:001/' *ee
sed -i 's/cuspidatum2x_cuspidatum/cuspidatum2x_cuspidatum:001/' *ee
sed -i 's/kuemmerlei2x_kuemmerlei/kuemmerlei2x_kuemmerlei:001/' *ee
sed -i 's/linariifolium2x_linariifolium/linariifolium2x_linariifolium:001/' *ee
sed -i 's/malcolmia2x_malcolmia/malcolmia2x_malcolmia:001/' *ee
sed -i 's/microstylum2x_microstylum/microstylum2x_microstylum:001/' *ee
sed -i 's/pectinatum2x_pectinatum/pectinatum2x_pectinatum:001/' *ee
sed -i 's/saxosum2x_saxosum/saxosum2x_saxosum:001/' *ee
sed -i 's/sylvestre2x_sylvestre/sylvestre2x_sylvestre:001/' *ee
sed -i 's/vitekii2x_vitekii/vitekii2x_vitekii:001/' *ee
sed -i 's/witmannii2x_witmannii/witmannii2x_witmannii:001/' *ee

sed -i 's/CIUC6R.h1_odoratum22chrom/CIUC6R.h1_odoratum22chrom:001/' *astralTree
sed -i 's/CIUC6R.h2_odoratum22chrom/CIUC6R.h2_odoratum22chrom:001/' *astralTree
sed -i 's/CIUC6R.h3_odoratum22chrom/CIUC6R.h3_odoratum22chrom:001/' *astralTree
sed -i 's/CIUC6R.h4_odoratum22chrom/CIUC6R.h4_odoratum22chrom:001/' *astralTree
sed -i 's/Maliscak9.h1_odoratum22chrom/Maliscak9.h1_odoratum22chrom:001/' *astralTree
sed -i 's/Maliscak9.h2_odoratum22chrom/Maliscak9.h2_odoratum22chrom:001/' *astralTree
sed -i 's/Maliscak9.h3_odoratum22chrom/Maliscak9.h3_odoratum22chrom:001/' *astralTree
sed -i 's/Maliscak9.h4_odoratum22chrom/Maliscak9.h4_odoratum22chrom:001/' *astralTree
sed -i 's/BREITs1.h1_odoratum4x/BREITs1.h1_odoratum4x:001/' *astralTree
sed -i 's/BREITs1.h2_odoratum4x/BREITs1.h2_odoratum4x:001/' *astralTree
sed -i 's/BREITs1.h3_odoratum4x/BREITs1.h3_odoratum4x:001/' *astralTree
sed -i 's/BREITs1.h4_odoratum4x/BREITs1.h4_odoratum4x:001/' *astralTree
sed -i 's/Kavlar2.h1_odoratum4x/Kavlar2.h1_odoratum4x:001/' *astralTree
sed -i 's/Kavlar2.h2_odoratum4x/Kavlar2.h2_odoratum4x:001/' *astralTree
sed -i 's/Kavlar2.h3_odoratum4x/Kavlar2.h3_odoratum4x:001/' *astralTree
sed -i 's/Kavlar2.h4_odoratum4x/Kavlar2.h4_odoratum4x:001/' *astralTree
sed -i 's/GLAVs1.h1_odoratum6x/GLAVs1.h1_odoratum6x:001/' *astralTree
sed -i 's/GLAVs1.h2_odoratum6x/GLAVs1.h2_odoratum6x:001/' *astralTree
sed -i 's/GLAVs1.h3_odoratum6x/GLAVs1.h3_odoratum6x:001/' *astralTree
sed -i 's/GLAVs1.h4_odoratum6x/GLAVs1.h4_odoratum6x:001/' *astralTree
sed -i 's/GLAVs1.h5_odoratum6x/GLAVs1.h5_odoratum6x:001/' *astralTree
sed -i 's/GLAVs1.h6_odoratum6x/GLAVs1.h6_odoratum6x:001/' *astralTree
sed -i 's/VGrad1.h1_odoratum6x/VGrad1.h1_odoratum6x:001/' *astralTree
sed -i 's/VGrad1.h2_odoratum6x/VGrad1.h2_odoratum6x:001/' *astralTree
sed -i 's/VGrad1.h3_odoratum6x/VGrad1.h3_odoratum6x:001/' *astralTree
sed -i 's/VGrad1.h4_odoratum6x/VGrad1.h4_odoratum6x:001/' *astralTree
sed -i 's/VGrad1.h5_odoratum6x/VGrad1.h5_odoratum6x:001/' *astralTree
sed -i 's/VGrad1.h6_odoratum6x/VGrad1.h6_odoratum6x:001/' *astralTree
sed -i 's/RET1s.h1_odoratumRetezat/RET1s.h1_odoratumRetezat:001/' *astralTree
sed -i 's/RET1s.h2_odoratumRetezat/RET1s.h2_odoratumRetezat:001/' *astralTree
sed -i 's/RET1s.h3_odoratumRetezat/RET1s.h3_odoratumRetezat:001/' *astralTree
sed -i 's/RET1s.h4_odoratumRetezat/RET1s.h4_odoratumRetezat:001/' *astralTree
sed -i 's/SCO3.h1_odoratumRetezat/SCO3.h1_odoratumRetezat:001/' *astralTree
sed -i 's/SCO3.h2_odoratumRetezat/SCO3.h2_odoratumRetezat:001/' *astralTree
sed -i 's/SCO3.h3_odoratumRetezat/SCO3.h3_odoratumRetezat:001/' *astralTree
sed -i 's/SCO3.h4_odoratumRetezat/SCO3.h4_odoratumRetezat:001/' *astralTree


exit
