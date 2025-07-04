#!/bin/bash



export SEQ   # fna file from 2.1.one_allopolyploid_plus_diploids


source /cvmfs/software.metacentrum.cz/modulefiles/5.1.0/loadmodules
module add perl-5.20.1-gcc
module add gcc-8.3.0
module add R-4.0.0-gcc
module add astral/5.7.7 
module add newick-utils
module add iqtree/1.6.12

# go to working directory
cd "$SCRATCHDIR"

cp "$SEQ" .

# PREVCLAF = original classification of polyploid homeolog
PREVCLAF="${SEQ%.fna}"
PREVCLAF=${PREVCLAF##*_}

# SAMPLE = label of polyploid sample
SAMPLE="${SEQ%_"$PREVCLAF".fna}"
SAMPLE=${SAMPLE##*_}

# NAME = name of sequence
NAME="${SEQ%_*.fna}"
NAME=${NAME##*>}


# we want 500 BS trees, but sometimes it fails, don't know why, so using this workaround
    let i++
    ./iqtree -s "$SEQ" -bo 500 -st DNA -AICc -nt $PBS_NUM_PPN -pre "$i" > /dev/null 2> /dev/null  # >> "$SEQ".BS.stats.log 2>&1
    cat "$i".boottrees > all.boottrees

    while true
    do  
        if [ $(wc -l all.boottrees | cut -d ' ' -f 1) -gt 499 ]; then
            break
        fi

        let i++
        ./iqtree -s "$SEQ" -bo 10 -st DNA -AICc -nt $PBS_NUM_PPN -pre "$i" > /dev/null 2> /dev/null
        cat "$i".boottrees >> all.boottrees
    done



#################
###    ASTRAL

# the same issue as previously, dataset contains missing samples, so we have to ensure, that the ASTRAL namemapfile is consistent with samples from fna file.

touch namemapfile
toPrint=""
if grep -q carniolicum "$SEQ"; then toPrint=$toPrint"carniolicum2x_carniolicum:"; fi
if grep -q  Podjelak2 "$SEQ"; then toPrint=$toPrint"Podjelak2_carniolicum,"; fi
if grep -q  Svibno7 "$SEQ"; then toPrint=$toPrint"Svibno7_carniolicum,"; fi
echo $toPrint >> namemapfile
toPrint=""
if grep -q crassistylum "$SEQ"; then toPrint=$toPrint"crassistylum2x_crassistylum:"; fi
if grep -q  652RDU2 "$SEQ"; then toPrint=$toPrint"652RDU2_crassistylum,"; fi
echo $toPrint >> namemapfile
toPrint=""
if grep -q cuspidatum "$SEQ"; then toPrint=$toPrint"cuspidatum2x_cuspidatum:"; fi
if grep -q  KONJ14 "$SEQ"; then toPrint=$toPrint"KONJ14_cuspidatum,"; fi
if grep -q  SVI9s "$SEQ"; then toPrint=$toPrint"SVI9s_cuspidatum,"; fi
echo $toPrint >> namemapfile
toPrint=""

if grep -q kuemmerlei "$SEQ"; then toPrint=$toPrint"kuemmerlei2x_kuemmerlei:"; fi
if grep -q  PRST2_kuemmerlei "$SEQ"; then toPrint=$toPrint"PRST2_kuemmerlei,"; fi
if grep -q  SHTRs3_kuemmerlei "$SEQ"; then toPrint=$toPrint"SHTRs3_kuemmerlei,"; fi
echo $toPrint >> namemapfile
toPrint=""

if grep -q linariifolium "$SEQ"; then toPrint=$toPrint"linariifolium2x_linariifolium:"; fi
if grep -q  DUCA8_linariifolium "$SEQ"; then toPrint=$toPrint"DUCA8_linariifolium,"; fi
if grep -q  Kavlar1lin_linariifolium "$SEQ"; then toPrint=$toPrint"Kavlar1lin_linariifolium,"; fi
echo $toPrint >> namemapfile
toPrint=""

if grep -q malcolmia "$SEQ"; then toPrint=$toPrint"malcolmia2x_malcolmia:"; fi
if grep -q  MalcolmiaI1_malcolmia "$SEQ"; then toPrint=$toPrint"MalcolmiaI1_malcolmia,"; fi
if grep -q  MalcolmiaII1_malcolmia "$SEQ"; then toPrint=$toPrint"MalcolmiaII1_malcolmia,"; fi
echo $toPrint >> namemapfile
toPrint=""

if grep -q microstylum "$SEQ"; then toPrint=$toPrint"microstylum2x_microstylum:"; fi
if grep -q  Alona1_microstylum "$SEQ"; then toPrint=$toPrint"Alona1_microstylum,"; fi
echo $toPrint >> namemapfile
toPrint=""

if grep -q pectinatum "$SEQ"; then toPrint=$toPrint"pectinatum2x_pectinatum:"; fi
if grep -q  PTER1_pectinatum "$SEQ"; then toPrint=$toPrint"PTER1_pectinatum,"; fi
if grep -q  TGTS1_pectinatum "$SEQ"; then toPrint=$toPrint"TGTS1_pectinatum,"; fi
echo $toPrint >> namemapfile
toPrint=""

if grep -q saxosum "$SEQ"; then toPrint=$toPrint"saxosum2x_saxosum:"; fi
if grep -q  BUILA1s_saxosum "$SEQ"; then toPrint=$toPrint"BUILA1s_saxosum,"; fi
if grep -q  CAKR14_saxosum "$SEQ"; then toPrint=$toPrint"CAKR14_saxosum,"; fi
echo $toPrint >> namemapfile
toPrint=""

if grep -q sylvestre "$SEQ"; then toPrint=$toPrint"sylvestre2x_sylvestre:"; fi
if grep -q  Dovlici1_sylvestre "$SEQ"; then toPrint=$toPrint"Dovlici1_sylvestre,"; fi
if grep -q  Sabotin11_sylvestre "$SEQ"; then toPrint=$toPrint"Sabotin11_sylvestre,"; fi
echo $toPrint >> namemapfile
toPrint=""

if grep -q vitekii "$SEQ"; then toPrint=$toPrint"vitekii2x_vitekii:"; fi
if grep -q  CRKV3_vitekii "$SEQ"; then toPrint=$toPrint"CRKV3_vitekii,"; fi
if grep -q  Neves1_vitekii "$SEQ"; then toPrint=$toPrint"Neves1_vitekii,"; fi
echo $toPrint >> namemapfile
toPrint=""

if grep -q witmannii "$SEQ"; then toPrint=$toPrint"witmannii2x_witmannii:"; fi
if grep -q  GEMHs1_witmannii "$SEQ"; then toPrint=$toPrint"GEMHs1_witmannii,"; fi
if grep -q  LESc_witmannii "$SEQ"; then toPrint=$toPrint"LESc_witmannii,"; fi
if grep -q  NASZ8_witmannii "$SEQ"; then toPrint=$toPrint"NASZ8_witmannii,"; fi
if grep -q  SUGA1sm_witmannii "$SEQ"; then toPrint=$toPrint"SUGA1sm_witmannii,"; fi
if grep -q  TUN2s_witmannii "$SEQ"; then toPrint=$toPrint"TUN2s_witmannii,"; fi
if grep -q  ZARN6s_witmannii "$SEQ"; then toPrint=$toPrint"ZARN6s_witmannii,"; fi
echo $toPrint >> namemapfile
toPrint=""



echo "$NAME":"$NAME"_"$PREVCLAF" >> namemapfile



while read TREE; do
echo $TREE > astralInputTree
java -jar ~/bin/astral.5.7.8/astral.5.7.8.jar -i astralInputTree -o TREE.astralTree --namemapfile namemapfile -t 0
cat TREE.astralTree >> astral.bsTrees
done < all.boottrees


sed -i 's/_apenninaxamporitana/_apxamp:001/g' astral.bsTrees
sed -i 's/_pratensis7x/_prat7x:001/g' astral.bsTrees
sed -i 's/_pratensisPyrenees/_pratPyrenees:001/g' astral.bsTrees

sed -i 's/carniolicum2x_carniolicum/carniolicum2x_carniolicum:001/' astral.bsTrees
sed -i 's/crassistylum2x_crassistylum/crassistylum2x_crassistylum:001/' astral.bsTrees
sed -i 's/cuspidatum2x_cuspidatum/cuspidatum2x_cuspidatum:001/' astral.bsTrees
sed -i 's/kuemmerlei2x_kuemmerlei/kuemmerlei2x_kuemmerlei:001/' astral.bsTrees
sed -i 's/linariifolium2x_linariifolium/linariifolium2x_linariifolium:001/' astral.bsTrees
sed -i 's/malcolmia2x_malcolmia/malcolmia2x_malcolmia:001/' astral.bsTrees
sed -i 's/microstylum2x_microstylum/microstylum2x_microstylum:001/' astral.bsTrees
sed -i 's/pectinatum2x_pectinatum/pectinatum2x_pectinatum:001/' astral.bsTrees
sed -i 's/saxosum2x_saxosum/saxosum2x_saxosum:001/' astral.bsTrees
sed -i 's/sylvestre2x_sylvestre/sylvestre2x_sylvestre:001/' astral.bsTrees
sed -i 's/vitekii2x_vitekii/vitekii2x_vitekii:001/' astral.bsTrees
sed -i 's/witmannii2x_witmannii/witmannii2x_witmannii:001/' astral.bsTrees

sed -i 's/CIUC6R.h1_odoratum22chrom/CIUC6R.h1_odoratum22chrom:001/' astral.bsTrees
sed -i 's/CIUC6R.h2_odoratum22chrom/CIUC6R.h2_odoratum22chrom:001/' astral.bsTrees
sed -i 's/CIUC6R.h3_odoratum22chrom/CIUC6R.h3_odoratum22chrom:001/' astral.bsTrees
sed -i 's/CIUC6R.h4_odoratum22chrom/CIUC6R.h4_odoratum22chrom:001/' astral.bsTrees
sed -i 's/Maliscak9.h1_odoratum22chrom/Maliscak9.h1_odoratum22chrom:001/' astral.bsTrees
sed -i 's/Maliscak9.h2_odoratum22chrom/Maliscak9.h2_odoratum22chrom:001/' astral.bsTrees
sed -i 's/Maliscak9.h3_odoratum22chrom/Maliscak9.h3_odoratum22chrom:001/' astral.bsTrees
sed -i 's/Maliscak9.h4_odoratum22chrom/Maliscak9.h4_odoratum22chrom:001/' astral.bsTrees
sed -i 's/BREITs1.h1_odoratum4x/BREITs1.h1_odoratum4x:001/' astral.bsTrees
sed -i 's/BREITs1.h2_odoratum4x/BREITs1.h2_odoratum4x:001/' astral.bsTrees
sed -i 's/BREITs1.h3_odoratum4x/BREITs1.h3_odoratum4x:001/' astral.bsTrees
sed -i 's/BREITs1.h4_odoratum4x/BREITs1.h4_odoratum4x:001/' astral.bsTrees
sed -i 's/Kavlar2.h1_odoratum4x/Kavlar2.h1_odoratum4x:001/' astral.bsTrees
sed -i 's/Kavlar2.h2_odoratum4x/Kavlar2.h2_odoratum4x:001/' astral.bsTrees
sed -i 's/Kavlar2.h3_odoratum4x/Kavlar2.h3_odoratum4x:001/' astral.bsTrees
sed -i 's/Kavlar2.h4_odoratum4x/Kavlar2.h4_odoratum4x:001/' astral.bsTrees
sed -i 's/GLAVs1.h1_odoratum6x/GLAVs1.h1_odoratum6x:001/' astral.bsTrees
sed -i 's/GLAVs1.h2_odoratum6x/GLAVs1.h2_odoratum6x:001/' astral.bsTrees
sed -i 's/GLAVs1.h3_odoratum6x/GLAVs1.h3_odoratum6x:001/' astral.bsTrees
sed -i 's/GLAVs1.h4_odoratum6x/GLAVs1.h4_odoratum6x:001/' astral.bsTrees
sed -i 's/GLAVs1.h5_odoratum6x/GLAVs1.h5_odoratum6x:001/' astral.bsTrees
sed -i 's/GLAVs1.h6_odoratum6x/GLAVs1.h6_odoratum6x:001/' astral.bsTrees
sed -i 's/VGrad1.h1_odoratum6x/VGrad1.h1_odoratum6x:001/' astral.bsTrees
sed -i 's/VGrad1.h2_odoratum6x/VGrad1.h2_odoratum6x:001/' astral.bsTrees
sed -i 's/VGrad1.h3_odoratum6x/VGrad1.h3_odoratum6x:001/' astral.bsTrees
sed -i 's/VGrad1.h4_odoratum6x/VGrad1.h4_odoratum6x:001/' astral.bsTrees
sed -i 's/VGrad1.h5_odoratum6x/VGrad1.h5_odoratum6x:001/' astral.bsTrees
sed -i 's/VGrad1.h6_odoratum6x/VGrad1.h6_odoratum6x:001/' astral.bsTrees
sed -i 's/RET1s.h1_odoratumRetezat/RET1s.h1_odoratumRetezat:001/' astral.bsTrees
sed -i 's/RET1s.h2_odoratumRetezat/RET1s.h2_odoratumRetezat:001/' astral.bsTrees
sed -i 's/RET1s.h3_odoratumRetezat/RET1s.h3_odoratumRetezat:001/' astral.bsTrees
sed -i 's/RET1s.h4_odoratumRetezat/RET1s.h4_odoratumRetezat:001/' astral.bsTrees
sed -i 's/SCO3.h1_odoratumRetezat/SCO3.h1_odoratumRetezat:001/' astral.bsTrees
sed -i 's/SCO3.h2_odoratumRetezat/SCO3.h2_odoratumRetezat:001/' astral.bsTrees
sed -i 's/SCO3.h3_odoratumRetezat/SCO3.h3_odoratumRetezat:001/' astral.bsTrees
sed -i 's/SCO3.h4_odoratumRetezat/SCO3.h4_odoratumRetezat:001/' astral.bsTrees




###    PhyloSD
#################


while read TREE; do
    echo $TREE > tree.bstree
    perl ./_reroot_tree.pl tree.bstree > tree.root.ph;
    sed -i 's/0.000000/0.000001/g' tree.root.ph;
    perl ./_check_lineages_polyploids.pl -t tree.root.ph >> log 
done < astral.bsTrees


grep "all$" log | head -n 1 > counts
grep "$NAME" log | grep -v ">" >> counts

head -n 501 counts > vybrane.counts

R CMD BATCH --no-save --no-restore PhyloSD.countBStrees.R count.log

mv t "$SEQ".counts   # that "t" is the name of R output.

mv counts "$SEQ".log



exit

