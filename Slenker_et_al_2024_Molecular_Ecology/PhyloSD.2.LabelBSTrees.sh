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
if grep -q EBalkan "$f"; then toPrint=$toPrint"EBalkan_EBalkan:"; fi
if grep -q  acraBAB6_EBalkan "$f"; then toPrint=$toPrint"acraBAB6_EBalkan,"; fi
if grep -q  acraC003.104_EBalkan "$f"; then toPrint=$toPrint"acraC003.104_EBalkan,"; fi
if grep -q  acraC013.101_EBalkan "$f"; then toPrint=$toPrint"acraC013.101_EBalkan,"; fi
if grep -q  acraC019.103_EBalkan "$f"; then toPrint=$toPrint"acraC019.103_EBalkan,"; fi
if grep -q  acraC149.8_EBalkan "$f"; then toPrint=$toPrint"acraC149.8_EBalkan,"; fi
if grep -q  acraC157.104_EBalkan "$f"; then toPrint=$toPrint"acraC157.104_EBalkan,"; fi
if grep -q  acraC161.104_EBalkan "$f"; then toPrint=$toPrint"acraC161.104_EBalkan,"; fi
if grep -q  acraGOL1_EBalkan "$f"; then toPrint=$toPrint"acraGOL1_EBalkan,"; fi
if grep -q  acraOSA10_EBalkan "$f"; then toPrint=$toPrint"acraOSA10_EBalkan,"; fi
echo $toPrint >> namemapfile

toPrint=""
if grep -q Dinaric "$f"; then toPrint=$toPrint"Dinaric_Dinaric:"; fi
if grep -q  acraC015.107_Dinaric "$f"; then toPrint=$toPrint"acraC015.107_Dinaric,"; fi
if grep -q  acraC016.103_Dinaric "$f"; then toPrint=$toPrint"acraC016.103_Dinaric,"; fi
echo $toPrint >> namemapfile

toPrint=""
if grep -q SharGramos "$f"; then toPrint=$toPrint"SharGramos_SharGramos:"; fi
if grep -q  acraC153.101_SharGramos "$f"; then toPrint=$toPrint"acraC153.101_SharGramos,"; fi
if grep -q  acraC191.10_SharGramos "$f"; then toPrint=$toPrint"acraC191.10_SharGramos,"; fi
if grep -q  acraC193.15_SharGramos "$f"; then toPrint=$toPrint"acraC193.15_SharGramos,"; fi
echo $toPrint >> namemapfile

toPrint=""
if grep -q vardousPindicola "$f"; then toPrint=$toPrint"vardousPindicola_vardousPindicola:"; fi
if grep -q  acrpC012.103_vardousPindicola "$f"; then toPrint=$toPrint"acrpC012.103_vardousPindicola,"; fi
if grep -q  acrpC154.110_vardousPindicola "$f"; then toPrint=$toPrint"acrpC154.110_vardousPindicola,"; fi
if grep -q  acrvC004.103_vardousPindicola "$f"; then toPrint=$toPrint"acrvC004.103_vardousPindicola,"; fi
echo $toPrint >> namemapfile

toPrint=""
if grep -q lazica "$f"; then toPrint=$toPrint"lazica_lazica:"; fi
if grep -q  lazTRI4_lazica "$f"; then toPrint=$toPrint"lazTRI4_lazica,"; fi
if grep -q  lazTRK2_lazica "$f"; then toPrint=$toPrint"lazTRK2_lazica,"; fi
echo $toPrint >> namemapfile

toPrint=""
if grep -q matthioli "$f"; then toPrint=$toPrint"matthioli_matthioli:"; fi
if grep -q  matC155.102_matthioli "$f"; then toPrint=$toPrint"matC155.102_matthioli,"; fi
if grep -q  matGRM9_matthioli "$f"; then toPrint=$toPrint"matGRM9_matthioli,"; fi
echo $toPrint >> namemapfile

toPrint=""
if grep -q rivularis "$f"; then toPrint=$toPrint"rivularis_rivularis:"; fi
if grep -q  rivALE4_rivularis "$f"; then toPrint=$toPrint"rivALE4_rivularis,"; fi
if grep -q  rivCAP9_rivularis "$f"; then toPrint=$toPrint"rivCAP9_rivularis,"; fi
if grep -q  rivPOB10_rivularis "$f"; then toPrint=$toPrint"rivPOB10_rivularis,"; fi
echo $toPrint >> namemapfile

toPrint=""
if grep -q anatolica "$f"; then toPrint=$toPrint"anatolica_anatolica:"; fi
if grep -q  ulDEL1_anatolica "$f"; then toPrint=$toPrint"ulDEL1_anatolica,"; fi
if grep -q  ulDEL3_anatolica "$f"; then toPrint=$toPrint"ulDEL3_anatolica,"; fi
if grep -q  ulUD4_anatolica "$f"; then toPrint=$toPrint"ulUD4_anatolica,"; fi
if grep -q  ulUD9_anatolica "$f"; then toPrint=$toPrint"ulUD9_anatolica,"; fi
if grep -q  ulULU3_anatolica "$f"; then toPrint=$toPrint"ulULU3_anatolica,"; fi
echo $toPrint >> namemapfile

echo "$NAME":"$NAME"_"$PREVCLAF" >> namemapfile




while read TREE; do
echo $TREE > astralInputTree
java -jar ~/bin/astral.5.7.8/astral.5.7.8.jar -i astralInputTree -o TREE.astralTree --namemapfile namemapfile -t 0
cat TREE.astralTree >> astral.bsTrees
done < all.boottrees


sed -i 's/_EBalkan/_EBalkan:001/g' astral.bsTrees
sed -i 's/_Dinaric/_Dinaric:001/g' astral.bsTrees
sed -i 's/_SharGramos/_SharGramos:001/g' astral.bsTrees
sed -i 's/_vardousPindicola/_vardousPindicola:001/g' astral.bsTrees
sed -i 's/_lazica/_lazica:001/g' astral.bsTrees
sed -i 's/_matthioli/_matthioli:001/g' astral.bsTrees
sed -i 's/_rivularis/_rivularis:001/g' astral.bsTrees
sed -i 's/_anatolica/_anatolica:001/g' astral.bsTrees
sed -i 's/_acrisPP/_acrisPP:001/g' astral.bsTrees




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

mv t "$SEQ".counts   # that "t" output from R..

mv counts "$SEQ".log



exit

