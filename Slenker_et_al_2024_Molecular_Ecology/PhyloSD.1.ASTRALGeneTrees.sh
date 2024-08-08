#!/bin/bash

# Input trees are bestTrees made by RAxML. We used consensus sequences of 2x samples, and phased sequences for polyploids. 
# however, dataset contains missing samples, so we have to check presence of each sample in tree

for f in *bestTree; do
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

if grep -q acraC018 "$f"; then
echo "acraC018.101.h1_acrisPP:acraC018.101.h1_acrisPP" >> namemapfile
echo "acraC018.101.h2_acrisPP:acraC018.101.h2_acrisPP" >> namemapfile
echo "acraC018.101.h3_acrisPP:acraC018.101.h3_acrisPP" >> namemapfile
echo "acraC018.101.h4_acrisPP:acraC018.101.h4_acrisPP" >> namemapfile
fi

if grep -q acraC018 "$f"; then
echo "acraC095.109.h1_acrisPP:acraC095.109.h1_acrisPP" >> namemapfile
echo "acraC095.109.h2_acrisPP:acraC095.109.h2_acrisPP" >> namemapfile
echo "acraC095.109.h3_acrisPP:acraC095.109.h3_acrisPP" >> namemapfile
fi

if grep -q acraC018 "$f"; then
echo "acraC192.3.h1_acrisPP:acraC192.3.h1_acrisPP" >> namemapfile
echo "acraC192.3.h2_acrisPP:acraC192.3.h2_acrisPP" >> namemapfile
echo "acraC192.3.h3_acrisPP:acraC192.3.h3_acrisPP" >> namemapfile
echo "acraC192.3.h4_acrisPP:acraC192.3.h4_acrisPP" >> namemapfile
fi

java -jar ~/bin/astral.5.7.8/astral.5.7.8.jar -i $f -o "$f".astralTree --namemapfile namemapfile -t 0
rm namemapfile; 
done

sed -i 's/_EBalkan/_EBalkan:001/g' *astralTree
sed -i 's/_Dinaric/_Dinaric:001/g' *astralTree
sed -i 's/_SharGramos/_SharGramos:001/g' *astralTree
sed -i 's/_vardousPindicola/_vardousPindicola:001/g' *astralTree
sed -i 's/_lazica/_lazica:001/g' *astralTree
sed -i 's/_matthioli/_matthioli:001/g' *astralTree
sed -i 's/_rivularis/_rivularis:001/g' *astralTree
sed -i 's/_anatolica/_anatolica:001/g' *astralTree
sed -i 's/_acrisPP/_acrisPP:001/g' *astralTree

exit
