#!/bin/bash


# !!!! run this script in separate jobs
# while read FILE; do
# export FILE
#     # !!! SUBMIT individual jobs 
#     qsub
# done < list


export PBS_NUM_PPN   # number of cores
export ASTRALNAMEMAPFILE # namemapfile maps samples and their species


module add raxml-8.2.8
module add jdk-8


./PhyloNetworks/scripts/raxml.pl --seqdir="$FILE" --raxmldir=raxml."$FILE" --astraldir=astral."$FILE" --nodoastral --numCores "$PBS_NUM_PPN" --numboot 100


cp $ASTRALNAMEMAPFILE astral."$FILE"/astral.nameMapFile

######
##    if you allowed missing samples in sequences, astral.nameMapFile could contain a sample, that is missing in the tree, thus that sample has to be removed
cp astral."$FILE"/astral.nameMapFile astral."$FILE"/astral.nameMapFile.samples
sed -i 's/,/\n/g' astral."$FILE"/astral.nameMapFile.samples
sed -i 's/:/:\n/g' astral."$FILE"/astral.nameMapFile.samples
sed -i 's/.*://' astral."$FILE"/astral.nameMapFile.samples
sed -i '/^$/d' astral."$FILE"/astral.nameMapFile.samples

while read SAMPLE; do
if grep -q "$SAMPLE" $FILE/"$FILE".fasta.nex; then
    echo "$SAMPLE is present"
else 
    echo "$SAMPLE MISSING, removing from astral nameMapFile"
    sed -i "s/$SAMPLE//" astral."$FILE"/astral.nameMapFile
fi
done < astral."$FILE"/astral.nameMapFile.samples

sed -i 's/,,/,/g' astral."$FILE"/astral.nameMapFile
sed -i 's/:,/:/' astral."$FILE"/astral.nameMapFile
sed -i 's/,$//' astral."$FILE"/astral.nameMapFile
sed -i '/:$/d' astral."$FILE"/astral.nameMapFile
sed -i 's/,,/,/g' astral."$FILE"/astral.nameMapFile
sed -i 's/,,/,/g' astral."$FILE"/astral.nameMapFile
sed -i 's/,,/,/g' astral."$FILE"/astral.nameMapFile
sed -i 's/,,/,/g' astral."$FILE"/astral.nameMapFile
sed -i 's/,,/,/g' astral."$FILE"/astral.nameMapFile
sed -i 's/,,/,/g' astral."$FILE"/astral.nameMapFile
sed -i 's/,,/,/g' astral."$FILE"/astral.nameMapFile
sed -i 's/,,/,/g' astral."$FILE"/astral.nameMapFile
sed -i 's/,,/,/g' astral."$FILE"/astral.nameMapFile
sed -i 's/,,/,/g' astral."$FILE"/astral.nameMapFile
sed -i 's/,,/,/g' astral."$FILE"/astral.nameMapFile
sed -i '/:,$/d' astral."$FILE"/astral.nameMapFile
sed -i 's/:,/:/g' astral."$FILE"/astral.nameMapFile

# done, astral.nameMapFile matches samples in sequences


java -jar ./PhyloNetworks/Astral_binary/astral.5.7.8.jar \
-i raxml."$FILE"/besttrees.tre -b astral."$FILE"/BSlistfiles \
-r 100 \
-o astral."$FILE"/astral.tre \
-a astral."$FILE"/astral.nameMapFile


exit

