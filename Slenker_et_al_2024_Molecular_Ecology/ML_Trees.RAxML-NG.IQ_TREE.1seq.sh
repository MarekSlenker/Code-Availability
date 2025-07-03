#!/bin/bash


export ALIGNMENT # input fasta sequence
export BOOTSTRAP="500"  # number of BS trees
export PBS_NUM_PPN="10" # number of CPU threads



echo "Loading modules"
source /cvmfs/software.metacentrum.cz/modulefiles/5.1.0/loadmodules
module add raxml-ng-8 
module add iqtree-1.6.12
module add R-4.0.0-gcc



# Going to working directory $SCRATCHDIR
cd "$SCRATCHDIR"/

cp "$ALIGNMENT" .

#======================================================================
#  iqtree model finder

iqtree -s "$ALIGNMENT" -m TESTONLY -st DNA -nt $PBS_NUM_PPN -pre iqtreeoutput -rcluster 10 -redo -quiet 


# find the best model
MODEL=$(grep "Best-fit model according to BIC:" iqtreeoutput.iqtree | cut -d ' ' -f 6)

echo "the best model: $MODEL" >> "${ALIGNMENT%.*}".log

# parse IQTREE output to RAxML input
NEWMODEL=$( ./ML_Trees.PhyhlogenyModelParser.sh $MODEL )

echo "the NEW model: $NEWMODEL" >> "${ALIGNMENT%.*}".log

[ -z "$NEWMODEL" ] && exit 1

#======================================================================
#  RAxML-NG 

raxml-ng --msa "$ALIGNMENT" --model "$NEWMODEL" --tree pars{5},rand{5} --blopt nr_safe --threads $PBS_NUM_PPN --prefix ${ALIGNMENT%.*} >> "${ALIGNMENT%.*}".log || ./raxml-ng --msa "$ALIGNMENT" --model "$NEWMODEL" --tree pars{5},rand{5} --blopt nr_safe --threads 1 --prefix ${ALIGNMENT%.*} >> "${ALIGNMENT%.*}".log


if [ $BOOTSTRAP -gt 1 ]; then
    for i in $(seq 1 "$PBS_NUM_PPN"); do
        raxml-ng --bootstrap --msa "$ALIGNMENT" --model ./*.raxml.bestModel --bs-trees $(($BOOTSTRAP / $PBS_NUM_PPN)) --blopt nr_safe --seed "$RANDOM" --threads 1 --prefix "${ALIGNMENT%.*}"_"$i" >> "${ALIGNMENT%.*}"_"$i".log &
        done
    wait
    cat *.raxml.bootstraps > allbootstraps.bootstraps  # merge all partial BSs
    
    ./raxml-ng --support --tree "${ALIGNMENT%.*}".raxml.bestTree --bs-trees allbootstraps.bootstraps --threads 1 --prefix "${ALIGNMENT%.*}"  >> "${ALIGNMENT%.*}".log
fi


exit
