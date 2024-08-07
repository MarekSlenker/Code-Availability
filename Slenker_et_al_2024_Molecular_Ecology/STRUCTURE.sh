#!/bin/bash




export NUMINDS="196" # number of individuals
export L="3869"      # number of markers
export MAINPARAMS    # mainparams file
export EXTRAPARAMS   # extraparams file

Ks=$(seq 1 1 7)   # inferred number of populations (groups)



for dataset in *.STRUCTURE; do
    export dataset

    for K in $Ks; do
        export K
        
        SEED=$(($RANDOM*$RANDOM))
        export SEED
        
        echo "submitting $dataset with K $K and SEED $SEED"
        
        # qsub  !!!!!!  do the following command in separate jobs !!!!!
        structure -K "$k" -L "$L" -N "$NUMINDS" -m "$MAINPARAMS" -e "$EXTRAPARAMS" -D $SEED -i "$dataset" -o ./${dataset%.*}_K$K.$SEED | tee ./${dataset%.*}_K$K.log
        
    done
done


exit

