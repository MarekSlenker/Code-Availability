#!/bin/bash

export ALIGNMENT    # PHYLIP format, without invariant sites
export PBS_NUM_PPN  # number of CPU threads
# GTR+ASC_FELS{176698}  = Felsenstein's method with total number of invariable sites w


source /cvmfs/software.metacentrum.cz/modulefiles/5.1.0/loadmodules
module load raxml-ng/8


raxml-ng \
--msa "$ALIGNMENT" \
--model "GTR+ASC_FELS{176698}" \
--tree pars{10},rand{10} \
--blopt nr_safe \
--threads $PBS_NUM_PPN \
--prefix ${ALIGNMENT%.*}_FELS

exit




