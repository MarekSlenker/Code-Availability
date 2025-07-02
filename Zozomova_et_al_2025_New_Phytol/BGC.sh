#!/bin/bash



prefix="RivMatHybr"


cd $SCRATCHDIR

cp $DATADIR/* .


module add python36-modules-gcc
module add hdf5/1.12.2-gcc-10.2.1-gfdwqr3


mcmc=100000   # Number of MCMC generations (integer)
burnin=50000  # MCMC burnin generations (integer)
thin=24       # Thin MCMC to every $thin generations (integer)
gt_uncert=0   # gt_uncert: 0 if not using gt_uncert, 1 if using it (boolean integer)
# error=0.001   # error: Error rate for genotype_uncertainty model (float)
linkage=0     # Linkage: 0 for no linkage model, 1 for linkage model (boolean integer)
# linkage_dist=0.5 # Maximum distance between loci, free recombination. Only used for linkage model (default=0.5)
stddev_gamma=0.05  # standard deviation for Gaussian proposal of cline parameter gamma (default=0.05)
stddev_zeta=0.05   # standard deviation for Gaussian proposal of cline parameter zeta (default=0.05)
stddev_eta_kappa=0.02  # standard deviation for Gaussian proposal of cline parameters eta and kappa (default=0.02)
mcmc_tuning=0.1   # MCMC tuning parameter, maximum deviate from uniform for proposed hybrid index  (default=0.1)

<PATH>/bin/BGC-Bayesian-genomic-clines/bgc \
	-a ${prefix}_p0in.txt \
		-b ${prefix}_p1in.txt \
		-h ${prefix}_admixedin.txt \
		-O 0 -x $mcmc -n $burnin -t $thin -N $gt_uncert \
		-E $error -m $linkage -q 1 -I 1 -p 1 \
		-F ${prefix}_mcmcout_${run} \
		-g $stddev_gamma -z $stddev_zeta \
		-e $stddev_eta_kappa -u $mcmc_tuning;


<PATH>/bin/BGC-Bayesian-genomic-clines/estpost -i ${prefix}_mcmcout_${run}.hdf5 \
	-p LnL -o ${prefix}_stat_lnl_${run} -s 2 -w 0

<PATH>/bin/BGC-Bayesian-genomic-clines/estpost -i ${prefix}_mcmcout_${run}.hdf5 \
	-p alpha -o ${prefix}_stat_a0_${run} -s 2 -w 0

<PATH>/bin/BGC-Bayesian-genomic-clines/estpost -i ${prefix}_mcmcout_${run}.hdf5 \
	-p beta -o ${prefix}_stat_b0_${run} -s 2 -w 0

<PATH>/bin/BGC-Bayesian-genomic-clines/estpost -i ${prefix}_mcmcout_${run}.hdf5 \
	-p hi -o ${prefix}_stat_hi_${run} -s 2 -w 0

<PATH>/bin/BGC-Bayesian-genomic-clines/estpost -i ${prefix}_mcmcout_${run}.hdf5 \
	-p gamma-quantile -o ${prefix}_stat_qa_${run} -s 2 -w 0

<PATH>/bin/BGC-Bayesian-genomic-clines/estpost -i ${prefix}_mcmcout_${run}.hdf5 \
	-p zeta-quantile -o ${prefix}_stat_qb_${run} -s 2 -w 0



cp -r $SCRATCHDIR $DATADIR

