#!/bin/bash


# make PolyRelatedness input file from VCF

bcftools view -h Cacris_sl.vcf.gz | tail -n 1 | sed 's/.*FORMAT\t//' > Cacris_sl.genotypes

bcftools query -f '[ %GT]\t\n' Cacris_sl.vcf.gz >> Cacris_sl.genotypes

sed -i 's/^ //' Cacris_sl.genotypes
sed -i 's/\t/ /g' Cacris_sl.genotypes

# add header IND a POP


# transpose the file
awk '
{ 
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END {    
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str" "a[i,j];
        }
        print str
    }
}' Cacris_sl.genotypes > Cacris_sl.genotypes.t


sed -i 's/\///g' Cacris_sl.genotypes.t
sed -i 's/\./9/g' Cacris_sl.genotypes.t


# 43485 is number of genotypes
for f in {1..43485}; do
echo "L$f"
done >> l

# merge and modify according to the PolyRelatedness input file:

//configuration
//#alleledigits(1~4)	#outputdigits(0~10)	#missingallele	#ambiguousallele	#nthreads(1~64)
1	8	9	8	4
//genotype	
Ind Pop L1 L2 L3 L4 L5 L6 L7 L8
243OSA4 Acris_SS 00 00 00 01 00 11 00 11
243OSA7 Acris_SS 00 00 99 00 00 11 00 11
C156_110 Pindicola 00 00 00 00 00 11 00 11 
C192_106 C192 99 99 99 0000 0000 0000 0000 0000
C192_108 C192 0000 0000 0000 0000 0000 1111 0000 1111 99 
//end of file


sed -i 's/ /\t/g' Cacris_sl.genotypes.t.polyrelatedness


###################################################
# run PolyRelatedness

# "e" for estimating the relatedness;
# 1 = Huang et al. 2014 MOM estimator;
# 0 for between all individuals or 1 for within population;
# see PDF for details
PolyRelatedness.out \
Cacris_sl.genotypes.t.polyrelatedness \
Cacris_sl.genotypes.t.polyrelatedness.OUT.txt e 1 0


###################################################
# extract rows with polyploids - for violin plots
grep "C243" Cacris_sl.genotypes.t.polyrelatedness.OUT.txt > C243.txt
grep "C192" Cacris_sl.genotypes.t.polyrelatedness.OUT.txt > C192.txt
.......




