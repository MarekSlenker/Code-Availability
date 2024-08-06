#!/bin/bash




export DATADIR   # Fasta sequences to phase
export DEDUPDIR  # fastq reads
export SAMPLE   # sample name
export PLOIDY   # sample ploidy



echo "Loading modules"
module add python-3.6.2-gcc
module add parallel-20160622
module add gatk-3.7
module load gatk-4.0.1.0
module add picard-2.22.1
module add bwa-0.7.17
module add samtools-1.10
module add bcftools-1.6


cd $SCRATCHDIR  # working directory

mkdir inFiles
# Copy needed data
cp "$DATADIR"/* inFiles
cp "$DEDUPDIR"/"$SAMPLE".* . && parallel -X bunzip2 -v ::: *.bz2 || exit 1
# copy scripts
cp Phasing.map_to_supercontigs_MK.sh .
cp Phasing.phase_alleles_with_WhatsHap_MK.sh .
cp Phasing.extract_phase_bcftools_MK.sh .



#####STEP ZERO

# Collect all sequences of SAMPLE
for i in $(ls -1 ./inFiles/*.fasta);
do
    awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' $i > temp # convert multiline fasta to single line fasta

    sed 's/-//g' temp > "$i"
    
    echo $i >> "$SAMPLE".fasta
    grep -A1 "$SAMPLE" "$i" |grep -v "$SAMPLE"  >> "$SAMPLE".fasta
done

sed -i 's/.\/inFiles\//>/g' "$SAMPLE".fasta
sed -i 's/.fasta//g' "$SAMPLE".fasta


# Remove empty sequences
awk 'BEGIN {RS = ">" ; FS = "\n" ; ORS = ""} {if ($2) print ">"$0}' "$SAMPLE".fasta > tempfile
mv tempfile "$SAMPLE".fasta




#####STEP ONE: Map reads
echo "#####STEP ONE: Map reads"

bash ./Phasing.map_to_supercontigs_MK.sh $SAMPLE "$SAMPLE".R1.fq "$SAMPLE".R2.fq $PLOIDY 


#####STEP TWO: Phase alleles with WhatsHap
echo "#####STEP TWO: Phase alleles with WhatsHap"
bash ./Phasing.phase_alleles_with_WhatsHap_MK.sh $SAMPLE $PLOIDY 



#####STEP THREE: Extract phased seqs
echo "#####STEP THREE: Extract phased seqs"

# make genelist
cut -f 2 $SAMPLE.whatshap.stats.tsv | sed '/chromosome/d' | sed '/ALL/d' > genelist.txt


bash ./Phasing.extract_phase_bcftools_MK.sh $SAMPLE $PLOIDY genelist.txt 


# remove linebreaks
for i in $(ls -1 *.fasta);
do
    cat $i | awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' > "$i"_x
    mv "$i"_x "$i"
done

cp $SAMPLE.fasta $SAMPLE.notPhased.fasta
while read PATTERN; do
    sed -i -e "/$PATTERN/,+1d" $SAMPLE.notPhased.fasta
done < genelist.txt



exit

