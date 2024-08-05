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
bash ./phase_alleles_with_WhatsHap_MK.sh $SAMPLE $PLOIDY 



#####STEP THREE: Extract phased seqs
echo "#####STEP THREE: Extract phased seqs"

# make genelist
cut -f 2 $SAMPLE.whatshap.stats.tsv | sed '/chromosome/d' | sed '/ALL/d' > genelist.txt


bash ./extract_phase_bcftools_MK.sh $SAMPLE $PLOIDY genelist.txt 



# Kontroloval som fazovanie, fazuje to viac menej uspesne...


# !!!!!
# OK, netusim preco, niektore sekvencie su duplikovane. a niektore zase chybaju, pravdepodobne bez variability
# tie co chybaju, tie nemaju varianty. TODO: tieto doplnim po 1 seq
# potrebujem inlinovat a spravit zoznam tych co su > do novych suborov
# tieto subory, aby boli rovnake:
# $SAMPLE.whatshap.gtf
# ./$SAMPLE.v1.phased.fasta
# ./$SAMPLE.v2.phased.fasta
# $SAMPLE.fasta.iupac.fasta

# nazvy seq su v "genelist.txt"
# cat genelist.txt | uniq > uniq_genelist.txt

# vyhodim linebraks
for i in $(ls -1 *.fasta);
do
    cat $i | awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' > "$i"_x
    mv "$i"_x "$i"
done

cp $SAMPLE.fasta $SAMPLE.notPhased.fasta
while read PATTERN; do
    
    # teraz nezmysel, vsetky su OK uniq
    # for((i=1;i<=$PLOIDY;i++)); do
    #     grep -m 1 -A 1 $PATTERN ./$SAMPLE.v$i.phased.fasta >> ./$SAMPLE.v$i.phased.uniq.fasta;
    # done
    
    # grep -m 1 -A 1 $PATTERN ./$SAMPLE.fasta.iupac.fasta >> ./$SAMPLE.fasta.iupac.uniq.fasta;
    # grep -m 1 $PATTERN ./$SAMPLE.whatshap.gtf >> ./$SAMPLE.whatshap.uniq.gtf;
    # a tie ktore nejaku homology
    
    sed -i -e "/$PATTERN/,+1d" $SAMPLE.notPhased.fasta
    
    #done < uniq_genelist.txt
done < genelist.txt



## OK, mam vsetky rovnake


# Copy results back to storage
echo "Copying results back to $RESULTSDIR"

# cp *iupac* *Nmasked* *uniq* *whatshap* *phased* *notPhased* "$SAMPLE".fasta *snps.toIupac.vc* *snps.N.vc* "$RESULTSDIR"/ # || export CLEAN_SCRATCH='false'

cp *iupac* *Nmasked* *uniq* *whatshap* *phased* *notPhased* "$SAMPLE".fasta *.vcf "$RESULTSDIR"/ # || export CLEAN_SCRATCH='false'



# cp ./$SAMPLE.vsetkySeq.fasta ./$SAMPLE.phased.haplonerate.fasta *.iupac.fasta *.notPhased.fasta *uniq.fasta *.gtf "$RESULTSDIR"/ || export CLEAN_SCRATCH='false'


# Clean-up of SCRATCH
if [ "$CLEAN_SCRATCH" != "false" ]; then rm -rf $SCRATCHDIR/*; fi

echo

exit















# !!!!!!!!!!!
# UZ som len skopiroval data, s ktorymi budem dalej robit lokalne
#  TU NASLEDUJE vytahanie fazovanych s haplonerate, ale to funguje len pre 2x.   ja to spravim lokalne inak
#  OK, upravil som haplonerate, berie to len 2, ale polyploidy viem zbehnut vzdy po 2, cize 4x zbehneme dvakar dve pary.






#####STEP THREE: Haplonerate
echo "#####STEP THREE: Haplonerate"

python3.6 haplonerate.py --edit delete ./$SAMPLE.whatshap.uniq.gtf ./$SAMPLE.v1.phased.uniq.fasta ./$SAMPLE.v2.phased.uniq.fasta --reference ./$SAMPLE.fasta.iupac.uniq.fasta --output ./$SAMPLE.phased.haplonerate.fasta


# ok, mas vsetky haplotypy v "./$SAMPLE.phased.haplonerate.fasta" a tie co haplotypy nemaju su v "$SAMPLE.notPhased.fasta". teraz ich zlucis a mas vsetky seq za vzorku

cat ./$SAMPLE.phased.haplonerate.fasta > ./$SAMPLE.vsetkySeq.fasta
cat $SAMPLE.notPhased.fasta >> ./$SAMPLE.vsetkySeq.fasta


#Removes line breaks from fasta file
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' ./$SAMPLE.vsetkySeq.fasta  > temp_file

# Sort by NAME
awk '{idx=index($0, " "); ID=$1; NAME=substr($0, idx+1); getline SEQ; printf("%s\t%s\t%s\n", NAME, ID, SEQ);}' temp_file | sort -k1d,1 | awk -F"\t" '{printf("%s \n%s\n", $2,$3)}' > ./$SAMPLE.vsetkySeq.fasta





















# Copy results back to storage.   Results are compressd!
cp -a STR_RES/*_f "$RESULTSDIR"/ #|| export CLEAN_SCRATCH=false


# Clean-up of SCRATCH
if [ "$CLEAN_SCRATCH" != "false" ]; then rm -rf $SCRATCHDIR/*; fi

exit

