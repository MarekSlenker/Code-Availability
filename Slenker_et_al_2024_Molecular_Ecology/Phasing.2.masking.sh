#!/bin/bash

SAMPLEPLOIDYLIST="samplePloidyList"

while read SAMPLE; do
    read PLOIDY
    echo "processing $SAMPLE"
    
    for((i=1;i<=$PLOIDY;i=i+1)); do
        echo "$i"
        python3 haplonerate_MK_bere1vrati1alelu-3x,5x,7x...py --edit ref ./$SAMPLE.whatshap.gtf ./$SAMPLE.v$i.phased.fasta --reference ./$SAMPLE.fasta.Nmasked.fasta --allele $i --output ./$SAMPLE.phased.haplonerate.Nmasked.$i.fasta
    done

done < $SAMPLEPLOIDYLIST




for i in $(ls -1 ./*.fasta);
do
    echo $i
    #Removes line breaks from fasta file
    awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' "$i" > temp_file
    
    # # Sort by NAME
    awk '{idx=index($0, " "); ID=$1; NAME=substr($0, idx+1); getline SEQ; printf("%s\t%s\t%s\n", NAME, ID, SEQ);}' temp_file | sort -k1d,1 | awk -F"\t" '{printf("%s \n%s\n", $2,$3)}' > "$i"
    
done




# merge phased and unPhased sequences
while read SAMPLE; do
    read PLOIDY
    
    echo "processing $SAMPLE"
    
    for((i=1;i<=$PLOIDY;i=i+1)); do
        echo "$i"
        cat "$SAMPLE".unPhased.fasta >> "$SAMPLE".phased.haplonerate.Nmasked."$i".fasta
    done
    
done < $SAMPLEPLOIDYLIST



# EXTRACT



read SAMPLE1 < $SAMPLEPLOIDYLIST

fastaLength=$(grep -c ">" "$SAMPLE1".phased.haplonerate.Nmasked.1.fasta)
doubleFastaLength=$((fastaLength+fastaLength))
echo $doubleFastaLength

mkdir RESDIR



while read SAMPLE; do
    read PLOIDY

    for((pld=1;pld<=$PLOIDY;pld++)); do  # pre kazdy jeden variant = ploidia
        echo "processing $SAMPLE, pld $pld/$PLOIDY"
        
        for((line=1;line<="$doubleFastaLength";line++)); do # pre kazdy gen vzorky = berem 2 riadky naraz
            
            exonName=$(sed "${line}q;d" "$SAMPLE".phased.haplonerate.Nmasked."$pld".fasta | sed 's/>//' | sed 's/ //')
            
            echo ">$SAMPLE-h$pld" >> ./RESDIR/"$exonName".Nmasked.fasta
            
            let "line++"
            
            echo $(sed "${line}q;d" "$SAMPLE".phased.haplonerate.Nmasked."$pld".fasta) >> ./RESDIR/"$exonName".Nmasked.fasta
        done
        
    done
    
done < $SAMPLEPLOIDYLIST




# "Phasing statistics"


printf "Individual\tphased\tunPhased\n" > Phasing_statistics.txt
while read SAMPLE; do
    read PLOIDY
    
    { printf '%s\t' "$SAMPLE"
        echo "$(($(grep -c ">" "$SAMPLE".phased.haplonerate.Nmasked.1.fasta)/2))" | xargs printf '%s\t%s'
        echo "$(grep -c ">" "$SAMPLE".unPhased.fasta)" | xargs printf '%s\t\n%s'
                
    } >> Phasing_statistics.txt    
done < $SAMPLEPLOIDYLIST
