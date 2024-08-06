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



# while read SAMPLE; do
#     read PLOIDY
#     echo "processing $SAMPLE"

#     grep -A 1 --no-group-separator "_h1" "$SAMPLE".phased.haplonerate.iupac.1.fasta | sed 's/supercontig//; s/_h1//' > "$SAMPLE".phased.haplonerate.iupac.h1.fasta
#     grep -A 1 --no-group-separator "_h2" "$SAMPLE".phased.haplonerate.iupac.1.fasta | sed 's/supercontig//; s/_h2//' > "$SAMPLE".phased.haplonerate.iupac.h2.fasta

#     grep -A 1 --no-group-separator "_h1" "$SAMPLE".phased.haplonerate.iupac.3.fasta | sed 's/supercontig//; s/_h1//' > "$SAMPLE".phased.haplonerate.iupac.h3.fasta
#     grep -A 1 --no-group-separator "_h2" "$SAMPLE".phased.haplonerate.iupac.3.fasta | sed 's/supercontig//; s/_h2//' > "$SAMPLE".phased.haplonerate.iupac.h4.fasta

#     grep -A 1 --no-group-separator "_h1" "$SAMPLE".phased.haplonerate.iupac.5.fasta | sed 's/supercontig//; s/_h1//' > "$SAMPLE".phased.haplonerate.iupac.h5.fasta
#     grep -A 1 --no-group-separator "_h2" "$SAMPLE".phased.haplonerate.iupac.5.fasta | sed 's/supercontig//; s/_h2//' > "$SAMPLE".phased.haplonerate.iupac.h6.fasta
#     ###
#     grep -A 1 --no-group-separator "_h1" "$SAMPLE".phased.haplonerate.Nmasked.1.fasta | sed 's/supercontig//; s/_h1//' > "$SAMPLE".phased.haplonerate.Nmasked.h1.fasta
#     grep -A 1 --no-group-separator "_h2" "$SAMPLE".phased.haplonerate.Nmasked.1.fasta | sed 's/supercontig//; s/_h2//' > "$SAMPLE".phased.haplonerate.Nmasked.h2.fasta

#     grep -A 1 --no-group-separator "_h1" "$SAMPLE".phased.haplonerate.Nmasked.3.fasta | sed 's/supercontig//; s/_h1//' > "$SAMPLE".phased.haplonerate.Nmasked.h3.fasta
#     grep -A 1 --no-group-separator "_h2" "$SAMPLE".phased.haplonerate.Nmasked.3.fasta | sed 's/supercontig//; s/_h2//' > "$SAMPLE".phased.haplonerate.Nmasked.h4.fasta

#     grep -A 1 --no-group-separator "_h1" "$SAMPLE".phased.haplonerate.Nmasked.5.fasta | sed 's/supercontig//; s/_h1//' > "$SAMPLE".phased.haplonerate.Nmasked.h5.fasta
#     grep -A 1 --no-group-separator "_h2" "$SAMPLE".phased.haplonerate.Nmasked.5.fasta | sed 's/supercontig//; s/_h2//' > "$SAMPLE".phased.haplonerate.Nmasked.h6.fasta
# done < barbaroides_samplePloidyList_bezOutgroups.txt


while read SAMPLE; do
    read PLOIDY
    echo "processing $SAMPLE"
    
    # grep -A 1 --no-group-separator "_h1" "$SAMPLE".phased.haplonerate.iupac.1.fasta | sed 's/_h1//' > "$SAMPLE".phased.haplonerate.iupac.h1.fasta
    # grep -A 1 --no-group-separator "_h2" "$SAMPLE".phased.haplonerate.iupac.1.fasta | sed 's/_h2//' > "$SAMPLE".phased.haplonerate.iupac.h2.fasta
    
    # grep -A 1 --no-group-separator "_h1" "$SAMPLE".phased.haplonerate.iupac.3.fasta | sed 's/_h1//' > "$SAMPLE".phased.haplonerate.iupac.h3.fasta
    # grep -A 1 --no-group-separator "_h2" "$SAMPLE".phased.haplonerate.iupac.3.fasta | sed 's/_h2//' > "$SAMPLE".phased.haplonerate.iupac.h4.fasta
    
    # grep -A 1 --no-group-separator "_h1" "$SAMPLE".phased.haplonerate.iupac.5.fasta | sed 's/_h1//' > "$SAMPLE".phased.haplonerate.iupac.h5.fasta
    # grep -A 1 --no-group-separator "_h2" "$SAMPLE".phased.haplonerate.iupac.5.fasta | sed 's/_h2//' > "$SAMPLE".phased.haplonerate.iupac.h6.fasta
    ###
    # grep -A 1 --no-group-separator "_h1" "$SAMPLE".phased.haplonerate.Nmasked.1.fasta | sed 's/_h1//' > "$SAMPLE".phased.haplonerate.Nmasked.h1.fasta
    # grep -A 1 --no-group-separator "_h2" "$SAMPLE".phased.haplonerate.Nmasked.1.fasta | sed 's/_h2//' > "$SAMPLE".phased.haplonerate.Nmasked.h2.fasta
    
    # grep -A 1 --no-group-separator "_h1" "$SAMPLE".phased.haplonerate.Nmasked.3.fasta | sed 's/_h1//' > "$SAMPLE".phased.haplonerate.Nmasked.h3.fasta
    # grep -A 1 --no-group-separator "_h2" "$SAMPLE".phased.haplonerate.Nmasked.3.fasta | sed 's/_h2//' > "$SAMPLE".phased.haplonerate.Nmasked.h4.fasta
    
    # grep -A 1 --no-group-separator "_h1" "$SAMPLE".phased.haplonerate.Nmasked.5.fasta | sed 's/_h1//' > "$SAMPLE".phased.haplonerate.Nmasked.h5.fasta
    # grep -A 1 --no-group-separator "_h2" "$SAMPLE".phased.haplonerate.Nmasked.5.fasta | sed 's/_h2//' > "$SAMPLE".phased.haplonerate.Nmasked.h6.fasta
    
    ###
    
    # cat "$SAMPLE".v1.phased.fasta >> "$SAMPLE".phased.chimeric.h1.fasta
    # cat "$SAMPLE".v2.phased.fasta >> "$SAMPLE".phased.chimeric.h2.fasta
    # cat "$SAMPLE".v3.phased.fasta >> "$SAMPLE".phased.chimeric.h3.fasta
    # cat "$SAMPLE".v4.phased.fasta >> "$SAMPLE".phased.chimeric.h4.fasta
    # cat "$SAMPLE".v5.phased.fasta >> "$SAMPLE".phased.chimeric.h5.fasta
    # cat "$SAMPLE".v6.phased.fasta >> "$SAMPLE".phased.chimeric.h6.fasta

    
    # for((i=1;i<=$PLOIDY;i=i+1)); do
    #     echo "$i"
    #     cat "$SAMPLE".phased.haplonerate.Nmasked."$PLOIDY".fasta >> "$SAMPLE".phased.haplonerate.Nmasked.fasta
    # done

done < $SAMPLEPLOIDYLIST



while read SAMPLE; do
    read PLOIDY
    
    echo "processing $SAMPLE"
    
    # cat "$SAMPLE".notPhased.fasta >> "$SAMPLE".phased.haplonerate.iupac.h1.fasta
    # cat "$SAMPLE".notPhased.fasta >> "$SAMPLE".phased.haplonerate.iupac.h2.fasta
    # cat "$SAMPLE".notPhased.fasta >> "$SAMPLE".phased.haplonerate.iupac.h3.fasta
    # cat "$SAMPLE".notPhased.fasta >> "$SAMPLE".phased.haplonerate.iupac.h4.fasta
    # cat "$SAMPLE".notPhased.fasta >> "$SAMPLE".phased.haplonerate.iupac.h5.fasta
    # cat "$SAMPLE".notPhased.fasta >> "$SAMPLE".phased.haplonerate.iupac.h6.fasta
    
    for((i=1;i<=$PLOIDY;i=i+1)); do
        echo "$i"
        cat "$SAMPLE".notPhased.fasta >> "$SAMPLE".phased.haplonerate.Nmasked."$i".fasta
    done

    # cat "$SAMPLE".notPhased.fasta >> "$SAMPLE".phased.haplonerate.Nmasked.h1.fasta
    # cat "$SAMPLE".notPhased.fasta >> "$SAMPLE".phased.haplonerate.Nmasked.h2.fasta
    # cat "$SAMPLE".notPhased.fasta >> "$SAMPLE".phased.haplonerate.Nmasked.h3.fasta
    # cat "$SAMPLE".notPhased.fasta >> "$SAMPLE".phased.haplonerate.Nmasked.h4.fasta
    # cat "$SAMPLE".notPhased.fasta >> "$SAMPLE".phased.haplonerate.Nmasked.h5.fasta
    # cat "$SAMPLE".notPhased.fasta >> "$SAMPLE".phased.haplonerate.Nmasked.h6.fasta
    
    # cat "$SAMPLE".notPhased.fasta >> "$SAMPLE".phased.chimeric.h1.fasta
    # cat "$SAMPLE".notPhased.fasta >> "$SAMPLE".phased.chimeric.h2.fasta
    # cat "$SAMPLE".notPhased.fasta >> "$SAMPLE".phased.chimeric.h3.fasta
    # cat "$SAMPLE".notPhased.fasta >> "$SAMPLE".phased.chimeric.h4.fasta
    # cat "$SAMPLE".notPhased.fasta >> "$SAMPLE".phased.chimeric.h5.fasta
    # cat "$SAMPLE".notPhased.fasta >> "$SAMPLE".phased.chimeric.h6.fasta
    
done < $SAMPLEPLOIDYLIST




# # AK si pridaval _h1 ... h6 to musi ist prec

# while read SAMPLE; do
#     read PLOIDY

#     echo "processing $SAMPLE"

#     sed -i 's/_h1/supercontig/' "$SAMPLE".phased.haplonerate.iupac.h1.fasta
#     sed -i 's/_h2/supercontig/' "$SAMPLE".phased.haplonerate.iupac.h2.fasta
#     sed -i 's/supercontig_h3/supercontig/' "$SAMPLE".phased.haplonerate.iupac.h3.fasta
#     sed -i 's/supercontig_h4/supercontig/' "$SAMPLE".phased.haplonerate.iupac.h4.fasta
#     sed -i 's/supercontig_h5/supercontig/' "$SAMPLE".phased.haplonerate.iupac.h5.fasta
#     sed -i 's/supercontig_h6/supercontig/' "$SAMPLE".phased.haplonerate.iupac.h6.fasta
#     ###
#     sed -i 's/_h1/supercontig/' "$SAMPLE".phased.haplonerate.Nmasked.h1.fasta
#     sed -i 's/_h2/supercontig/' "$SAMPLE".phased.haplonerate.Nmasked.h2.fasta
#     sed -i 's/supercontig_h3/supercontig/' "$SAMPLE".phased.haplonerate.Nmasked.h3.fasta
#     sed -i 's/supercontig_h4/supercontig/' "$SAMPLE".phased.haplonerate.Nmasked.h4.fasta
#     sed -i 's/supercontig_h5/supercontig/' "$SAMPLE".phased.haplonerate.Nmasked.h5.fasta
#     sed -i 's/supercontig_h6/supercontig/' "$SAMPLE".phased.haplonerate.Nmasked.h6.fasta

# done < barbaroides_samplePloidyList_bezOutgroups.txt




# EXTRACT



read SAMPLE1 < $SAMPLEPLOIDYLIST

fastaLength=$(grep -c ">" "$SAMPLE1".phased.haplonerate.Nmasked.1.fasta)
doubleFastaLength=$((fastaLength+fastaLength))
echo $doubleFastaLength
# grep ">" "$SAMPLE1".v1.vsetky.fasta | sed s'/>//' > genes.txt

mkdir RESDIR



while read SAMPLE; do
    read PLOIDY
    
    # viem ze mam PLOIDY suborov
    
    for((pld=1;pld<=$PLOIDY;pld++)); do  # pre kazdy jeden variant = ploidia
        echo "processing $SAMPLE, pld $pld/$PLOIDY"
        
        for((line=1;line<="$doubleFastaLength";line++)); do # pre kazdy gen vzorky = berem 2 riadky naraz
            
            exonName=$(sed "${line}q;d" "$SAMPLE".phased.haplonerate.Nmasked."$pld".fasta | sed 's/>//' | sed 's/ //')
            
            # echo $(sed "${line}q;d" "$SAMPLE".v"$pld".vsetky.fasta)_h"$pld" >> ./RESDIR/"$exonName".fasta
            
            # echo ">$SAMPLE-h$pld" >> ./RESDIR/"$exonName".iupac.fasta
            echo ">$SAMPLE-h$pld" >> ./RESDIR/"$exonName".Nmasked.fasta
            # echo ">$SAMPLE-h$pld" >> ./RESDIR/"$exonName".chimeric.fasta
            
            let "line++"
            
            # echo $(sed "${line}q;d" "$SAMPLE".phased.haplonerate.iupac.h"$pld".fasta) >> ./RESDIR/"$exonName".iupac.fasta
            echo $(sed "${line}q;d" "$SAMPLE".phased.haplonerate.Nmasked."$pld".fasta) >> ./RESDIR/"$exonName".Nmasked.fasta
            # echo $(sed "${line}q;d" "$SAMPLE".phased.chimeric.h"$pld".fasta) >> ./RESDIR/"$exonName".chimeric.fasta
        done
        
    done
    
done < $SAMPLEPLOIDYLIST




# "Phasing statistics"


printf "Individual\tfazovane\tNEfazovane\n" > Phasing_statistics.txt
while read SAMPLE; do
    read PLOIDY
    
    { printf '%s\t' "$SAMPLE"
        echo "$(($(grep -c ">" "$SAMPLE".phased.haplonerate.Nmasked.1.fasta)/2))" | xargs printf '%s\t%s'
        echo "$(grep -c ">" "$SAMPLE".notPhased.fasta)" | xargs printf '%s\t\n%s'
                
    } >> Phasing_statistics.txt    
done < $SAMPLEPLOIDYLIST
