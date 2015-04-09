#!/bin/bash

lib=$1

## run only once!
#cat ../out_bins/$lib/orphan_1.fasta  | blastall -p blastn -d ../db/$lib\_all_amplicons -e 0.001 -m 8 -v 1 -b 2 -F F > $lib\_1_blast-compact.out
#cat ../out_bins/$lib/orphan_2.fasta | blastall -p blastn -d ../db/$lib\_all_amplicons -e 0.001 -m 8 -v 1 -b 2 -F F > $lib\_2_blast-compact.out


echo "*** Library $lib ***"
echo "Total orphan reads:"
total_reads=`cat $lib\_1_blast-compact.out | cut -f1 | uniq | wc -l` 
echo $total_reads


echo "Reads with BLAST-mm > 80:"
reads_with_long_enough_mm=`cat $lib\_1_blast-compact.out | awk '$4>80' | cut -f1 | uniq | wc -l`

echo "Potentially mistakenly orphan reads:"
echo $reads_with_long_enough_mm

echo "Calibration reads among the potentially mistakenly orphan reads:"
echo `cat $lib\_1_blast-compact.out | awk '$4>80' | grep 'D_' | cut -f1 | uniq | wc -l`
echo "";

#../assign_blast_output_into_bins.pl ../../All_NGS_library_information_A.txt A_1_blast-compact.out A_2_blast-compact.out
