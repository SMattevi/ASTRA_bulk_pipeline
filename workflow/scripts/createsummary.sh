bam=$1
bamQC=$2
threads=$3
tech=$4

## compute flagstats and indxstats for both filtered and not filtered mapping results
samtools flagstat -@ $threads $bam > results/$tech/$tech.flagstat.txt
samtools idxstats -@ $threads $bam > results/$tech/$tech.idxstat.txt

samtools flagstat -@ $threads $bamQC > results/$tech/$tech.flagstat.MAPQ30.txt  
samtools idxstats -@ $threads $bamQC > results/$tech/$tech.idxstat.MAPQ30.txt

## create summary file
samtools view -@ $threads -q 5 -f 0x2 $bam > tmp_sam_file.sam
total_uniq_mapped=$( wc -l tmp_sam_file.sam | cut -d ' ' -f1 )  ## number of unique mapped reads

total_uniq_mapped=$((${total_uniq_mapped}/2))

total_pairs=$(grep 'paired in' results/$tech/$tech.flagstat.txt | cut -d ' ' -f1) 
total_pairs=$((${total_pairs} / 2)) 
total_pairs_mapped=$(grep 'with itself and mate mapped'  results/$tech/$tech.flagstat.txt | cut -d ' ' -f1)
total_pairs_mapped=$((${total_pairs_mapped} / 2)) 
total_mito_mapped=$(grep MT results/$tech/$tech.idxstat.txt | cut -f3)
total_mito_unmapped=$(grep MT results/$tech/$tech.idxstat.txt | cut -f4)
total_mito=$((${total_mito_mapped} / 2 + ${total_mito_unmapped} / 2))
total_mito_mapped=$((${total_mito_mapped} / 2))
total_dups=$(grep 'primary duplicates' results/$tech/$tech.flagstat.txt | cut -d ' ' -f1)
total_dups=$((${total_dups} / 2))



total_pairs_MAPQH=$(grep 'with itself and mate mapped'  results/$tech/$tech.flagstat.MAPQ30.txt | cut -d ' ' -f1)
total_pairs_MAPQH=$((${total_pairs_MAPQH}/2)) 
total_mito_MAPQH=$(grep MT results/$tech/$tech.idxstat.MAPQ30.txt | cut -f3)
total_mito_MAPQH=$((${total_mito_MAPQH}/2))
total_dups_MAPQH=$(grep 'primary duplicates' results/$tech/$tech.flagstat.MAPQ30.txt | cut -d ' ' -f1)
total_dups_MAPQH=$((${total_dups_MAPQH}/2))

rm results/$tech/$tech.idxstat.txt 
rm results/$tech/$tech.flagstat.txt 

rm results/$tech/$tech.flagstat.MAPQ30.txt 
rm results/$tech/$tech.idxstat.MAPQ30.txt 

#print to file
echo "Total_Pairs    $total_pairs" > results/$tech/alignment/$tech.MappingStats 
echo "Total_Pairs_Mapped    $total_pairs_mapped" >> results/$tech/alignment/$tech.MappingStats 
echo "Total_Uniq_Mapped    $total_uniq_mapped" >> results/$tech/alignment/$tech.MappingStats 
#echo "Total_Mito    $total_mito" >> results/$tech/alignment/$tech.MappingStats 
echo "Total_Mito_Mapped    $total_mito_mapped" >> results/$tech/alignment/$tech.MappingStats 
echo "Total_Dups    $total_dups" >> results/$tech/alignment/$tech.MappingStats 


echo "Total_Pairs_MAPQ30    $total_pairs_MAPQH" >> results/$tech/alignment/$tech.MappingStats 
echo "Total_Mito_MAPQ30    $total_mito_MAPQH" >> results/$tech/alignment/$tech.MappingStats 
echo "Total_Dups_MAPQ30    $total_dups_MAPQH" >> results/$tech/alignment/$tech.MappingStats 

rm tmp_sam_file.sam
