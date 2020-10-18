#!/bin/bash

#CHR_list=('chr1' 'chr2' 'chr3' 'chr4' 'chr5' 'chr6' 'chr7' 'chr8' 'chr9' 'chr10' 'chr11' 'chr12' 'chr13' 'chr14' 'chr15' 'chr16' 'chr17' 'chr18' 'chr19' 'chr20' 'chr21' 'chr22' 'chrX')
CHR_list=('chr22' 'chr21' 'chr20' 'chr19' 'chr18' 'chr17' 'chr16' 'chr15' 'chr14' 'chr13' 'chr12' 'chr11' 'chr10' 'chr9' 'chr8' 'chr7' 'chr6' 'chr5' 'chr4' 'chr3' 'chr2' 'chr1' 'chrX')

#input_graph_file="./output/graph_variant_oct18_HCC1954_all_"
input_graph_file="./output/graph_variant_may9_RPE1_tenx_"
input_scaffold_file="./output/hap_full_scaffold_may10_RPE1_22_"
#string_id="oct18_HCC1954_all"
string_id="may10_RPE1_22"
nohup_file="nohup_recover_"$string_id"_"

technology="tenx"

for chr_choice in "${CHR_list[@]}"
do
        loop_nohup_file=$nohup_file$chr_choice".out"
	loop_scaffold_file=$input_scaffold_file$chr_choice".dat"
	loop_graph_file=$input_graph_file$chr_choice".dat"
        echo $loop_nohup_file
	echo $loop_graph_file
	echo $loop_scaffold_file
        echo $string_id
	#SECONDS=0
        nohup time ./mlinker recover -c $chr_choice -g $loop_graph_file -i $loop_scaffold_file -n $string_id -e $technology > $loop_nohup_file &
        #elapsedseconds=$SECONDS
        #echo "elapsed_seconds: "$elapsedseconds >> $loop_nohup_file
done



