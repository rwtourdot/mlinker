#!/bin/bash

#CHR_list=('chr1' 'chr2' 'chr3' 'chr4' 'chr5' 'chr6' 'chr7' 'chr8' 'chr9' 'chr10' 'chr11' 'chr12' 'chr13' 'chr14' 'chr15' 'chr16' 'chr17' 'chr18' 'chr19' 'chr20' 'chr21' 'chr22' 'chrX')
CHR_list=('chr22' 'chr21' 'chr20' 'chr19' 'chr18' 'chr17' 'chr16' 'chr15' 'chr14' 'chr13' 'chr12' 'chr11' 'chr10' 'chr9' 'chr8' 'chr7' 'chr6' 'chr5' 'chr4' 'chr3' 'chr2' 'chr1' 'chrX')

#input_vcf_file="./vcf_data/RPE-1.hets.chr1-X.BA.SNP.rmcent.vcf"
#input_vcf_file="./vcf_data/RPE-1.hets.chr1-X.BA.SNP.sites.rmcent.vcf"
input_vcf_file="./vcf_data/RPE-1.hets.chr1-X.BA.SNP.sites.snp_filter.vcf"
input_bam_file="/czlab/Data/RPE-1_10X/RPE-1_bam/phased_possorted_bam.bam"
tech="tenx"

string_id="may9_RPE1"
nohup_file="nohup_extract_"$string_id"_"$tech"_"


for chr_choice in "${CHR_list[@]}"
do
        loop_nohup_file=$nohup_file$chr_choice".out"
        echo $loop_nohup_file
        echo $string_id
        nohup time ./mlinker extract -c $chr_choice -i $input_bam_file -v $input_vcf_file -e $tech -n $string_id > $loop_nohup_file &
done



