#!/usr/bin/env python3

import pysam
import sys
import gzip
#import re
#import os
import vcf


###########################################################################
chrom_list = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX"]
start = 0
end = 300000000
#vcf_file_path="/czlab/tourdot/BL1954_variants/BL1954_PCRFree.hets.recalibrated.copy.vcf.gz"
#vcf_file_path="/czlab/Data/HCC1954_HaplotypeCalling/v.2/HCC1954_GT/HCC1954_illumina_hets.BA.GT.vcf.gz"
vcf_file_path="/czlab/Data/HCC1954_HaplotypeCalling/v.2/HCC1954_GT/HCC1954_PCRFree_hets.BA.GT.vcf.gz"
#cov_output_file = "./../output/het_coverage_june11_illumina_"
cov_output_file = "./../output/het_coverage_june11_PCRFree_"

###########################################################################
#sname = "HCC1954"
#sname = "HCC1954_BL"
#sname = "HCC1954BL"
sname = "BL"

sname2 = sname
#sname2 = "BL"

###########################################################################
def construct_cov_file(vcf_file_path,chr_choice,start,end,sample_name,cfile):
	f = open(cfile,'w')
	base_order = ['I',"D","G","C","A","T"]
	vcffile = vcf.Reader(open(vcf_file_path,"r"))
	for rec in vcffile.fetch(chr_choice,start,end):
		#record = rec.genotype(sample_name)
		record = rec.samples[1]
		#record = rec.samples[0]
		if "AD" in rec.FORMAT.split(":") and record.data[1] != None:
			if len(rec.ALT[0]) == 1 and len(rec.REF) == 1:				#and record.is_het and record.is_variant:
				base_cov = {}; base_cov["I"] = 0; base_cov["D"] = 0; base_cov["G"] = 0; base_cov["C"] = 0; base_cov["A"] = 0; base_cov["T"] = 0;
				ref_base = str(rec.REF); var_base = str(rec.ALT[0]); position = str(rec.POS-1)
				base_cov[ref_base] = record.data[1][0]; base_cov[var_base] = record.data[1][1];
				output_str = position+"_"+ref_base+"_"+var_base+"\t"+position+"\t"+ref_base+":"+var_base+"\t"
				sum_bases = 0;
				for item in base_order:
					output_str = output_str + str(item) + "|" + str(base_cov[item]) + "\t"
					sum_bases += base_cov[item];
				output_str = output_str + str(sum_bases)
				#print cfile,chr_choice,output_str  # chr_choice + "_" +
				print(cfile,chr_choice,output_str)  # chr_choice + "_" +
				f.write(output_str+"\n")
	f.close()


###########################################################################
for chr_choice in chrom_list:
	#cfile = cov_output_file + sname + "_" + chr_choice + ".dat"
	cfile = cov_output_file + sname2 + "_" + chr_choice + ".dat"
	construct_cov_file(vcf_file_path,chr_choice,start,end,sname,cfile)


###########################################################################



                                #print(rec.samples)
                                #print(rec.samples[0])
                                #print(rec.samples[1])
                                #print(rec.samples[1].data[1])

