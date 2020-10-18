#!/usr/bin/env python3

import pandas as pd
import numpy as np
import sys
import os
import matplotlib.pyplot as plt
import copy
import pysam
import vcf


chrom_list = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX"]
#chrom_list = ["chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22"]  #,"chrX"
vcf_file_path = "./../sample_data/RPE-1.hets.chr1-X.BA.SNP.rmcent.vcf"
tenx_solution = "./../output/hap_solution_jan22_tenx_"
spin_cutoff = -200

#tenx_solution = "./../RPE1/hap_solution_jan22_tenx_"

###############################################
class pos_dictionary:
    def __init__(self,hap,rbase,vbase,block,chrom,deltaE,switchE):
        self.hap = hap
        self.rbase = rbase
        self.vbase = vbase
        self.block = block
        self.chrom = chrom
        self.spin_flipE = deltaE
        self.block_flipE = switchE

###############################################
class tenx_block:
    def __init__(self,pos,rbase,vbase):
        self.pos_array = [pos]
        self.rbase_array = [rbase]
        self.vbase_array = [vbase]
        self.n = 1

    def add_to_block(self,pos,rbase,vbase):
        self.pos_array.append(pos)
        self.rbase_array.append(rbase)
        self.vbase_array.append(vbase)
        self.n += 1

    def calc_bounds(self):
        self.min_pos = min(self.pos_array)
        self.max_pos = max(self.pos_array)
        self.mid_pos = int((min(self.pos_array)+max(self.pos_array))/2)

###############################################################
def construct_vcf_dict(vcf_file_path,chr_dict,spin_cutoff):
        vcffile = vcf.Reader(open(vcf_file_path,"r"))   #vcffile = pysam.VariantFile(vcf_file_path,index_filename=vcf_file_path+".idx")
        vcf_writer = vcf.Writer(open('./RPE-1.hets.chr1-X.BA.SNP.rmcent.sfcut.vcf', 'w'), vcffile)
        variants_list = []
        remove_counter = 0;
        keep_counter = 0;
        for rec in vcffile:
            if rec.CHROM in chr_dict:
                if (rec.POS-1) in chr_dict[rec.CHROM]:  #for tenx_rec in chr_dict[rec.CHROM]:
                    remove = False;  #print(chr_dict[rec.CHROM][rec.POS-1].spin_flipE)
                    if abs(chr_dict[rec.CHROM][rec.POS-1].spin_flipE) < abs(spin_cutoff):
                        remove = True;
                        remove_counter += 1;  #print("remove",chr_dict[rec.CHROM][tenx_rec].spin_flipE)
                    if not remove:
                        keep_counter += 1;
                        vcf_writer.write_record(rec)
        print("removed: ",remove_counter)
        print("keep: ",keep_counter)
        vcf_writer.close()
        #vcffile.close()
        return variants_list


###############################################
def load_tenx_hap_file(haplotype_file):   #,chr_choice
    fid = open(haplotype_file, 'rb')
    pos_dict = {}; block_dict = {}
    for line in fid:
        print(line.rstrip().split("\t"))
        index,pos,ref_string,var_string,hap,ref_count,var_count,deltaE,switchE,block,span = line.rstrip().split("\t")
        rbase,vbase = ref_string.split("_")[-1],var_string.split("_")[-1]
        chrom = ref_string.split("_")[0]
        pos_dict[int(pos)] = pos_dictionary(int(hap),rbase,vbase,int(block),chrom,float(deltaE),float(switchE))
        if int(block) not in block_dict:
            block_dict[int(block)] = tenx_block(int(pos),rbase,vbase)
        else:
            block_dict[int(block)].add_to_block(int(pos),rbase,vbase)
    return pos_dict,block_dict

###############################################
###############################################
###############################################

chr_dict = {}
for chr in chrom_list:
    haplotype_file = tenx_solution + chr + ".dat"
    pos_dict,block_dict = load_tenx_hap_file(haplotype_file)
    chr_dict[chr] = pos_dict


vlist = construct_vcf_dict(vcf_file_path,chr_dict,spin_cutoff)
