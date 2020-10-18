#!/usr/bin/env python3

import os
import sys
sys.path.append("/czlab/tourdot/packages/PyVCF/")
import vcf

hg38_banding = "./../hg38_info/chromosome_banding_hg38.txt"
#vcf_file_input = "./../sample_data/RPE-1.hets.chr1-X.BA.SNP.vcf"
#vcf_file_output = "./../sample_data/RPE-1.hets.chr1-X.BA.SNP.rmcent.vcf"
#vcf_file_input = "./../sample_data/BL1954_PCRFree.hets.recalibrated.vcf"
#vcf_file_output = "./../sample_data/BL1954_PCRFree.hets.recalibrated.rmcent.vcf"
vcf_file_input = "./../sample_data/HG002_GRCh38_GIAB_v.3.3.2_highconf_triophased.vcf"
vcf_file_output = "./../sample_data/HG002_GRCh38_GIAB_v.3.3.2_highconf_triophased.rmcent.vcf"

def construct_vcf_dict(vcf_file_path,vcf_file_output,hg38_centromere_dict):
        vcffile = vcf.Reader(open(vcf_file_path,"r"))   #vcffile = pysam.VariantFile(vcf_file_path,index_filename=vcf_file_path+".idx")
        vcf_writer = vcf.Writer(open(vcf_file_output,'w'), vcffile)
        variants_list = []
        for rec in vcffile:
                #print rec.CHROM,rec.POS,rec.REF,rec.ALT
                #print rec.CHROM,hg38_centromere_dict[rec.CHROM]
                remove = False;
                nregions1 = len(hg38_centromere_dict[rec.CHROM]["start"])
                for j in range(nregions1):
                        start_loop,end_loop = hg38_centromere_dict[rec.CHROM]["start"][j],hg38_centromere_dict[rec.CHROM]["end"][j]
                        if int(rec.POS) >= start_loop and int(rec.POS) <= end_loop:
                                remove = True;
                if not remove:
                        print(rec.CHROM,rec.POS,rec.REF,rec.ALT)
                        print("to keep")
                        vcf_writer.write_record(rec)
        return variants_list


##############################################################################
def load_centromere_positions(hg38_centromere_file):
        lfile = open(hg38_centromere_file, 'r')
        lines = lfile.readlines()
        hg38_centromere_dict = {}
        for l in lines:
                chromosome,start,end,arm,gram = l[:-2].split("\t")
                if gram == "gva" or gram == "stal" or gram == "ace":
                        #print(chromosome,start,end)
                        if chromosome not in hg38_centromere_dict:
                                hg38_centromere_dict[chromosome] = {}
                                hg38_centromere_dict[chromosome]["start"] = [int(start)]
                                hg38_centromere_dict[chromosome]["end"] = [int(end)]
                        else:
                                hg38_centromere_dict[chromosome]["start"].append(int(start))
                                hg38_centromere_dict[chromosome]["end"].append(int(end))
        return hg38_centromere_dict

##############################################################################
hg38_centromere_dict = load_centromere_positions(hg38_banding)
vlist = construct_vcf_dict(vcf_file_input,vcf_file_output,hg38_centromere_dict)



