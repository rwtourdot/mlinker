#!/usr/bin/env python3

import pandas as pd
import numpy as np
import sys
import os
import matplotlib.pyplot as plt

###############################################
class variant_node:
    def __init__(self,hash,p,ref_base,base):
        self.hash = hash
        self.pos = int(p)
        self.ref_base = ref_base
        self.base = base
        if ref_base == base:
            self.variant = False
        else:
            self.variant = True
        self.links = []
        self.link_depth = []

    def add_link(self,link,num):
        self.links.append(link)
        self.link_depth.append(num)


###############################################
def load_tenx_hap_file(haplotype_file):   #,chr_choice
    fid = open(haplotype_file, 'rb')
    hap_dict = {}
    ref_base_dict = {}
    alt_base_dict = {}
    block_dict = {}
    for line in fid:
        print(line.rstrip().split("\t"))
        index,pos,ref_string,var_string,hap,ref_count,var_count,deltaE,switchE,block,span = line.rstrip().split("\t")
        rbase,vbase = ref_string.split("_")[-1],var_string.split("_")[-1]
        chrom = ref_string.split("_")[0]
        hap_dict[int(pos)] = int(hap);
        ref_base_dict[int(pos)] = rbase;
        alt_base_dict[int(pos)] = vbase;
        block_dict[int(pos)] = int(block);
    return hap_dict,ref_base_dict,alt_base_dict,block_dict


###############################################
def load_hic_phase_file(hic_link_file):
    fid = open(hic_link_file, 'rb')
    link_dict = {}
    for line in fid:
        pos1,het_string1,pos2,het_string2,count,reads = line.rstrip().split("\t")
        #print pos1,het_string1,pos2,het_string2,count
        nreads = len(reads.split(","))
        rbase1,bbase1 = het_string1.split("_")[1],het_string1.split("_")[2]
        rbase2,bbase2 = het_string2.split("_")[1],het_string2.split("_")[2]
        if het_string1 not in link_dict:
            link_dict[het_string1] = variant_node(het_string1,pos1,rbase1,bbase1)
            link_dict[het_string1].add_link(het_string2,nreads)
        else:
            link_dict[het_string1].add_link(het_string2,nreads)
        if het_string2 not in link_dict:
            link_dict[het_string2] = variant_node(het_string2,pos2,rbase2,bbase2)
            link_dict[het_string2].add_link(het_string1,nreads)
        else:
            link_dict[het_string2].add_link(het_string1,nreads)
    return link_dict


###############################################
def link_10x_blocks(link_dict,hap_dict,ref_base_dict,alt_base_dict,block_dict):
    correct_delta = []; incorrect_delta = []; block_phase = {}
    for link in link_dict:
        hash1,pos1,base1 = link_dict[link].hash,link_dict[link].pos,link_dict[link].base
        for i,link2 in enumerate(link_dict[link].links):
            depth_links = link_dict[link].link_depth[i]
            hash2,pos2,base2 = link_dict[link2].hash,link_dict[link2].pos,link_dict[link2].base
            if pos1 in hap_dict and pos2 in hap_dict:
                hap1,rbase1,vbase1,block1 = hap_dict[pos1],ref_base_dict[pos1],alt_base_dict[pos1],block_dict[pos1]
                hap2,rbase2,vbase2,block2 = hap_dict[pos2],ref_base_dict[pos2],alt_base_dict[pos2],block_dict[pos2]
                delta_p = abs(pos2-pos1)
                link_hap1 = 1 if base1 == rbase1 else -1;
                link_hap2 = 1 if base2 == rbase2 else -1;
                print("p1",pos1,block1,rbase1,vbase1,"p2",pos2,block2,rbase2,vbase2,"\t10xhap\t",hap1,hap2,"\tlink\t",base1,base2,link_hap1,link_hap2,"\t",depth_links,delta_p)
                lphase = depth_links if link_hap1*hap1 == link_hap2*hap2 else -depth_links;
                if block1 in block_phase:
                    if block2 in block_phase[block1]:
                        block_phase[block1][block2].append(lphase)
                    else:
                        block_phase[block1][block2] = [lphase]
                else:
                    block_phase[block1] = {}
                    block_phase[block1][block2] = [lphase]
                if block2 in block_phase:
                    if block1 in block_phase[block2]:
                        block_phase[block2][block1].append(lphase)
                    else:
                        block_phase[block2][block1] = [lphase]
                else:
                    block_phase[block2] = {}
                    block_phase[block2][block1] = [lphase]
    return block_phase


###############################################
###############################################
###############################################
chr_choice = "chr12"

###############################################
tenx_haplotype = "./../output/hap_solution_dec14_tenx_" + chr_choice + ".dat"
hap_dict,ref_base_dict,alt_base_dict,block_dict = load_tenx_hap_file(tenx_haplotype)  #,chr_choice

###############################################
hic_link_file = "./../output/hic_links_dec14_K562_69_" + chr_choice + ".dat"
link_dict = load_hic_phase_file(hic_link_file)
block_phase = link_10x_blocks(link_dict,hap_dict,ref_base_dict,alt_base_dict,block_dict)


###############################################
block_phase_list = sorted(list(block_phase.keys()))
block_matrix = np.zeros((len(block_phase_list),len(block_phase_list)))

for block1 in block_phase:
    for block2 in block_phase[block1]:
        index1 = block_phase_list.index(block1)
        index2 = block_phase_list.index(block2)
        #print block1,block2,index1,index2,sum(block_phase[block1][block2]),len(block_phase[block1][block2])
        if index1 == index2:
            print("intra block: ",index1,sum(block_phase[block1][block2]))
        block_matrix[index1][index2] = sum(block_phase[block1][block2])


#im = plt.imshow(block_matrix, cmap='RdBu_r')
#plt.colorbar(im, extend='both')
#plt.clim(-1, 1);
#plt.show()





















#correct_delta = []
#incorrect_delta = []
#block_phase = {}
#for link in link_dict:
#    hash1,pos1,base1 = link_dict[link].hash,link_dict[link].pos,link_dict[link].base
#    for i,link2 in enumerate(link_dict[link].links):
#        depth_links = link_dict[link].link_depth[i]
#        hash2,pos2,base2 = link_dict[link2].hash,link_dict[link2].pos,link_dict[link2].base
#        if pos1 in hap_dict and pos2 in hap_dict:
#            hap1,rbase1,vbase1,block1 = hap_dict[pos1],ref_base_dict[pos1],alt_base_dict[pos1],block_dict[pos1]
#            hap2,rbase2,vbase2,block2 = hap_dict[pos2],ref_base_dict[pos2],alt_base_dict[pos2],block_dict[pos2]
#            delta_p = abs(pos2-pos1)
#            link_hap1 = 1 if base1 == rbase1 else -1;
#            link_hap2 = 1 if base2 == rbase2 else -1;
#            print "p1",pos1,block1,rbase1,vbase1,"p2",pos2,block2,rbase2,vbase2,"\t10xhap\t",hap1,hap2,"\tlink\t",base1,base2,link_hap1,link_hap2,"\t",depth_links,delta_p
#            lphase = depth_links if link_hap1*hap1 == link_hap2*hap2 else -depth_links;
#            if block1 in block_phase:
#                if block2 in block_phase[block1]:
#                    block_phase[block1][block2].append(lphase)
#                else:
#                    block_phase[block1][block2] = [lphase]
#            else:
#                block_phase[block1] = {}
#                block_phase[block1][block2] = [lphase]
#            if block2 in block_phase:
#                if block1 in block_phase[block2]:
#                    block_phase[block2][block1].append(lphase)
#                else:
#                    block_phase[block2][block1] = [lphase]
#            else:
#                block_phase[block2] = {}
#                block_phase[block2][block1] = [lphase]

#print "correct",len(correct_delta)
#print "incorrect",len(incorrect_delta)

#print correct_delta

#plt.hist(correct_delta, bins=100)
#plt.hist(incorrect_delta, bins=100)
#plt.xlim(0,100000)
#ax = plt.gca()
#ax.set_yscale('log')
#plt.show()








        #chrom,pos1,ref_base,alt_base,filter,type,allele,ref,alt,cov,refcov,altcov,refmaj,altmaj,totref,totalt,all_im_bal = line.rstrip().split("\t")
        #if filter == "PASS":
        #    pos1_int_shift = int(pos1) - 1 # shift position
        #    if type == "SNP":
        #        #print chrom,pos1,ref_base,alt_base,filter,type,allele,ref,alt,cov,refcov,altcov,refmaj,altmaj,totref,totalt,all_im_bal
        #   al_sign = np.sign(int(hap))

                    #print chrom,pos1,ref_base,alt_base,filter,type,int(allele),al_sign,all_im_bal

#haplotype_file = "./RPE-1_Haplotype.txt"


            #if hap1 == hap2:
            #    print "same hap"


#for link in link_dict:
#    pos1 = link_dict[link].pos
#    print pos1
#    if pos1 in hap_dict:
#        print "in hap dict"

#print len(hap_dict)
#print len(link_dict)/2




#for hpos in hap_dict:
#    print hap_dict[hpos]
