#!/usr/bin/env python3

import os
import sys
import numpy as np

##########################################################
hap_solution_file = sys.argv[1] #'hap_solution_tenx_chr3_oct_20_tumor.dat'
substring_id = hap_solution_file[13:]
output_file = "block_" + substring_id

print("### output_file ### ",output_file)
##########################################################
lines = [line.split() for line in open(hap_solution_file)]
#fout = open(output_file, 'w')

##########################################################
block_cutoff = -10

energy_array = []
flip_array = []
##########################################################
counter = 0
for l in lines:
    index,position = int(l[0]),int(l[1])
    reference_allele,variant_allele = str(l[2]),str(l[3])
    hap = int(l[4])
    depth_reference,depth_variant = int(l[5]),int(l[6])
    spin_error,flip_error = float(l[7]),float(l[8])
    span = int(l[9])
    print(spin_error,flip_error)
    energy_array.append(spin_error)
    flip_array.append(flip_error)
    #fout.write(str(index)+"\t"+str(position)+"\t"+str(reference_allele)+"\t"+str(variant_allele)+"\t"+str(hap)+"\t"+str(depth_reference)+"\t"+str(depth_variant)+"\t"+str(spin_error)+"\t"+str(flip_error)+"\t"+str(span)+"\t"+str(counter)+"\n")
    if flip_error > block_cutoff:
        counter += 1


choice_bins=range(-100, 100)
hist, bin_edges = np.histogram(energy_array,bins=choice_bins)  #, density=True)
hist2, bin_edges2 = np.histogram(flip_array,bins=choice_bins)
for b in range(len(hist)):
    print(bin_edges[b],hist[b],hist2[b])
#print bin_edges
#print hist
#print energy_array


    #print index,"\t",position,"\t",reference_allele,"\t",variant_allele,"\t",hap,"\t",depth_reference,"\t",depth_variant,"\t",spin_error,"\t",flip_error,"\t",span,"\t",counter
