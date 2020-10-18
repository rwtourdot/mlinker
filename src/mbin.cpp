#include "mbin.h"

void mbin::initialize_mbin( int bin_n ) {
    bin_index = bin_n;
    num_mbins = 0;
}

void mbin::set_energies( int bin_ind, double switch_prob ) {
    merged_bins.push_back(bin_ind);
    cn_switch_prob.push_back(switch_prob);
    num_mbins += 1;
}


