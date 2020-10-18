#include "bin_reference.h"

void gen_bin::find_unique_bx() {
    num_bx = connected_bx.size();
    std::set<std::string> set_temp(connected_bx.begin(),connected_bx.end());
    unique_bx.resize(set_temp.size());
    std::copy(set_temp.begin(),set_temp.end(),unique_bx.begin());
    num_unique_bx = unique_bx.size();
}

void gen_bin::add_connection( std::string chr, int binname ) { 
    chr_connections[chr][binname] += 1; 
};

void gen_bin::add_connected_read( std::string bxname ) { connected_bx.push_back(bxname); }

void gen_bin::set_values( int position ) { pos = position; }

void cn_bin::initialize( int pos ) { 
    start_pos = pos;
    num_hets = 0;
    switch_spin = 1; 
}

void cn_bin::flip_hap() { switch_spin = switch_spin*-1; }

void cn_bin::add_het( int hpos, int hap, int block, std::string rbase, std::string vbase, int rcov, int vcov ) {
    het_pos.push_back(hpos);
    het_hap.push_back(hap);
    het_block.push_back(block);
    ref_base.push_back(rbase); var_base.push_back(vbase);
    ref_cov.push_back(rcov);   var_cov.push_back(vcov);
    num_hets += 1; 
}

void cn_bin::cov_calc() {
    std::set<int> set_block(het_block.begin(),het_block.end());
    num_blocks = set_block.size();
    int i = 0;
    hapA_tot = 0;
    hapB_tot = 0;
    for (auto k : het_hap) {
	cout << k << "  " << ref_cov[i] << "  " << var_cov[i] << endl;
	if (k == 1) { hapA_tot += ref_cov[i]; hapB_tot += var_cov[i]; }
	else if (k == -1) { hapA_tot += var_cov[i]; hapB_tot += ref_cov[i]; }
	i++;	
    }
    if (num_hets > 0 && hapA_tot > 0) { hapA_avg = (double)hapA_tot/(double)num_hets; }
    if (num_hets > 0 && hapB_tot > 0) { hapB_avg = (double)hapB_tot/(double)num_hets; }
    if (num_hets == 0 || hapA_tot == 0) { hapA_avg = 0.0; }
    if (num_hets == 0 || hapB_tot == 0) { hapB_avg = 0.0; }
}

void cn_bin::enter_coverage( double hap1_cn, double hap2_cn, int block ) {
     hapA_avg = hap1_cn;
     hapB_avg = hap2_cn;
     num_hets = 10;
     num_blocks = 1;
     het_block.push_back(block);
}




