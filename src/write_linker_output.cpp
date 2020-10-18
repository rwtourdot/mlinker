#include "write_linker_output.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////
void write_het_coverage( std::unordered_map<std::string,variant_node>& var_dict, std::string coverageFile, std::string chr_choice ) {
        ofstream ofile; ofile.open(coverageFile);
	vector<char> base_order = {'I','D','G','C','A','T'};
        for (auto& it : var_dict) {
		if ( it.second.var ) {
                ofile << it.first << "\t" << it.second.pos << "\t";
                ofile << it.second.ref_base << ":" << it.second.var_base << "\t";
		for(auto base : base_order) {
			//std::cout << base << "|" << it.second.base_dict[base] << "\n";
			ofile << base << "|" << it.second.base_dict[base] << "\t";
		}
                //for (auto& it2 : it.second.base_dict) { ofile << it2.first << "|" << it2.second << "\t"; }
		ofile << it.second.total_bases << "\t";
                ofile << endl;
		}
        }
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
void write_het_bx_coverage( std::unordered_map<std::string,variant_node>& var_dict, std::string coverageFile, std::string chr_choice ) {
	// print each bx tag
        ofstream ofile; ofile.open(coverageFile);
	vector<char> base_order = {'I','D','G','C','A','T'};
        for (auto& it : var_dict) {
                if ( it.second.var ) {
                char rbase = it.second.ref_base[0]; char vbase = it.second.var_base[0];
                std::set<std::string> ref_bx_set = it.second.base_dict_set[rbase];
                std::set<std::string> var_bx_set = it.second.base_dict_set[vbase];
                ofile << it.first << "\t" << it.second.pos << "\t";
                ofile << it.second.ref_base << ":" << it.second.var_base << "\t";
                for(auto base : base_order) {
                        //std::cout << base << "|" << it.second.base_dict[base] << "\n";
                        ofile << base << "|" << it.second.base_dict[base] << "\t";
                }
                //for (auto& it2 : it.second.unique_dict) { ofile << it2.first << "|" << it2.second << "\t"; }
                ofile << it.second.unique_total << "\t";
		ofile << it.second.ref_base << ":";
		for ( std::set<std::string>::iterator it4=ref_bx_set.begin(); it4!=ref_bx_set.end(); ++it4 ) {
			ofile << *it4;
			if (std::distance(ref_bx_set.begin(), it4) < (ref_bx_set.size()-1) ) { ofile << ","; }
		}
		ofile << "\t";
		ofile << it.second.var_base << ":";
                for ( std::set<std::string>::iterator it5=var_bx_set.begin(); it5!=var_bx_set.end(); ++it5 ) {
			ofile << *it5;
			if (std::distance(var_bx_set.begin(), it5) < (var_bx_set.size()-1) ) { ofile << ","; }
		}
                ofile << endl;
                }
        }
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
void write_het_coverage_filter( std::unordered_map<std::string,variant_node>& var_dict, std::string coverageFile ) {
        ofstream ofile; ofile.open(coverageFile);
	vector<char> base_order = {'I','D','G','C','A','T'};
        for (auto& it : var_dict) {
		//cout << it.first << "  " << it.second.filter << "  " << it.second.var << endl;
                if ( it.second.var ) {
		if ( it.second.filter ) {
        		ofile << it.first << "\t" << it.second.pos << "\t";
        		ofile << it.second.ref_base << ":" << it.second.var_base << "\t";
                	for(auto base : base_order) {
                        	//std::cout << base << "|" << it.second.base_dict[base] << "\n";
                        	ofile << base << "|" << it.second.base_dict[base] << "\t";
                	}
        		//for (auto& it2 : it.second.base_dict) { ofile << it2.first << "|" << it2.second << "\t"; }
        		ofile << it.second.total_bases << "\t";
        		ofile << endl;
		}
                }
        }
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
void write_hic_links( std::unordered_map<std::string,variant_node>& var_dict, std::string outputFile, std::string chr_choice ) {
        ofstream ofile; ofile.open(outputFile);
        for (auto& it : var_dict) {
                int pos1 = it.second.pos;
                std::string base1 = it.first;
                if (it.second.unique_total > 0) {
                for (auto it2 : it.second.connections) {  //int pos2 = 100;  //var_dict[it2.first].pos;
                        int pos2 = var_dict[it2.first].pos;
			int nconnections = it2.second;
                        std::string base2 = it2.first;
			std::vector<std::string> v1 = var_dict[it2.first].connected_reads_long_form;
			std::vector<std::string> v2 = var_dict[it.first].connected_reads_long_form;
                        ofile << pos1 << "\t" << base1 << "\t" << pos2 << "\t" << base2 << "\t" << nconnections << "\t";
			bool first = true; ofile << "[";
			for (int j = 0; j < v1.size(); j++) {
				if ( std::find(v2.begin(), v2.end(), v1[j]) != v2.end() ) {
					if (!first) { ofile << ','; }
					first = false;
					ofile << v1[j];
				}
			}
			ofile << "]" << endl;
                }
                }
        }
        ofile.close();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
void write_scaffold( std::string chr_choice, std::string outputFile, block_dictionary& bl_dict, coord_dictionary& pdict, int centromere_pos ) {
        ofstream ofile; ofile.open(outputFile);
	for (int i = 0; i < pdict.num_paired; i++) {
		int pos1 = pdict.double_positions[i];
                int hap1 = pdict.haplotype[i];
                int block1 = pdict.block[i];
		if ( std::find( bl_dict.subset_blocks.begin(), bl_dict.subset_blocks.end(), block1 ) != bl_dict.subset_blocks.end() )  {
                double deltaE = pdict.deltaE[i];
                double switchE = pdict.switchE[i];
                std::string ref_hash = pdict.ref_handle[i];
                std::string alt_hash = pdict.alt_handle[i];
		int block_hap1 = bl_dict.subset_haplotype[bl_dict.block_map[block1]];
		std::string arm1 = "p"; if ( pos1 > centromere_pos ) { arm1 = "q"; };
		ofile << chr_choice << "\t" << arm1 << "\t" << pos1 << "\t" << ref_hash << "\t" << alt_hash << "\t" << hap1*block_hap1 << "\t" << "block:" << "\t" << block_hap1 << "\t" << block1 << "\t" << "tenx_energy(spin/block):" << "\t" << deltaE << "\t" << switchE << endl;
		}
	}
  ofile.close();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
void write_recovered( std::string outputFile, std::map<int,recovered_node>& recovered_map ) {
  ofstream ofile; ofile.open(outputFile);
  ofile << "#position" << "\t" << "scaffold_haplotype" << "\t" << "ref_numA" << "\t" << "ref_numB" << "\t" << "alt_numA" << "\t" << "alt_numB" << "\t" << "ref_numA_uniq_pos" << "\t" << "ref_numB_uniq_pos" << "\t" << "alt_numA_uniq_pos" << "\t" << "alt_numB_uniq_pos" << "\t" << "ref_numA_uniq_hash" << "\t" << "ref_numB_uniq_hash" << "\t" << "alt_numA_uniq_hash" << "\t" << "alt_numB_uniq_hash"<< "\t" << endl;
	for (auto& it : recovered_map) {
    recovered_node loop_node = it.second;  //cout << loop_node.pos << "\t" << loop_node.hap << endl;
    ofile << loop_node.pos << "\t" << loop_node.hap << "\t" << loop_node.ref_numA << "\t" << loop_node.ref_numB << "\t" << loop_node.alt_numA << "\t" << loop_node.alt_numB << "\t" << loop_node.upos_ref_numA << "\t" << loop_node.upos_ref_numB << "\t" << loop_node.upos_alt_numA << "\t" << loop_node.upos_alt_numB << "\t" << loop_node.ref_numA_hashcall << "\t" << loop_node.ref_numB_hashcall << "\t" << loop_node.alt_numA_hashcall << "\t" << loop_node.alt_numB_hashcall << endl;
	}
  ofile.close();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
void write_scaffold_energy( std::string chr_choice, std::string outputFile, block_dictionary& bl_dict, int centromere_pos ) {
        ofstream ofile; ofile.open(outputFile);
        for (int j = 0; j < bl_dict.length_subset; j++) {
                int block = bl_dict.subset_blocks[j];
                int nhets = bl_dict.num_hets[block];
                int blhap = bl_dict.subset_haplotype[j];
                double deltaE = bl_dict.subset_deltaE[j];
                double switchE = bl_dict.subset_switchE[j];
                int min_pos = *min_element(bl_dict.bl_pos[block].begin(), bl_dict.bl_pos[block].end());
                int max_pos = *max_element(bl_dict.bl_pos[block].begin(), bl_dict.bl_pos[block].end());
		std::string arm1 = "p"; if ( min_pos > centromere_pos ) { arm1 = "q"; };
                //cout << j << "\t" << block << "\t" << min_pos << "\t" << max_pos << "\t" << nhets << "\t" << blhap << "\t" << deltaE << "\t" << switchE << endl;
                ofile << j << "\t" << block << "\t" << arm1 << "\t" << min_pos << "\t" << max_pos << "\t" << nhets << "\t" << blhap << "\t" << deltaE << "\t" << switchE << endl;
        }
        ofile.close();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
void write_link_network( std::unordered_map<std::string,variant_node>& var_dict, std::string outputFile, std::string chr_choice ) {   // cout << " writing network file " << endl;
        ofstream ofile; ofile.open(outputFile);
        for (auto& it : var_dict) {
                ofile << chr_choice << "_" << it.first << "  " << it.second.pos << "  ["; //<< endl;
                bool first = true;
                for (auto it2 : it.second.connections) { if (!first) { ofile << ','; } first = false; ofile << chr_choice << "_" << it2.first; }
                ofile << "] [";
                bool first_num = true;
                for (auto it2 : it.second.connections) { if (!first_num) { ofile << ','; } first_num = false; ofile << it2.second; }
                ofile << "] {"; std::string base_output;
                for (auto it2 : it.second.base_dict) {
                        if(!base_output.empty()) { base_output += ","; }
                        string ssbase(1,it2.first);
                        base_output += ssbase + ":" + std::to_string(it2.second);
                }
                ofile << base_output << "} paired " << std::boolalpha << it.second.paired << endl;
        }
        ofile.close();
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
void write_variant_link_list( std::unordered_map<std::string,variant_node>& var_dict, read_graph& rgraph, std::string outputFile, std::string chr_choice ) {
        ofstream ofile; ofile.open(outputFile);
        for (auto& it : var_dict) {
                std::string variant_string = it.first;
                std::vector<std::string> hash_connections = it.second.connected_reads_long_form;
                std::vector<std::string> hash_readnames = it.second.connected_readnames;
                for (int j = 0; j < hash_connections.size(); j++) {
                        ofile << variant_string << "\t" << chr_choice << "\t" << it.second.pos << "\t" << hash_connections[j] << "\t" << hash_readnames[j] << endl;
                }
        }
        ofile.close();
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
void write_hash_link_list( std::unordered_map<std::string,variant_node>& var_dict, read_graph& rgraph, std::string outputFile, std::string chr_choice ) {
        ofstream ofile; ofile.open(outputFile);
        for (auto& it : rgraph) {
        	std::vector<std::string> het_connections = rgraph[it.first].connected_strings;
        	std::vector<std::string> het_readnames = rgraph[it.first].connected_readnames;
        	for (int j = 0; j < het_connections.size(); j++) {
                	ofile << it.first << "\t" << chr_choice << "\t" << het_connections[j] << "\t" << het_readnames[j] << endl;
                }
        }
        ofile.close();
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
//void write_het_links( std::unordered_map<std::string,variant_node>& var_dict, std::string outputFile, std::string chr_choice ) {
//        ofstream ofile; ofile.open(outputFile);
//        for (auto& it : var_dict) {
//		std::string variant1 = chr_choice + "_" + it.first;
//		for (auto it2 : it.second.connections) {
//			std::string variant2 = chr_choice + "_" it2.first;
//		}
//	}
//	ofile.close();
//}

////////////////////////////////////////////////////////////////////////////////////////////////////////
void write_bx_links( read_graph &rgraph, std::string outputFile ) {
        ofstream ofile; ofile.open(outputFile);
        for (auto& it : rgraph) {
		ofile << it.first << "\t" << it.second.num_cnx << "\t";
		std::vector<std::string> het_connections = rgraph[it.first].connected_strings;
		bool first_num = true;
		for (int j = 0; j < het_connections.size(); j++) {
			if (!first_num) { ofile << ','; }
			first_num = false;
			ofile << het_connections[j];
		}
		ofile << endl;
	}
	ofile.close();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
void write_hap_solution_loh( std::unordered_map<std::string,variant_node>& var_dict, std::string hapsolutionFile, coord_dictionary& pdict, std::string chr_choice ) {
        ofstream ofile; ofile.open(hapsolutionFile);
        std::unordered_map<int,int> opposite; opposite[0] = 1; opposite[1] = 0;
        for (int i = 0; i < pdict.num_paired; i++) {
                int p = pdict.sorted_paired_positions[i];
                int m = pdict.ref_index[p];
                int pair = pdict.paired[p];
                std::string ref_hash,alt_hash;
                int base_number_ref,base_number_alt;
                if (pair == 1) {
                        ref_hash = chr_choice + "_" + pdict.paired_dict[p][opposite[m]];
                        alt_hash = chr_choice + "_" + pdict.paired_dict[p][m];
                        base_number_ref = pdict.base_number[p][opposite[m]];
                        base_number_alt = pdict.base_number[p][m];
                } else {
                        int var = pdict.ref_index[p];
                        if (var == 1) {
                                ref_hash = chr_choice + "_" + std::to_string(p) + "_l_h";
                                alt_hash = chr_choice + "_" + pdict.unpaired_dict[p]; //pdict.paired_dict[p][m];
                                base_number_ref = 0;
                                base_number_alt = pdict.base_number[p][0];
                        } else {
                                ref_hash = chr_choice + "_" + pdict.unpaired_dict[p]; //pdict.paired_dict[p][opposite[m]];
                                alt_hash = chr_choice + "_" + std::to_string(p) + "_l_h";
                                base_number_ref = pdict.base_number[p][0];
                                base_number_alt = 0;
                        }
                }
                //cout << i << "\t" << p << "\t" << ref_hash << "\t" << alt_hash << "\t" << pdict.haplotype[i] << "\t" << base_number_ref << "\t" << base_number_alt << "\t" << pdict.deltaE[i] << "\t" << pdict.switchE[i] << "\t" << pdict.block[i] << "\t" << pdict.span_bound[i] << endl;
                ofile << i << "\t" << p << "\t" << ref_hash << "\t" << alt_hash << "\t" << pdict.haplotype[i] << "\t" << base_number_ref << "\t" << base_number_alt << "\t" << pdict.deltaE[i] << "\t" << pdict.switchE[i] << "\t" << pdict.block[i] << "\t" << pdict.span_bound[i] << endl;
        }
        ofile.close();
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
void write_hap_solution( std::unordered_map<std::string,variant_node>& var_dict, std::string hapsolutionFile, coord_dictionary& pdict, std::string chr_choice ) {
        ofstream ofile; ofile.open(hapsolutionFile);
        std::unordered_map<int,int> opposite; opposite[0] = 1; opposite[1] = 0;
        for (int i = 0; i < pdict.num_paired; i++) {
                int p = pdict.sorted_paired_positions[i];
                int m = pdict.ref_index[p];
                ofile << i << "\t" << p << "\t" << chr_choice << "_" << pdict.paired_dict[p][opposite[m]] << "\t" << chr_choice << "_" << pdict.paired_dict[p][m] << "\t" << pdict.haplotype[i] << "\t" << pdict.base_number[p][opposite[m]] << "\t" << pdict.base_number[p][m] << "\t" << pdict.deltaE[i] << "\t" << pdict.switchE[i] << "\t" << pdict.block[i] << "\t" << pdict.span_bound[i] << endl;
        }
        ofile.close();
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
void write_bin_matrix( std::unordered_map<int,gen_bin>& bx_map, std::string outputFile, std::string contig_name, int binsize ) {
        ofstream ofile; ofile.open(outputFile);
        for (auto& it : bx_map) {
                for (auto& it2 : it.second.connections) {
			ofile << contig_name << "\t" << it.first*binsize << "\t" << "-" << "\t" << contig_name << "\t" << it2.first*binsize << "\t" << "+" << "\t" << it2.second << endl;
		}
        }
        ofile.close();
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
void write_bin_bx_cov( std::unordered_map<int,gen_bin>& bx_map, std::string outputFile, std::string contig_name, int binsize ) {
        ofstream ofile; ofile.open(outputFile);
        for (auto& it : bx_map) {
        	ofile << contig_name << "\t" << it.first*binsize << "\t" << it.second.num_bx << "\t" << it.second.num_unique_bx << endl;
        }
        ofile.close();
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
void write_bin_matrix_genome( full_map& chr_map, std::string outputFile, int binsize ) {
        ofstream ofile; ofile.open(outputFile);
        for (auto& it : chr_map) {
                for (auto& it2 : it.second) {
                	for (auto& it3 : it2.second.chr_connections) {
				for (auto& it4 : it3.second) {
                        		ofile << it.first << "\t" << it2.first*binsize << "\t" << "-" << "\t" << it3.first << "\t" << it4.first*binsize << "\t" << "+" << "\t" << it4.second << endl;
				}
                	}
        	}
	}
        ofile.close();
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
void write_block_phasing_matrix( block_dictionary& bl_dict, std::string outputFile ) {
        ofstream ofile; ofile.open(outputFile);
	for (int i = 0; i < bl_dict.length_subset; i++) {
		for (int j = 0; j < bl_dict.length_subset; j++) {
			if (bl_dict.subset_matrix(i,j) != 0) {
				ofile << i << "\t" << j << "\t" << bl_dict.block_map_inverted[i] << "\t" << bl_dict.block_map_inverted[j] <<  "\t" << bl_dict.subset_haplotype[i] << "\t" << bl_dict.subset_haplotype[j] << "\t" << bl_dict.subset_matrix(i,j) << endl;
			}
		}
	}
        ofile.close();
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
void write_phased_sv( std::vector<sv_entry>& sv_list, std::string outputFile ) {
	ofstream ofile; ofile.open(outputFile);
        for (int i = 0; i < sv_list.size(); i++) {
                ofile << sv_list[i].chr_one << ":" << sv_list[i].pos_one << "  " << sv_list[i].chr_two << ":" << sv_list[i].pos_two << endl;
                for (auto& it : sv_list[i].hap_end1) { ofile << "0 " << sv_list[i].chr_one << " " << it.first << "  " << it.second << "  " << sv_list[i].num_end1[it.first] << endl; }
                for (auto& it : sv_list[i].hap_end2) { ofile << "1 " << sv_list[i].chr_two << " " << it.first << "  " << it.second << "  " << sv_list[i].num_end2[it.first] << endl; }
                ofile << "######" << endl;
        }
        ofile.close();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
void write_cn_phased( cn_map& chromosome_map, std::string outputFile, std::vector<int>& good_bins, std::vector<int>& merged_bins ) {
	ofstream ofile; ofile.open(outputFile);
	cout << "write output " << endl;
	int k = 0;
	int p = 0;
	for (auto ind : good_bins) {
		for (int j = 0; j < chromosome_map[ind].num_hets; j++) {
			int hap = chromosome_map[ind].switch_spin*chromosome_map[ind].het_hap[j];
			int orig_hap = chromosome_map[ind].het_hap[j];
			int hapA_cov = 0;
			int hapB_cov = 0;
			int orig_hapA_cov = 0;
			int orig_hapB_cov = 0;
			if (hap == 1) {  hapA_cov = chromosome_map[ind].ref_cov[j]; hapB_cov = chromosome_map[ind].var_cov[j]; }
			if (hap == -1) { hapB_cov = chromosome_map[ind].ref_cov[j]; hapA_cov = chromosome_map[ind].var_cov[j]; }
                        if (orig_hap == 1) {  orig_hapA_cov = chromosome_map[ind].ref_cov[j]; orig_hapB_cov = chromosome_map[ind].var_cov[j]; }
                        if (orig_hap == -1) { orig_hapB_cov = chromosome_map[ind].ref_cov[j]; orig_hapA_cov = chromosome_map[ind].var_cov[j]; }
			ofile << k << "\t" << ind << "\t" << chromosome_map[ind].het_pos[j] << "\t" << chromosome_map[ind].ref_base[j] << "\t" << chromosome_map[ind].var_base[j] << "\t" << hap << "\t" << hapA_cov <<  "\t" << hapB_cov << "\t" << orig_hap << "\t" << orig_hapA_cov <<  "\t" << orig_hapB_cov << "\t" << chromosome_map[ind].het_block[j] << "\t" << chromosome_map[ind].switch_spin << "\t" << merged_bins[p] << endl;
			k += 1;
		}
		p+=1;
	}
	cout << "number of bins: " << p << endl;
	cout << "number of hets: " << k << endl;
        ofile.close();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
void write_cn_phased_bins( cn_map& chromosome_map, std::string outputFile, std::vector<int>& good_bins, std::vector<int>& merged_bins ) {
        ofstream ofile; ofile.open(outputFile);
        cout << "write output " << endl;
	int nbins = good_bins.size();
	int m = 0;
	double min_switch_prob = 0;
        for (auto ind : good_bins) {
		if (nbins > (m+1)) {
			int next = good_bins[m+1];
			min_switch_prob = calc_switch(chromosome_map,ind,next);
		}
                double hapA_cov,hapB_cov;
                int output_hap = chromosome_map[ind].switch_spin;
                if (output_hap == 1) {  hapA_cov = chromosome_map[ind].hapA_avg; hapB_cov = chromosome_map[ind].hapB_avg; }
                if (output_hap == -1) { hapA_cov = chromosome_map[ind].hapB_avg; hapB_cov = chromosome_map[ind].hapA_avg; }
                ofile << ind << "\t" << output_hap << "\t" << hapA_cov << "\t" << hapB_cov << "\t" << min_switch_prob << "\t" << merged_bins[m] << endl;
                m += 1;
        }
        ofile.close();
}













                        //cout << pos1 << "\t" << base1 << "\t" << pos2 << "\t" << base2 << "\t" << nconnections << endl;
                        //cout << pos1 << "\t" << pos2 << "\t" << " reads1" << endl;
                        //for (int j = 0; j < v1.size(); j++) { cout << v1[j] << " " << v1[j].length() << endl;}
                        //for ( std::vector<std::string>::iterator it5 = v1.begin(); it5 != v1.end(); ++it5 ) { cout << *it5 << " " << *it5.length() << endl; }
                        //cout << pos1 << "\t" << pos2 << "\t" << " reads2" << endl;
                        //for (int k = 0; k < v2.size(); k++) { cout << v2[k] << " " << v2[k].length() << endl;}
                        //for ( std::vector<std::string>::iterator it6 = v2.begin(); it6 != v2.end(); ++it6 ) { cout << *it6 << " " << *it6.length() << endl; }







                        //bool first = true;
                        //cout << "===" << endl;
                        //for (auto rn : read_names) { if (!first) { ofile << ','; } first = false; ofile << rn; }  //cout << rn << ",";
