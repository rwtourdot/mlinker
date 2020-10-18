#include "read_linker_output.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////
bool is_file_exist(std::string filename ) {
	std::ifstream infile(filename);
	return infile.good();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
std::string split_string_first( std::string s, std::string delimiter, int choice ) {
	std::string return_string;
	if ( choice == 0 ) { return_string = s.substr(0, s.find(delimiter)); }
	if ( choice == 1 ) { return_string = s.substr(s.find(delimiter)+1, s.length()); }
	return return_string;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
void read_bin_haplotype_cn_data( std::string hap_cn_file, cn_map& chromosome_map, int binsize, std::vector<int>& bin_array ) {
        if (!is_file_exist(hap_cn_file)) { cout << " file not found " << hap_cn_file << endl; exit(1); }
        std::ifstream inputFile(hap_cn_file); std::string line;
	int k = 0;
        while (getline(inputFile,line)) {
		cout << line << endl;
		std::istringstream ss(line);
		std::string chrom_contig,pos_start_str,pos_end_str,hap1_cn_str,hap2_cn_str,h1_mean_frac_str;
		ss >> chrom_contig >> pos_start_str >> pos_end_str >> hap1_cn_str >> hap2_cn_str >> h1_mean_frac_str;
		if ( chrom_contig != "CONTIG" ) {
                        int pos_start = std::stoi(pos_start_str);
                        int pos_end = std::stoi(pos_end_str);
                        double hap1_cn = std::stod(hap1_cn_str);
                        double hap2_cn = std::stod(hap2_cn_str);
                        double h1_frac = std::stod(h1_mean_frac_str);
			if ( pos_end % binsize == 0) { //cout << j << endl;
				int bin = pos_end/binsize;
				cn_bin loop_bin;
				loop_bin.initialize(bin);
				chromosome_map[bin] = loop_bin;
				bin_array.push_back(bin);
				chromosome_map[bin].enter_coverage(hap1_cn,hap2_cn,bin);
				cout << pos_end << "\t" << bin << endl;
			}
		}
        }
	return;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
void read_het_coverage( std::string coverageFile, variant_graph& vgraph ) {
	if (!is_file_exist(coverageFile)) { cout << " file not found " << coverageFile << endl; exit(1); }
        std::ifstream inputFile(coverageFile); std::string line;
	while (getline(inputFile,line)) {
		std::istringstream ss(line);
		int pos,tot_cov;
		std::string het_name,ref_var,dbase,abase,ibase,tbase,cbase,gbase;
		ss >> het_name >> pos >> ref_var >> ibase >> dbase >> gbase >> cbase >> abase >> tbase >> tot_cov;
		//cout << coverageFile << "\t" << pos << " " << abase << " " << tbase << " " << cbase << " " << gbase << endl;
		std::string ref = split_string_first(ref_var,":",0);
		std::string var = split_string_first(ref_var,":",1);
		int anum = std::stoi(split_string_first(abase,"|",1));
		int tnum = std::stoi(split_string_first(tbase,"|",1));
		int cnum = std::stoi(split_string_first(cbase,"|",1));
		int gnum = std::stoi(split_string_first(gbase,"|",1));
		int dnum = std::stoi(split_string_first(dbase,"|",1));
		int inum = std::stoi(split_string_first(ibase,"|",1));
		std::string variant_id = std::to_string(pos)+"_"+ref+"_"+var;
		std::string reference_id = std::to_string(pos)+"_"+ref+"_"+ref;
		//cout << pos << " " << anum << " " << tnum << " " << cnum << " " << gnum << endl;
		//cout << variant_id << " " << reference_id << " " << het_name << "  " << pos << "  " << ref_var << "  " << ref << "  " << var << "  " << tot_cov << "  " << anum << endl;
                variant_node v_node1,v_node2;
                vgraph[variant_id] = v_node1;
		vgraph[reference_id] = v_node2;
                vgraph[variant_id].reset_values(pos,true,var,ref,anum,tnum,cnum,gnum,dnum,inum);
                vgraph[reference_id].reset_values(pos,false,var,ref,anum,tnum,cnum,gnum,dnum,inum);
		vgraph[variant_id].paired = true;
		vgraph[reference_id].paired = true;
	}
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<vcf_entry> read_het_coverage_vvec( std::string coverageFile, int chromosome ) {
	if (!is_file_exist(coverageFile)) { cout << " file not found " << coverageFile << endl; exit(1); }
	int i = 0;
	std::vector<vcf_entry> entry;
        std::ifstream inputFile(coverageFile); std::string line;
	while (getline(inputFile,line)) {
                std::istringstream ss(line);
                int pos,tot_cov;
                std::string het_name,ref_var,dbase,abase,ibase,tbase,cbase,gbase;
                ss >> het_name >> pos >> ref_var >> ibase >> dbase >> gbase >> cbase >> abase >> tbase >> tot_cov;
                //cout << coverageFile << "\t" << pos << " " << abase << " " << tbase << " " << cbase << " " << gbase << endl;
                std::string ref = split_string_first(ref_var,":",0);
                std::string var = split_string_first(ref_var,":",1);
                int anum = std::stoi(split_string_first(abase,"|",1));
                int tnum = std::stoi(split_string_first(tbase,"|",1));
                int cnum = std::stoi(split_string_first(cbase,"|",1));
                int gnum = std::stoi(split_string_first(gbase,"|",1));
                int dnum = std::stoi(split_string_first(dbase,"|",1));
                int inum = std::stoi(split_string_first(ibase,"|",1));
                std::string variant_id = std::to_string(pos)+"_"+ref+"_"+var;
                std::string reference_id = std::to_string(pos)+"_"+ref+"_"+ref;
		entry.push_back(vcf_entry());
		entry[i].pos = pos;    entry[i].chromosome_id = chromosome;
		entry[i].ref_base = ref;   entry[i].var_base = var;
		entry[i].bounded = true;
		i += 1;
	}
        return entry;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
void read_hap_solution( std::string hapfile, coord_dictionary& pdict ) {
        if (!is_file_exist(hapfile)) { cout << " file not found " << hapfile << endl; exit(1); }
	std::unordered_map<int,int> index_map;
	for (int i = 0; i < pdict.num_paired; i++) { index_map[pdict.sorted_paired_positions[i]] = i; }
        std::ifstream inputfile(hapfile); std::string line;
        while (getline(inputfile,line)) {
                std::istringstream ss(line);
                int index,pos,hap,hapa_cov,hapb_cov,block,width;
                std::string ref_name,het_name,ref_var,dbase,abase,ibase,tbase,cbase,gbase,spine,switche;
                ss >> index >> pos >> ref_name >> het_name >> hap >> hapa_cov >> hapb_cov >> spine >> switche >> block >> width;
                std::string chrom = split_string_first(ref_name,"_",0);
		int pind = index_map[pos];
		pdict.haplotype[pind] = hap;
		pdict.block[pind] = block;
		pdict.deltaE[pind] = std::stod(spine);
		pdict.switchE[pind] = std::stod(switche);
		pdict.reload_bool[pind] = true;   //cout << pdict.reload_bool[pind] << endl;
                cout << index << "  " << pind << "  " << pos << "  " << ref_name << "  " << chrom << "  " << hap << "  " << spine << "  " << switche << " "  << block << endl;
        }
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
void read_hap_solution_initialize( std::string hapfile, coord_dictionary& pdict ) {
        if (!is_file_exist(hapfile)) { cout << " file not found " << hapfile << endl; exit(1); }
        //std::unordered_map<int,int> index_map;
        //for (int i = 0; i < pdict.num_paired; i++) { index_map[pdict.sorted_paired_positions[i]] = i; }
	int num_paired = 0;
        std::ifstream inputfile(hapfile); std::string line;
        while (getline(inputfile,line)) {
                std::istringstream ss(line);
                int index,pos,hap,hapa_cov,hapb_cov,block,width;
                std::string ref_name,het_name,ref_var,dbase,abase,ibase,tbase,cbase,gbase;
		double spine,switche;
                ss >> index >> pos >> ref_name >> het_name >> hap >> hapa_cov >> hapb_cov >> spine >> switche >> block >> width;
                std::string chrom = split_string_first(ref_name,"_",0);
                std::string ref_hash = split_string_first(ref_name,"_",1);
                std::string alt_hash = split_string_first(het_name,"_",1);
		//ref_name.back()
		//het_name.back()
		pdict.ref_handle.push_back(ref_hash);
		pdict.alt_handle.push_back(alt_hash);
		pdict.double_positions.push_back(pos);
		pdict.all_positions.push_back(pos);
        	pdict.deltaE.push_back(spine);
        	pdict.switchE.push_back(switche);
        	pdict.span_bound.push_back(0);
       		pdict.block.push_back(block);
       		pdict.haplotype.push_back(hap);
        	//pdict.within_filter.push_back(false);
		num_paired += 1;
                //cout << index << "  " << pos << "  " << ref_name << "  " << chrom << "  " << hap << "  " << spine << "  " << switche << " "  << block << endl;
        }
	pdict.num_paired = num_paired;
	pdict.num_total = num_paired;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
void read_scaffold( std::string scaffile, coord_dictionary& pdict ) {
        if (!is_file_exist(scaffile)) { cout << " file not found " << scaffile << endl; exit(1); }
        //std::unordered_map<int,int> index_map;
        //for (int i = 0; i < pdict.num_paired; i++) { index_map[pdict.sorted_paired_positions[i]] = i; }
        int num_paired = 0;
        std::ifstream inputfile(scaffile); std::string line;
        while (getline(inputfile,line)) {
                std::istringstream ss(line);
                int pos,hap,block_hap,block,width;
                std::string ref_name,alt_name,arm,chrom,block_str,tenx_str;
                double tenx_flipE,tenx_blockE;
                ss >> chrom >> arm >> pos >> ref_name >> alt_name >> hap >> block_str >> block_hap >> block >> tenx_str >> tenx_flipE >> tenx_blockE;
                pdict.ref_handle.push_back(ref_name);
                pdict.alt_handle.push_back(alt_name);
                pdict.double_positions.push_back(pos);
                pdict.all_positions.push_back(pos);
                pdict.deltaE.push_back(tenx_flipE);
                pdict.switchE.push_back(tenx_blockE);
                pdict.span_bound.push_back(0);
                pdict.block.push_back(block);
                pdict.haplotype.push_back(hap);
                //pdict.within_filter.push_back(false);
                num_paired += 1;
                //cout << index << "  " << pos << "  " << ref_name << "  " << chrom << "  " << hap << "  " << spine << "  " << switche << " "  << block << endl;
        }
        pdict.num_paired = num_paired;
        pdict.num_total = num_paired;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
void read_sv_file( std::string svfile, std::vector<sv_entry>& sv_list ) {
        if (!is_file_exist(svfile)) { cout << " file not found " << svfile << endl; exit(1); }
        std::ifstream inputfile(svfile); std::string line;
	int i=0;
        while (getline(inputfile,line)) {
		sv_entry sv_temp;
                std::istringstream ss(line);
                int raindex,pos1,str1,pos2,str2,totalcount;
                std::string chr1,chr2;
                ss >> raindex >> chr1 >> pos1 >> str1 >> chr2 >> pos2 >> str2 >> totalcount;
		if ( chr1.length() > 1 ) {
			if ( totalcount > min_sv_reads ) {
				sv_temp.set_locations(chr1,pos1,chr2,pos2,str1,str2,totalcount);
				sv_list.push_back(sv_temp);
                		//cout << chr1 << "  " << pos1 << "  " << str1 << "  " << chr2 << "  " << pos2 << "  " << str2 << "  " << totalcount << endl;
			}
		}
		//if ( i > 40 ) { break; }
		i++;
        }
};


////////////////////////////////////////////////////////////////////////////////////////////////////////
void read_hic_links_file( std::string hic_links_file, hic_link& hlink ) {
        if (!is_file_exist(hic_links_file)) { cout << " file not found " << hic_links_file << endl; exit(1); }
        std::ifstream inputfile(hic_links_file); std::string line;
        int i=0;
        while (getline(inputfile,line)) {  //sv_entry sv_temp;
                std::istringstream ss(line);
                int pos1,pos2,depth;
		std::string var1,var2,link_readnames;
                ss >> pos1 >> var1 >> pos2 >> var2 >> depth >> link_readnames;
		//cout << pos1 << "\t" << pos2 << endl;
		hlink.add_link(pos1,pos2,var1,var2,depth);  //anchor_pos1.push_back(pos1);
                i++;
        }
	hlink.nlinks = i;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
void read_variant_graph_file( std::string graphFile, std::string chromosome, variant_graph& vgraph, std::unordered_map<std::string,read_tree>& rgraph ) {
	if (!is_file_exist(graphFile)) { cout << " file not found " << graphFile << endl; exit(1); }
        std::ifstream inputFile(graphFile); std::string line;
        int i=0;
        while (getline(inputFile,line)) {
                std::istringstream ss(line);
                int pos;
		std::string het_hash,chromosome,hash,readname;
                ss >> het_hash >> chromosome >> pos >> hash >> readname;
		//cout << het_hash << "\t" << chromosome << "\t" << hash << "\t" << readname << endl;
		if (rgraph.find(hash) != rgraph.end()) { rgraph[hash].add_connection(het_hash,readname); }
		else { rgraph[hash].set_name(hash); rgraph[hash].add_connection(het_hash,readname); }
		//############################

		//
		bool present1 = get_tag_value(het_hash,vgraph);
		bool is_variant = false;
		std::string temp_string = split_string_first(het_hash,"_",1);
		std::string ref_base = split_string_first(temp_string,"_",0);
		std::string read_base = split_string_first(temp_string,"_",1);
		//cout << pos << "\t" << ref_base << "\t" << read_base << endl;
		char ref_base_char = ref_base[0];
		char read_base_char = read_base[0];
		if (ref_base != read_base) { is_variant = true; }
                if (!present1) {
                        variant_node v_node;
                        vgraph[het_hash] = v_node;
                        vgraph[het_hash].set_values(pos,is_variant,read_base,ref_base);
                }
		vgraph[het_hash].add_base_to_dict(read_base_char,hash);
		vgraph[het_hash].add_connected_read(hash,readname);
	}
};



                //int pind = index_map[pos];
                //cout << index << endl;


                //pdict.haplotype[pind] = hap;
                //pdict.block[pind] = block;
                //pdict.deltaE[pind] = std::stod(spine);
                //pdict.switchE[pind] = std::stod(switche);
                //pdict.reload_bool[pind] = true;   //cout << pdict.reload_bool[pind] << endl;

//void read_hap_solution_file( std::string hap_solutionFile, coord_dictionary& pdict ) {
//	if (!is_file_exist(hap_solutionFile)) { cout << " file not found " << hap_solutionFile << endl; exit(1); }
//	std::ifstream inputFile(hap_solutionFile); std::string line;
//        int i=0;
//        while (getline(hap_solutionFile,line)) {
//	}
//}

        //std::string ref_base,var_base;
        //int rec_position,rec_id,len_ref;
        //bcf_srs_t *sr = bcf_sr_init();
        //bcf_sr_add_reader(sr,input_vcf_file.c_str());
        //bcf1_t *line;
        //int i = 0;
        //while(bcf_sr_next_line(sr)) {
        //        line = bcf_sr_get_line(sr,0);
        //        if (line->rid == chromosome) {
        //                if ((DEBUG == 1) && (i > het_cutoff)) { break; }
        //                ref_base = line->d.allele[0];
        //                var_base = line->d.allele[1];
        //                rec_position = line->pos;
        //                rec_id = line->rid;
        //                len_ref = line->rlen;
        //                entry.push_back(vcf_entry());
        //                entry[i].pos = rec_position;    entry[i].chromosome_id = rec_id;
        //                entry[i].ref_base = ref_base;   entry[i].var_base = var_base;
        //                entry[i].bounded = true;
        //                i += 1;
        //        }
        //}
        //return entry;
//};
