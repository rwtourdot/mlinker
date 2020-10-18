#include "run_scaffold.h"

/////////////// namespaces ///////////////////
namespace opt {
        static std::string chr_choice = "chr20";
        static std::string technology = "tenx";
        static std::string input_hap_sol_file = "./output/hap_solution_trial_chr4.dat";
        static std::string input_hic_graph_file = "./output/graph_variant_sep25_hic_chr4.dat";
        static std::string id_string = "default";
	static int min_het = 5;
	static int min_links = 10;
	static int max_hic_phasing = 10000000;
        static int start_bound = 0;
        static int end_bound = 300000000;
        static int tenx_block_cutoff = -700;
        static bool cutoff_defined = false;
        static bool region_defined = false;
};

/////////////// structures ///////////////////
static const char* shortopts = "ho:i:g:c:e:n:";
static const struct option longopts[] = {
        { "help",        no_argument, NULL, 'h' },
        { "hap-file",    no_argument, NULL, 'i' },
        { "hic-file",    no_argument, NULL, 'g' },
        { "chr-choice",  no_argument, NULL, 'c' },
        { "tenx_blockE_cut",  no_argument, NULL, 'e' },
        { "id_string",   no_argument, NULL, 'n' }
};

///////////////////////////////////////////////////////////////////////////////////////////////////////
static const char *SCAFFOLD_USAGE_MESSAGE =
"Usage: mlinker scaffold [OPTION] -i /output/hap_solution_*chr4.dat -g /output/graph_hic_file_*chr4.dat -c chr4 \n\n"
"\n"
"********* this command is set up for hg38 *********"
"\n"
"  Options\n"
"  -i,      input haplotype solution file path \n"
"  -g,      input hic graph file \n"
"  -c,      chromosome name ( chr4 ) or contig name depending on bam \n"
"  -e,      tenx blockE cutoff (-700) \n"
"  -n,      id string for output files \n"
"\n";

///////////////////////////////////////////////////////////////////////////////////////////////////////
static void parse_scaffold_options( int argc, char** argv ) {
        bool die = false; //bool vcf_load = false; //bool cov_load = false;
        if (argc < 2) {
          std::cerr << "\n" << SCAFFOLD_USAGE_MESSAGE;
          exit(1);
        }
        if (string(argv[1]) == "help" || string(argv[1]) == "--help") {
                std::cerr << "\n" << SCAFFOLD_USAGE_MESSAGE;
                exit(1);
        }
        for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
                std::istringstream arg(optarg != NULL ? optarg : "");
                switch (c) {
                case 'h': die = true; break;
                case 'i': arg >> opt::input_hap_sol_file; break;
                case 'g': arg >> opt::input_hic_graph_file; break;
                case 'c': arg >> opt::chr_choice; break;
                case 'n': arg >> opt::id_string; break;
                case 'e': opt::cutoff_defined = true; arg >> opt::tenx_block_cutoff; break;
                }
        }
        //if (!opt::vcf_load && !opt::cov_load) { die = true; }
        if (die) {
          std::cerr << "\n" << SCAFFOLD_USAGE_MESSAGE;
          exit(1);
        }
        if (opt::tenx_block_cutoff > 0) { int temp_bc = opt::tenx_block_cutoff; opt::tenx_block_cutoff = -temp_bc; }
        //parse_region_string( opt::chr_choice );
        cout << endl;
        cout << "############### running scaffold phaser ############### " << endl;
        cout << "== chromosome === " << opt::chr_choice << endl;
        cout << "== input hap sol file  === " << opt::input_hap_sol_file << endl;
        cout << "== input hic graph file  === " << opt::input_hic_graph_file << endl;
        cout << "== id string  === " << opt::id_string << endl;
	if (opt::cutoff_defined) { cout << "== blockE cutoff === " << opt::tenx_block_cutoff << endl; }
        cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx " << endl;
        cout << endl;
};

void init_block_matrix( int& num_block, map_matrix<int>& block_matrix, map_matrix<int>& total_matrix, block_dictionary& bl_dict, coord_dictionary& pdict, hic_link hlink, variant_graph& hic_vgraph, int centromere_pos ) {
        static const std::size_t length = pdict.num_paired;
        std::unordered_map<int,int> pos_to_block;
        std::unordered_map<int,int> pos_to_index;
        int loop_p;
        for (int i = 0; i < length; i++) {
                loop_p = pdict.double_positions[i];
                num_block = pdict.block[i];
                pos_to_block[loop_p] = num_block;
                pos_to_index[loop_p] = i;  
		cout << loop_p << "\t" << num_block << "\t" << i << endl;
        }
        map_matrix<int> temp_block_matrix(num_block);
        map_matrix<int> temp_total_matrix(num_block);
        // ############################### fill block matrix #################################
        for (int i = 0; i < length; i++) {
                int pos1 = pdict.double_positions[i];
		std::string arm1 = "p"; if ( pos1 > centromere_pos ) { arm1 = "q"; }
                int hap1 = pdict.haplotype[i];
                int block1 = pdict.block[i];
                std::string ref_hash = pdict.ref_handle[i];
                std::string alt_hash = pdict.alt_handle[i];
		cout << pos1 << " " << ref_hash << " " << alt_hash << endl;
                bl_dict.add_bl_het(block1,pos1,hap1,ref_hash,alt_hash);
		for (auto& it : hic_vgraph[ref_hash].connections) { //for (auto it2 : it.second.connections)
			if ( hic_vgraph[ref_hash].unique_total > 0 ) {
			int alt1 = 1;
			int pos2 = std::stoi(split_string_first(it.first,"_",0));
			std::string hic_ref_base = split_string_first(split_string_first(it.first,"_",1),"_",0);
			std::string hic_read_base = split_string_first(split_string_first(it.first,"_",1),"_",1);
			std::string arm2 = "p"; if ( pos2 > centromere_pos ) { arm2 = "q"; }
			cout << "in loop " << ref_hash << " " << it.first << endl;
			int depth = it.second;
			int diff = std::abs(pos2-pos1);
                        int block2 = pos_to_block[pos2];
                        int hap2 = pdict.haplotype[pos_to_index[pos2]];
                        int alt2 = -1; if (hic_ref_base == hic_read_base) { alt2 = 1; }
			if ( block1 != 0 & block2 != 0 ) {  
                                int link_hap_old = alt1*hap1*alt2*hap2;
                                int link_hap = alt1*hap1*alt2*hap2*depth;
                                temp_total_matrix.add_to(block1,block2,link_hap);
                                temp_total_matrix.add_to(block2,block1,link_hap);
                                if ( diff < opt::max_hic_phasing & block1 != block2 & arm1 == arm2 ) {
                                        temp_block_matrix.add_to(block1,block2,link_hap);
                                        temp_block_matrix.add_to(block2,block1,link_hap);
					bl_dict.bl_nlinks[block1] += 1;
					bl_dict.bl_nlinks[block2] += 1;
				}
			} }
		}
                for (auto& it : hic_vgraph[alt_hash].connections) { //for (auto it2 : it.second.connections)
			if ( hic_vgraph[alt_hash].unique_total > 0 ) {
                        int alt1 = -1;
                        int pos2 = std::stoi(split_string_first(it.first,"_",0));
                        std::string hic_ref_base = split_string_first(split_string_first(it.first,"_",1),"_",0);
                        std::string hic_read_base = split_string_first(split_string_first(it.first,"_",1),"_",1);
			std::string arm2 = "p"; if ( pos2 > centromere_pos ) { arm2 = "q"; }
                        int depth = it.second;
                        int diff = std::abs(pos2-pos1);
                        int block2 = pos_to_block[pos2];
                        int hap2 = pdict.haplotype[pos_to_index[pos2]];
                        int alt2 = -1; if (hic_ref_base == hic_read_base) { alt2 = 1; }
                        if ( block1 != 0 & block2 != 0 ) {  //cout << pos2 << "\t" << block1 << "\t" << block2 << endl;
                                int link_hap_old = alt1*hap1*alt2*hap2;
                                int link_hap = alt1*hap1*alt2*hap2*depth;
                                temp_total_matrix.add_to(block1,block2,link_hap);
                                temp_total_matrix.add_to(block2,block1,link_hap);
                                if ( diff < opt::max_hic_phasing & block1 != block2 & arm1 == arm2 ) {
                                        //cout << block1 << "\t" << block2 << "\t" << ref_hash << "\t" << j << "\t" << link_hap << "\t" << link_hap_old << endl;
                                        temp_block_matrix.add_to(block1,block2,link_hap);
                                        temp_block_matrix.add_to(block2,block1,link_hap);
                                        bl_dict.bl_nlinks[block1] += 1;
                                        bl_dict.bl_nlinks[block2] += 1;
                                }
                        } }
                }
        }
	total_matrix = temp_total_matrix;
	block_matrix = temp_block_matrix;
};

void pq_phasing(int num_block, map_matrix<int> total_matrix, block_dictionary& bl_dict, coord_dictionary& pdict, hic_link hlink, int centromere_pos ) {
        std::unordered_map<int,int> pos_to_index;
        for (int i = 0; i < pdict.num_paired; i++) { pos_to_index[pdict.double_positions[i]] = i; }
	int flip = 0; int dont_flip = 0;
        for (int i = 0; i < pdict.num_paired; i++) {
                int pos1 = pdict.double_positions[i];
                int hap1 = pdict.haplotype[i];
                int block1 = pdict.block[i];
                std::string ref_hash = pdict.ref_handle[i];
                std::string alt_hash = pdict.alt_handle[i];
		std::string arm1 = "p"; if ( pos1 > centromere_pos ) { arm1 = "q"; }
		if ( std::find( bl_dict.subset_blocks.begin(), bl_dict.subset_blocks.end(), block1 ) != bl_dict.subset_blocks.end() )  {
		//cout << pos1 << "\t" << block1 << "\t" << arm1 << "\t" << bl_dict.block_map[block1] << "\t" << bl_dict.length_subset << endl;
		int block_hap1 = bl_dict.subset_haplotype[bl_dict.block_map[block1]];
		if ( hlink.hic_link_map.find(ref_hash) != hlink.hic_link_map.end() ) {
                        int alt1 = 1;
                        for( int j = 0; j < hlink.hic_link_map[ref_hash].size(); j++ ) {
                                //cout << ref_hash << "\t" << j << "\t" << hlink.hic_link_map[ref_hash][j] << endl;
                                std::string hic_link_string = hlink.hic_link_map[ref_hash][j];
                                int depth = hlink.hic_link_depth_map[ref_hash][j];
                                int pos2 = std::stoi(split_string_first(hic_link_string,"_",0));  //int diff = std::abs(pos2-pos1);
                                std::string hic_ref_base = split_string_first(split_string_first(hic_link_string,"_",1),"_",0);
                                std::string hic_read_base = split_string_first(split_string_first(hic_link_string,"_",1),"_",1);
				std::string arm2 = "p"; if ( pos2 > centromere_pos ) { arm2 = "q"; }
                                int block2 = pdict.block[pos_to_index[pos2]];
				if ( std::find( bl_dict.subset_blocks.begin(), bl_dict.subset_blocks.end(), block2 ) != bl_dict.subset_blocks.end() )  {
					int block_hap2 = bl_dict.subset_haplotype[bl_dict.block_map[block2]];
                                	int hap2 = pdict.haplotype[pos_to_index[pos2]];
                                	int alt2 = -1; if (hic_ref_base == hic_read_base) { alt2 = 1; }
					//cout << hic_link_string << "\t" << pos2 << "\t" << hic_ref_base << "\t" << hic_read_base << "\t" << alt2 << endl;
					int link_hap = alt1*block_hap1*hap1*alt2*block_hap2*hap2;
					//int link_hap = alt1*block_hap1*hap1*alt2*block_hap2*hap2*depth;
					if (arm1 != arm2) {
						//cout << ref_hash << "\t" << pos1 << "\t" << arm1 << "\t" << block1 << "\t" << pos2 << "\t"<< arm2 << "\t" << block2 << "\t"<< link_hap << endl;
						if (link_hap == -1) { flip += depth; }
						if (link_hap == 1) { dont_flip += depth; }
					}
				}
			}
		}
		if ( hlink.hic_link_map.find(alt_hash) != hlink.hic_link_map.end() ) {
                        int alt1 = -1;
                        for( int j = 0; j < hlink.hic_link_map[alt_hash].size(); j++ ) {
                                //cout << ref_hash << "\t" << j << "\t" << hlink.hic_link_map[ref_hash][j] << endl;
                                std::string hic_link_string = hlink.hic_link_map[alt_hash][j];
                                int depth = hlink.hic_link_depth_map[alt_hash][j];
                                int pos2 = std::stoi(split_string_first(hic_link_string,"_",0));  //int diff = std::abs(pos2-pos1);
                                std::string hic_ref_base = split_string_first(split_string_first(hic_link_string,"_",1),"_",0);
                                std::string hic_read_base = split_string_first(split_string_first(hic_link_string,"_",1),"_",1);
                                std::string arm2 = "p"; if ( pos2 > centromere_pos ) { arm2 = "q"; }
                                int block2 = pdict.block[pos_to_index[pos2]];
                                if ( std::find( bl_dict.subset_blocks.begin(), bl_dict.subset_blocks.end(), block2 ) != bl_dict.subset_blocks.end() )  {
                                        int block_hap2 = bl_dict.subset_haplotype[bl_dict.block_map[block2]];
                                        int hap2 = pdict.haplotype[pos_to_index[pos2]];
                                        int alt2 = -1; if (hic_ref_base == hic_read_base) { alt2 = 1; }
					//cout << hic_link_string << "\t" << pos2 << "\t" << hic_ref_base << "\t" << hic_read_base << "\t" << alt2 << endl;
                                        int link_hap = alt1*block_hap1*hap1*alt2*block_hap2*hap2;
                                        //int link_hap = alt1*block_hap1*hap1*alt2*block_hap2*hap2*depth;
                                        if (arm1 != arm2) {
                                                //cout << alt_hash << "\t" << pos1 << "\t" << arm1 << "\t" << block1 << "\t" << pos2 << "\t"<< arm2 << "\t" << block2 << "\t" << link_hap << endl;
						if (link_hap == -1) { flip += depth; }
						if (link_hap == 1) { dont_flip += depth; }
                                        }
                                }
                        }
                }
		}
	}
	cout << "flip: " << flip << endl;
	cout << "dont flip: " << dont_flip << endl;
	///////// flip pq arms if needed
	if (flip > dont_flip) {
		cout << "flip pq haplotype" << endl;
        	for (int l = 0; l < bl_dict.length_subset; l++) {
			int block_position = bl_dict.bl_min[bl_dict.block_map_inverted[l]];
			if ( block_position > centromere_pos ) {
				//cout << "flip block: " << bl_dict.block_map_inverted[l] << "\t" << l << "\t" << block_position << endl;
				bl_dict.subset_haplotype[l] = bl_dict.subset_haplotype[l]*-1;
			}
        	}
	}
	//if (flip > dont_flip) {
	//	cout << "flip pq haplotype" << endl;
	//	for (int i = 0; i < pdict.num_paired; i++) {
	//		int pos1 = pdict.double_positions[i];   //int hap1 = pdict.haplotype[i];
        //        	int block1 = pdict.block[i];
	//		if ( pos1 > centromere_pos ) {
	//			bl_dict.subset_haplotype[bl_dict.block_map[block1]] = bl_dict.subset_haplotype[bl_dict.block_map[block1]]*-1;
	//		}
	//	}
	//}
}

void create_hic_links( variant_graph& hic_vgraph, hic_link& hlink ) {
	for (auto& it : hic_vgraph) { it.second.unique_hash(); }
        int i=0;
	for (auto& it : hic_vgraph) {
                int pos1 = it.second.pos;
                std::string var1 = it.first;   
		//cout << pos1 << "\t" << var1 << "\t" << it.second.unique_total << endl;
                if (it.second.unique_total > 0) {
                for (auto it2 : it.second.connections) {  //int pos2 = 100;  //var_dict[it2.first].pos;
                        int pos2 = hic_vgraph[it2.first].pos;
                        int depth = it2.second;
                        std::string var2 = it2.first;
			//cout << pos1 << "\t" << var1 << "\t" << pos2 << "\t" << var2 << "\t" << depth << endl;
			hlink.add_link(pos1,pos2,var1,var2,depth);  //anchor_pos1.push_back(pos1);
                	i++;
		} }
	}
        hlink.nlinks = i;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
void run_scaffold( int argc, char** argv ) {
        parse_scaffold_options(argc, argv);
        std::string scaffoldsolutionFile = "./output/hap_full_scaffold_" + opt::id_string + "_" + opt::chr_choice + ".dat";
        std::string scaffoldEnergyFile = "./output/block_scaffold_energy_" + opt::id_string + "_" + opt::chr_choice + ".dat";
        std::string bmatrixFile = "./output/block_phase_matrix_" + opt::id_string + "_" + opt::chr_choice + ".dat";
        cout << "- output scaffold file: " << scaffoldsolutionFile << endl;
        cout << "- output energy file: " << scaffoldEnergyFile << endl;
        cout << "- output matrix file: " << bmatrixFile << endl;
        cout << endl;
	std::unordered_map<std::string,int> hg38_pq = {{"chr1",123400000},{"chr2",93900000},{"chr3",90900000},{"chr4",50000000},{"chr5",48800000},{"chr6",59800000},{"chr7",60100000},{"chr8",45200000},{"chr9",43000000},{"chr10",39800000},{"chr11",53400000},{"chr12",35500000},{"chr13",17700000},{"chr14",17200000},{"chr15",19000000},{"chr16",36800000},{"chr17",25100000},{"chr18",18500000},{"chr19",26200000},{"chr20",28100000},{"chr21",12000000},{"chr22",15000000},{"chrX",61000000},{"chrY",10400000}};
	//std::unordered_map<std::string,int> hg19_pq = {{"chr1",125000000},{"chr2",93300000},{"chr3",91000000},{"chr4",50400000},{"chr5",48400000},{"chr6",61000000},{"chr7",59900000},{"chr8",45600000},{"chr9",49000000},{"chr10",40200000},{"chr11",53700000},{"chr12",35800000},{"chr13",17900000},{"chr14",17600000},{"chr15",19000000},{"chr16",36600000},{"chr17",24000000},{"chr18",17200000},{"chr19",26500000},{"chr20",27500000},{"chr21",13200000},{"chr22",14700000},{"chrX",60600000},{"chrY",12500000}};
	int centromere_pos = hg38_pq[opt::chr_choice];
        coord_dictionary pdict;
	block_dictionary bl_dict;
	// ############################### load hap solution and hic links ###########################
	read_hap_solution_initialize( opt::input_hap_sol_file, pdict );
	//
	hic_link hlink;
	read_graph hic_rgraph;
	variant_graph hic_vgraph;
	read_variant_graph_file( opt::input_hic_graph_file, opt::chr_choice, hic_vgraph, hic_rgraph );
	link_hashes( hic_vgraph, hic_rgraph );
	calc_coverage_unique_hash( hic_vgraph );
	create_hic_links( hic_vgraph, hlink );
	//
	//read_hic_links_file( opt::input_hic_links_file, hlink );
	call_blocks(pdict,opt::tenx_block_cutoff);
	int num_block;
        map_matrix<int> block_matrix;
        map_matrix<int> total_matrix;
	init_block_matrix( num_block, block_matrix, total_matrix, bl_dict, pdict, hlink, hic_vgraph, centromere_pos );
	// ############################### initialize block matrix #################################
	bl_dict.subset_initialize( num_block, opt::min_het, opt::min_links, block_matrix );
	bl_dict.subset_hap_random_initialization();
	//#############################################
	//std::vector<int> hic_limit_loop = {100000,500000,1000000,5000000,10000000};
	std::vector<int> hic_limit_loop = {50000,100000,500000,1000000,5000000,10000000};
	//std::vector<int> hic_limit_loop = {10000000};
	solver_recursive_hic( bl_dict, hic_limit_loop, block_matrix );  //limit_matrix
	//#############################################
        cout << " solver finished " << endl;
	pq_phasing( num_block, total_matrix, bl_dict, pdict, hlink, centromere_pos );
	cout << " phasing pq arm chr: " << opt::chr_choice << endl;
	cout << " centromere position: " << centromere_pos << endl;
	//#############################################
        ///////////// write haplotype output
        cout << " writing scaffold file: " << scaffoldsolutionFile << endl;
	write_scaffold( opt::chr_choice, scaffoldsolutionFile, bl_dict, pdict, centromere_pos );
        cout << " writing energy file: " << scaffoldEnergyFile << endl;
	write_scaffold_energy( opt::chr_choice, scaffoldEnergyFile, bl_dict, centromere_pos );
	write_block_phasing_matrix( bl_dict, bmatrixFile );
        return;
};
