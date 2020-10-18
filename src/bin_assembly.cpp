#include "bin_assembly.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////
void contig_node::initialize( std::unordered_map<int,gen_bin> gen_map ) {
    std::vector<std::string> temp_total_bx;
    int maximum_bin = 0;
    for (auto& it : gen_map) {
        if (it.first > maximum_bin) { maximum_bin = it.first; }
        for (int l=0; l < it.second.num_unique_bx; l++) {
            num_each_bx[it.second.unique_bx[l]] += 1;
            reverse_bx_map[it.second.unique_bx[l]].push_back(it.first);
            temp_total_bx.push_back(it.second.unique_bx[l]);
        }
    }       //cout << " maximum " << maximum_bin << endl;
    num_bins = maximum_bin;
    std::set<std::string> temp_cnx_bx(temp_total_bx.begin(),temp_total_bx.end());
    cnx_bx = temp_cnx_bx;
    num_bx = cnx_bx.size();
}


////////////////////////////////////////////////////////////////////////////////////////////////////////
static void init_bins(std::unordered_map<int,gen_bin> &gen_map, int contig_length, int binsize) {
        int num_bins = std::floor(contig_length/binsize);
        for (int i=0; i < num_bins; i++) {
                gen_bin bxnode;
                bxnode.set_values((int)i*binsize);
                gen_map[i] = bxnode;
        }
	return;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
static void init_bins_genome( full_map& chr_map, std::unordered_map<std::string,int>& contig_dict, int binsize) {	
	for (auto& it : contig_dict) {
		bx_map gen_map;
		int num_bins = std::floor(it.second/binsize);
        	for (int i=0; i < num_bins; i++) {
                	gen_bin bxnode;
               		bxnode.set_values((int)i*binsize);
                	gen_map[i] = bxnode;
        	}
		chr_map[it.first] = gen_map;
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
static void link_bins_genome( full_map& chr_map, std::unordered_map<std::string,read_tree> &bx_graph ) {
        int i = 0;
        clock_t t;
        t = clock();
	for (auto& it : chr_map) {
        	for (auto& it2 : it.second) {     //cout << it.first << endl;
                	for (int l=0; l < it2.second.connected_bx.size(); l++) {
                        	//cout << it.first << " bin connected " <<  it.second.connected_bx[l] << "   " << bx_graph[it.second.connected_bx[l]].connected_bins.size() << endl;
                        	for (int m=0; m < bx_graph[it2.second.connected_bx[l]].connected_ints.size(); m++) {
                                	int connected_bin = bx_graph[it2.second.connected_bx[l]].connected_ints[m];    //.connected_het_strings[m];
					std::string connected_chrom = bx_graph[it2.second.connected_bx[l]].connected_chr[m];
					chr_map[it.first][it2.first].add_connection(connected_chrom,connected_bin);
                                	//if (connected_bin != it2.first) {	 
					//		it2.second.add_connection(connected_bin); 
					//}
                        	}
                	}
                	i++;
        	}
	}
        //for (auto& it : var_dict) { it.second.count_connections(); }
        t = clock() - t;
        cout << "connected hashes ----- time: " << t << endl;
        return;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
static void create_contig( std::unordered_map<int,gen_bin> &gen_map, int contig_length ) {
        for (auto& it : gen_map) { it.second.find_unique_bx(); }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
static std::unordered_map<std::string,int> select_contigs( std::string contig_str, std::unordered_map<std::string,int>& contig_dict ) {
	std::unordered_map<std::string,int> subset_dict;
	//std::vector<std::string> all_chrom = {"chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX"};
	std::vector<std::string> all_chrom = {"chr17","chr18","chr19","chr20","chr21","chr22","chrX"};
        cout << endl;
        cout << "choice.chrom.sizes file:" << endl;
	for (auto& it : contig_dict) {
        	if (contig_str == "all") {
			if ( std::find(all_chrom.begin(), all_chrom.end(), it.first) != all_chrom.end() ) {
				subset_dict[it.first] = it.second;
				cout << it.first << "  " << it.second << endl;
			}
		}
		else {
			if ( it.first == contig_str ) {
				subset_dict[it.first] = it.second;
        			cout << it.first << "  " << it.second << endl;
			}
		}	
	}
        cout << endl;
	return subset_dict;	
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
void create_contig_dict(std::string input_bam_file,std::unordered_map<std::string,int> &contig_dict, std::unordered_map<std::string,contig_node> &contig_bx, std::string contig_str, int binsize, std::string output_file, bx_map& gen_map, full_map& chr_map ) {
	read_graph bx_graph;
	std::unordered_map<std::string,int> subset_dict;
	subset_dict = select_contigs( contig_str, contig_dict );
	init_bins_genome( chr_map, subset_dict, binsize );
        init_bins(gen_map,contig_dict[contig_str],binsize);
        connect_assembly_bam_pileup(input_bam_file,gen_map,bx_graph,binsize,chr_map,subset_dict);
       	//link_bins(gen_map,bx_graph);
	link_bins_genome(chr_map,bx_graph);
        //write_bin_matrix(gen_map,output_file,contig_str,binsize);
        //frac_matrix(bx_map,bx_graph,output_frac_file);
        create_contig(gen_map,contig_dict[contig_str]);
      	//contig_bx[it.first].initialize(bx_map);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
void create_contig_dict_nolinks(std::string input_bam_file,std::unordered_map<std::string,int> &contig_dict, std::unordered_map<std::string,contig_node> &contig_bx, std::string contig_str, int binsize, std::string output_file, bx_map& gen_map, full_map& chr_map ) {
        read_graph bx_graph;
        std::unordered_map<std::string,int> subset_dict;
        subset_dict = select_contigs( contig_str, contig_dict );
        init_bins_genome( chr_map, subset_dict, binsize );
        init_bins(gen_map,contig_dict[contig_str],binsize);
        connect_assembly_bam_pileup(input_bam_file,gen_map,bx_graph,binsize,chr_map,subset_dict);
        create_contig(gen_map,contig_dict[contig_str]);
}












////////////////////////////////////////////////////////////////////////////////////////////////////////
//static void link_bins( std::unordered_map<int,gen_bin> &gen_map, std::unordered_map<std::string,read_tree> &bx_graph ) {
//        int i = 0;
//        clock_t t;
//        t = clock();
//        for (auto& it : gen_map) {     //cout << it.first << endl;
//                for (int l=0; l < it.second.connected_bx.size(); l++) {
                        //cout << it.first << " bin connected " <<  it.second.connected_bx[l] << "   " << bx_graph[it.second.connected_bx[l]].connected_bins.size() << endl;
//                        for (int m=0; m < bx_graph[it.second.connected_bx[l]].connected_ints.size(); m++) {
//                                int connected_bin = bx_graph[it.second.connected_bx[l]].connected_ints[m];    //.connected_het_strings[m];
//                                if (connected_bin != it.first) { it.second.add_connection(connected_bin); }
//                        }
//                }
//                i++;
//        }
        //for (auto& it : var_dict) { it.second.count_connections(); }
//        t = clock() - t;
//        cout << "connected hashes ----- time: " << t << endl;
//        return;
//};






