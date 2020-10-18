#include "variant_site.h"

void variant_node::add_base_to_dict( char base, std::string readhash_long ) { 
    base_dict[base] += 1;
    total_bases += 1;
    if (readhash_long != "nohash" ) {
    	base_dict_tags[base].push_back(readhash_long) ; 
    	base_dict_set[base].insert(readhash_long) ;
    } 
}

void variant_node::add_connected_read( std::string readhash_long, std::string readname ) { 
    connected_reads_long_form.push_back(readhash_long);
    connected_readnames.push_back(readname); 
}

//Functions added by Greg 02/21/2020
void variant_node::add_bx_tag( std::string bx_tag ) { 
    connected_bx_tags.push_back(bx_tag); 
}

void variant_node::get_bx() { 
    unique_dict['A'] = 0; unique_dict['T'] = 0;
    unique_dict['C'] = 0; unique_dict['G'] = 0;
    unique_dict['D'] = 0; unique_dict['I'] = 0;
    unique_total = 0;
    for (auto& it2 : base_dict_tags) { 
	char base = it2.first;
	std::vector<std::string> bx_tags_base = base_dict_tags[base];		//base_dict_tags[base] returns the read hashes for reads supporting that base
	std::set<std::string> unique_tags(bx_tags_base.begin(),bx_tags_base.end());			//make a set of read hashes for reads supporting that base -- removes duplicate reads. Use this for bx_tags
		for (auto elem : unique_tags) {						//loop added by Greg 02/21/2020
					add_bx_tag(elem);
		}
	unique_dict[base] = unique_tags.size();								//get the number of unique reads supporting that base
	//cout << pos << "  " << base << "  " << unique_tags.size() << endl;
	//std::set<std::string>::iterator it;
	//for (it = unique_tags.begin(); it != unique_tags.end(); it++) { cout << *it << endl; }
    }
    for (auto& it2 : unique_dict) { unique_total += it2.second; }		//go through each base in unique_dict and add the number of reads supporting that base to unique_total
	connected_readnames = connected_bx_tags;
	connected_reads_long_form = connected_bx_tags;
}

void variant_node::unique_hash() { 
    unique_dict['A'] = 0; unique_dict['T'] = 0;
    unique_dict['C'] = 0; unique_dict['G'] = 0;
    unique_dict['D'] = 0; unique_dict['I'] = 0;
    unique_total = 0;
    for (auto& it2 : base_dict_tags) { 
	char base = it2.first;
	std::vector<std::string> bx_tags_base = base_dict_tags[base];		//base_dict_tags[base] returns the read hashes for reads supporting that base
	std::set<std::string> unique_tags(bx_tags_base.begin(),bx_tags_base.end());			//make a set of read hashes for reads supporting that base -- removes duplicate reads. Use this for bx_tags
	unique_dict[base] = unique_tags.size();								//get the number of unique reads supporting that base
	//cout << pos << "  " << base << "  " << unique_tags.size() << endl;
	//std::set<std::string>::iterator it;
	//for (it = unique_tags.begin(); it != unique_tags.end(); it++) { cout << *it << endl; }
    }
    for (auto& it2 : unique_dict) { unique_total += it2.second; }		//go through each base in unique_dict and add the number of reads supporting that base to unique_total
}

//Here, defines the number of connected variants (num_hets) and the number of reads containing the variant (num_reads)
//connections is defined in variant_site.h [std::unordered_mat<std::string,int> connections;] 
void variant_node::count_connections() {
    num_hets = connections.size();
	num_reads = connected_reads_long_form.size();
    //num_reads = connected_bx_tags.size();
}

void variant_node::make_unique() {
		std::set<std::string> unique_tags_long(connected_reads_long_form.begin(), connected_reads_long_form.end());
		std::vector<std::string> unique_tags_long_vec(unique_tags_long.begin(), unique_tags_long.end());
		connected_reads_long_form = unique_tags_long_vec;
		std::set<std::string> unique_tags(connected_readnames.begin(), connected_readnames.end());
		std::vector<std::string> unique_tags_vec(unique_tags.begin(), unique_tags.end());
		connected_readnames = unique_tags_vec;
}

void variant_node::add_connection( std::string hashname ) { connections[hashname] += 1; };




void variant_node::reset_values( int position, bool variant, std::string variant_base, std::string reference_base, int abase, int tbase, int cbase, int gbase, int dbase, int ibase ) {
    pos = position;                 // int
    var = variant;                  // bool
    var_base = variant_base;        // std::string
    ref_base = reference_base;      // std::string
    base_dict['A'] = abase;
    base_dict['T'] = tbase;
    base_dict['C'] = cbase;
    base_dict['G'] = gbase;
    base_dict['D'] = dbase;
    base_dict['I'] = ibase;
    total_bases = abase + tbase + cbase + gbase + dbase + ibase;
    filter = false;
    if (variant == true) { variant_id = std::to_string(position) + "_" + reference_base + "_" + variant_base; }
    else {                 variant_id = std::to_string(position) + "_" + reference_base + "_" + reference_base; }
}

void variant_node::set_values( int position, bool variant, std::string variant_base, std::string reference_base ) {
    pos = position;                 // int
    var = variant;                  // bool
    var_base = variant_base;        // std::string
    ref_base = reference_base;      // std::string
    connected_reads_long_form.reserve(100);
    base_dict['A'] = 0; 	
    base_dict['T'] = 0; 
    base_dict['C'] = 0; 
    base_dict['G'] = 0; 
    base_dict['D'] = 0; 
    base_dict['I'] = 0;
    total_bases = 0;
    filter = false;
    if (variant == true) { variant_id = std::to_string(position) + "_" + reference_base + "_" + variant_base; }
    else {                 variant_id = std::to_string(position) + "_" + reference_base + "_" + reference_base; }
}

	//for (int i=0; i < bx_tags_base.size(); i++) { cout << bx_tags_base[i] << " " << endl; }
