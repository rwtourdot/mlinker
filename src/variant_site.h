#ifndef VARIANT_SITE_H
#define VARIANT_SITE_H

//////////////// c++ include //////////////////
#include <string>
#include <vector>
#include <unordered_map>
#include <set>
#include <iostream>
using namespace std;

/////////////// structures ///////////////////
struct vcf_entry {
    bool bounded;
    int pos,chromosome_id;
    std::string ref_base,var_base;
};

///////////////////////////////////////////////
class variant_node {
public:
    bool var;  // check out what bool variant is
    bool paired;
    bool filter;
    int pos;
    int num_reads,num_hets;
    int total_bases;
    int unique_total;
    std::string var_base,ref_base,variant_id;
    std::vector<std::string> connected_reads_long_form,connected_hets;
    std::vector<std::string> connected_readnames;
	  std::vector<std::string> connected_bx_tags;			//Added by Greg 02/21/2020
    std::unordered_map<char,int> base_dict;  // have to add a function to add bases
    std::unordered_map<char,int> unique_dict;
    std::unordered_map<char,std::vector<std::string>> base_dict_tags;
    std::unordered_map<char,std::set<std::string>> base_dict_set;
    std::unordered_map<std::string,int> connections;
    void set_values(int,bool,std::string,std::string);
    void reset_values(int,bool,std::string,std::string,int,int,int,int,int,int);
    void add_connected_read(std::string,std::string);
    void add_connection(std::string);
    void count_connections();
    void add_base_to_dict(char,std::string);
    void unique_hash();
	void add_bx_tag(std::string);
	void get_bx(); 			//Added by Greg 02/21/2020
	void make_unique();
};

//////////////// definitions //////////////////
typedef std::unordered_map<std::string,variant_node> variant_graph;
typedef std::vector<vcf_entry> vcf_vector;

#endif  // VARIANT_SITE_H
