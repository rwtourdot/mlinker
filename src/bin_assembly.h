#ifndef BIN_ASSEMBLY_H
#define BIN_ASSEMBLY_H

//////////////// c++ include //////////////////
#include <string>
#include <vector>
#include <unordered_map>
#include <set>
#include <cmath>
#include <iostream>
using namespace std;

//////////////// linker include ////////////////
#include "bin_reference.h"
#include "read_tree.h"
#include "read_bam.h"
#include "write_linker_output.h"

/////////////// definitions ///////////////////
//#define binsize 10000                   // 10 kb
#define minimum_contig 1000000          // 1 Mb

/////////////////////
class contig_node {
    std::string id;
public:
    int num_bins;
    int num_bx;
    std::unordered_map<std::string,int> num_each_bx;
    std::unordered_map<std::string,std::vector<int>> reverse_bx_map;
    std::set<std::string> cnx_bx;
    void initialize( std::unordered_map<int,gen_bin> );
};

/////////////////////
//typedef std::unordered_map<int,gen_bin> bx_map;
//typedef std::unordered_map<std::string,read_tree> bx_graph;
//typedef std::unordered_map<std::string,contig_node> contig_bx;

/////////////// functions /////////////////////
static void init_bins_genome( full_map& chr_map, std::unordered_map<std::string,int>& contig_dict, int binsize);
static std::unordered_map<std::string,int> select_contigs( std::string contig_str, std::unordered_map<std::string,int>& contig_dict );
static void init_bins( std::unordered_map<int,gen_bin> &gen_map, int contig_length );
static void link_bins( std::unordered_map<int,gen_bin> &gen_map, std::unordered_map<std::string,read_tree>& bx_graph );
static void link_bins_genome( full_map& chr_map, std::unordered_map<std::string,read_tree> &bx_graph );
void create_contig_dict( std::string input_bam_file, std::unordered_map<std::string,int> &contig_dict, std::unordered_map<std::string,contig_node> &contig_bx, std::string contig_str, int binsize, std::string output_file, bx_map& gen_map, full_map& chr_map );
static void create_contig( std::unordered_map<int,gen_bin> &gen_map, int contig_length );
void create_contig_dict_nolinks( std::string input_bam_file, std::unordered_map<std::string,int> &contig_dict, std::unordered_map<std::string,contig_node> &contig_bx, std::string contig_str, int binsize, std::string output_file, bx_map& gen_map, full_map& chr_map );

#endif  // BIN_ASSEMBLY_H
