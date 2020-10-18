#ifndef BLOCK_DICT_H
#define BLOCK_DICT_H

//////////////// c++ include //////////////////
#include <string>
#include <vector>
#include <unordered_map>
#include <set>
#include <iostream>
#include <algorithm>
using namespace std;

//////////////// linker include ////////////////
#include "map_matrix.h"

///////////////////////////////////////////////
class block_dictionary {
public:
    vector<double> subset_deltaE,subset_switchE;
    vector<int> subset_haplotype,subset_blocks;
    //bool var;  // check out what bool variant is
    int length_subset;
    map_matrix<int> subset_matrix;
    map_matrix<int> distance_matrix;
    //std::string var_base,ref_base,variant_id;
    std::unordered_map<int,std::vector<std::string>> bl_ref;  // have to add a function to add bases
    std::unordered_map<int,std::vector<std::string>> bl_alt;  // have to add a function to add bases
    std::unordered_map<int,std::vector<int>> bl_pos;  // have to add a function to add bases
    std::unordered_map<int,std::vector<int>> bl_hap;  // have to add a function to add bases
    std::unordered_map<int,int> num_hets;  // have to add a function to add bases
    std::unordered_map<int,int> block_map;
    std::unordered_map<int,int> block_map_inverted;
    std::unordered_map<int,int> bl_max;
    std::unordered_map<int,int> bl_min;
    std::unordered_map<int,int> bl_nlinks;
//    std::unordered_map<int,std::unordered_map<std::string,std::vector<link_anchor>>> bl_cnxs;  // have to add a function to add bases
    void add_bl_het(int,int,int,std::string,std::string);
    //void add_bl_connection(int,std::string,int,int,bool);
    void subset_initialize(int num_blocks, int min_hets, int min_links, map_matrix<int> block_matrix);
    void scale_matrix(int scale_limit, map_matrix<int> block_matrix);
    void subset_hap_random_initialization();
};

#endif  // BLOCK_DICT_H



