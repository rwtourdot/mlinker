#ifndef HIC_LINKS_H
#define HIC_LINKS_H

//////////////// c++ include //////////////////
#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <iostream>
using namespace std;

//////////////// linker include ////////////////
#include "map_matrix.h"
#include "variant_site.h"
//#include "read_linker_output.h"

///////////////////// dictionary with general info
class hic_link {
public:
    vector<int> anchor_pos1,anchor_pos2,link_depth;
    vector<std::string> anchor_var1,anchor_var2;
    vector<std::string> anchor_ref_base1,anchor_ref_base2;
    vector<std::string> anchor_read_base1,anchor_read_base2;
    vector<bool> anchor_var_bool1,anchor_var_bool2;
    std::unordered_map<int,std::vector<int>> hic_pos_map;
    std::unordered_map<std::string,std::vector<std::string>> hic_link_map;
    std::unordered_map<std::string,std::vector<int>> hic_link_depth_map;
    //std::unordered_map<int,std::vector<std::string>> hic_pos_var;
    int nlinks,nc;
    void add_link(int,int,std::string,std::string,int);
};

std::string split_string_first_hic_link( std::string, std::string, int );

//////////////// definitions //////////////////
//std::unordered_map<std::string,int> contig_dict;

#endif  // HIC_LINKS_H
