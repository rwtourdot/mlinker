#ifndef BIN_REFERENCE_H
#define BIN_REFERENCE_H

//////////////// c++ include //////////////////
#include <string>
#include <vector>
#include <unordered_map>
#include <set>
#include <iostream>
#include <sstream>
using namespace std;

/////////////// structures ///////////////////
class gen_bin {
public:
    int pos;
    int maxb;
    int num_unique_bx,num_bx;
    std::vector<std::string> unique_bx;
    std::vector<std::string> connected_bx;     // connected_hets;
    std::unordered_map<int,int> connections;
    std::unordered_map<std::string,std::unordered_map<int,int>> chr_connections;
    void set_values(int);
    void add_connected_read(std::string);
    void add_connection(std::string,int);
    void find_unique_bx();
};

class cn_bin {
public:
    int start_pos,end_pos,avg_pos;
    int num_hets;
    int num_blocks;
    int switch_spin;
    int hapA_tot,hapB_tot;
    double hapA_avg,hapB_avg;
    std::vector<std::string> ref_base,var_base;
    std::vector<int> ref_cov,var_cov;
    std::vector<int> het_pos,het_hap;
    std::vector<int> het_block;
    void initialize(int);
    void add_het(int,int,int,std::string,std::string,int,int);
    void cov_calc();
    void flip_hap();
    void enter_coverage(double,double,int);
};

//////////////// structures //////////////////
typedef std::unordered_map<int,gen_bin> bx_map;
typedef std::unordered_map<std::string,bx_map> full_map;
typedef std::unordered_map<int,cn_bin> cn_map;

#endif  // BIN_REFERENCE_H
