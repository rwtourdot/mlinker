#ifndef MBIN_H
#define MBIN_H

/////////////// c++ include //////////////////
#include <string>
#include <vector>
#include <unordered_map>
using namespace std;

////////////// linker include /////////////// 
//#include <bin_reference.h>

/////////////// structures //////////////////
class mbin {
public:
    int num_mbins;
    int bin_index;
    std::vector<int> merged_bins;
    std::vector<double> cn_switch_prob;
    double minimum_prob;
    void initialize_mbin(int);
    //void merge_bins(int,int);
    void set_energies(int,double); 
};

#endif  // MBIN_H

