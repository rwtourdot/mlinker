#ifndef BUILD_HAP_MATRIX_H
#define BUILD_HAP_MATRIX_H

//////////////// c++ include //////////////////
#include <iostream>
#include <functional>
#include <math.h>
#include <numeric>
#include <string>
#include <vector>
#include <unordered_map>
#include <set>
using namespace std;

//////////////// linker include ////////////////
#include "coord_dict.h"
#include "map_matrix.h"
#include "variant_site.h"
#include "read_tree.h"

///////////// link matrix cut ////
#define minimum_link_number 2           // {val} is the number of link hashes (reads,bx tags) which need to span two hets

/////////////// functions /////////////////////
void link_hashes( std::unordered_map<std::string,variant_node>& var_dict, std::unordered_map<std::string,read_tree>& read_graph );
void prune_graph( std::unordered_map<std::string,variant_node>& var_dict );
void calculate_link_fractions( int i, int j, map_matrix_vector& expanded, map_matrix<int>& nmatrix, map_matrix<int>& span_matrix, map_matrix<double>& corr_matrix, map_matrix<double>& diff_matrix );
void link_matrix_calculations( map_matrix_vector& expanded, map_matrix<int>& nmatrix, map_matrix<int>& span_matrix, map_matrix<double>& corr_matrix, map_matrix<double>& diff_matrix );
ptrdiff_t get_index_var( vector<int> paired_pos, int find_p );
void initialize_pdict( std::unordered_map<std::string,variant_node>& var_dict, coord_dictionary& pdict, bool loh_mode );
void initialize_solver( std::unordered_map<std::string,variant_node>& var_dict, coord_dictionary& pdict, map_matrix<double>& diff_matrix, map_matrix<int>& num_matrix_second );

void initialize_solver_loh( std::unordered_map<std::string,variant_node>& var_dict, coord_dictionary& pdict, map_matrix<double>& diff_matrix, map_matrix<int>& num_matrix_second );
void calc_coverage_unique_hash( std::unordered_map<std::string,variant_node>& var_dict );

#endif  // BUILD_HAP_MATRIX_H
