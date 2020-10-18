#ifndef MC_SOLVER_H
#define MC_SOLVER_H

//////////////// c++ include //////////////////
#include <string>
#include <vector>
#include <unordered_map>
#include <iostream>
#include <list> 
using namespace std;

//////////////// linker include ////////////////
#include "coord_dict.h"
#include "block_dict.h"
#include "map_matrix.h"
#include "sub_matrix.h"

///////////// solver_cutoffs //////
#define solver_loops 30                 // 10  // {val} the number of spin flip block flip loops
//#define solver_loops 5                 // 10  // {val} the number of spin flip block flip loops
#define solver_loops_hic 10                 // 10  // {val} the number of spin flip block flip loops
#define pos_diff_cutoff 100000          // {val} if the maximum delta genome distance - band width

/////////////// functions /////////////////////
void solver( std::unordered_map<std::string,variant_node>& var_dict, coord_dictionary& pdict,map_matrix<double> diff_matrix,map_matrix<int> num_matrix_second ); ///// fix this
static void set_switch_energy( coord_dictionary& pdict, vector<bool>& flip_maxima );
static void get_flip_positions( coord_dictionary& pdict, vector<bool>& flip_maxima );
static void energy_sum_min_max( coord_dictionary& pdict );
static void single_spin_flip_map( coord_dictionary& pdict, map_matrix<double>& diff_matrix, map_matrix<int>& nmatrix );
static void block_flip_map( coord_dictionary& pdict, map_matrix<double>& diff_matrix, map_matrix<int>& nmatrix, vector<bool>& flip_maxima );
static double energy_function_opt( submatrix_opt tempsub, vector<int> haplotype );
static double dot_product( vector<int> a, vector<double> b );
static vector<double> dot_product_matrix( vector<int> a, std::vector< std::vector<double> > b );
void call_blocks(coord_dictionary& pdict, int switch_cutoff);

/////////////////
void solver_recursive( std::unordered_map<std::string,variant_node>& var_dict, coord_dictionary& pdict, map_matrix<double> diff_matrix, map_matrix<int> num_matrix_second );
void solver_recursive_hic( block_dictionary& bdict, std::vector<int> hic_limit_loop, map_matrix<int> block_matrix );

static void single_spin_flip_recursive( coord_dictionary& pdict, map_matrix<double>& diff_matrix, map_matrix<int>& nmatrix );
static void block_flip_recursive( coord_dictionary& pdict, map_matrix<double>& diff_matrix, map_matrix<int>& nmatrix );
static void length_cutoff_nmatrix( coord_dictionary& pdict, map_matrix<double>& diff_matrix, map_matrix<int>& nmatrix, map_matrix<int>& nmatrix2 );

static void single_spin_flip_recursive_hic( block_dictionary& bdict );
static void block_flip_recursive_hic( block_dictionary& bdict );
static void energy_sum_min_max_hic( block_dictionary& bdict );
static void block_flip_brute_force_hic( block_dictionary& bdict );
//static void block_flip_brute_force2_hic( block_dictionary& bdict );

static void block_phasing_quick_flip( block_dictionary& bdict, vector<int>& graph_flip_list );
static void traverse_graph( std::unordered_map<int,std::unordered_map<int,int>> high_correlation_graph, vector<int>& traverse_list, int start_node, int subset_length );

#endif  // MC_SOLVER_H
