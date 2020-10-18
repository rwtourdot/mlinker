#ifndef RUN_SCAFFOLD_H
#define RUN_SCAFFOLD_H

//////////////// c++ include //////////////////
#include <getopt.h>
#include <string>
#include <vector>
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <algorithm>
using namespace std;
//using namespace opt;

//////////////// linker include ////////////////
#include "read_bam.h"
#include "read_vcf.h"
#include "read_linker_output.h"
#include "map_matrix.h"
#include "sub_matrix.h"
#include "mc_solver.h"
#include "coord_dict.h"
#include "variant_site.h"
#include "build_hap_matrix.h"
#include "write_linker_output.h"
#include "hic_links.h"
#include "block_dict.h"

//////////////// definitions //////////////////

/////////////// functions /////////////////////
//void init_block_matrix( int& num_block, map_matrix<int>& block_matrix, map_matrix<int>& total_matrix, block_dictionary& bl_dict, coord_dictionary& pdict, hic_link hlink );
void init_block_matrix( int& num_block, map_matrix<int>& block_matrix, map_matrix<int>& total_matrix, block_dictionary& bl_dict, coord_dictionary& pdict, hic_link hlink, variant_graph& hic_vgraph, int centromere_pos );
void pq_phasing(int num_block, map_matrix<int> total_matrix, block_dictionary& bl_dict, coord_dictionary& pdict, hic_link hlink, int centromere_pos );
void run_scaffold( int argc, char** argv);
static void parse_scaffold_options( int argc, char** argv );
void create_hic_links( variant_graph& hic_vgraph, hic_link& hlink );
//static void parse_region_string( std::string samtools_region );

#endif  // RUN_SCAFFOLD_H
