#ifndef RUN_VARIANT_FILTERING_H
#define RUN_VARIANT_FILTERING_H

//////////////// c++ include //////////////////
#include <getopt.h>
#include <string>
#include <vector>
#include <map>
#include <iostream>
using namespace std;
//using namespace opt;

//////////////// linker include ////////////////
#include "read_vcf.h"
#include "read_bam.h"
#include "read_linker_output.h"
#include "write_linker_output.h"
#include "read_tree.h"
#include "variant_site.h"
#include "coord_dict.h"

////////////////////////////////
//#define filter_fraction 0.1                
#define filter_fraction 0.2                
//#define filter_fraction 0.25               
//#define min_total_cov 30
#define min_total_cov 5
#define max_total_cov 500

//////////////// definitions //////////////////
//std::unordered_map<std::string,variant_node> variant_graph;
//std::unordered_map<std::string,read_tree> read_graph;
//std::vector<vcf_entry> vcf_vector;
//std::map<std::string,int> map_str_int;

/////////////// functions /////////////////////
static void parse_filtering_options( int argc, char** argv );
void run_variant_filtering( int argc, char** argv );
void filter_het_sites(coord_dictionary& pdict, coord_dictionary& pdict2, variant_graph& vgraph, variant_graph& vgraph2);

#endif  // RUN_VARIANT_FILTERING_H

