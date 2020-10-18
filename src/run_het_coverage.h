#ifndef RUN_HET_COVERAGE_H
#define RUN_HET_COVERAGE_H

//////////////// c++ include //////////////////
#include <getopt.h>
#include <string>
#include <vector>
#include <iostream>
using namespace std;
//using namespace opt;

//////////////// linker include ////////////////
#include "read_bam.h"
#include "read_vcf.h"
#include "write_linker_output.h"
#include "read_tree.h"
#include "variant_site.h"
#include "build_hap_matrix.h"

//////////////// definitions //////////////////
//std::unordered_map<std::string,variant_node> variant_graph;
//std::unordered_map<std::string,read_tree> read_graph;
//std::vector<vcf_entry> vcf_vector;
//std::map<std::string,int> map_str_int;

/////////////// functions /////////////////////
static void parse_het_coverage_options( int argc, char** argv );
void run_het_coverage( int argc, char** argv );

#endif  // RUN_HET_COVERAGE_H
