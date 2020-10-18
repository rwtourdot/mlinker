#ifndef RUN_EXTRACT_HASH_H
#define RUN_EXTRACT_HASH_H

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

//////////////// definitions //////////////////

/////////////// functions /////////////////////
void run_extract_hash( int argc, char** argv);
static void parse_extract_hash_options( int argc, char** argv );
static void parse_region_string( std::string samtools_region );

#endif  // RUN_EXTRACT_HASH_H
