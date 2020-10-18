#ifndef RUN_SOLVE_HAP_H
#define RUN_SOLVE_HAP_H

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
void run_solve_hap( int argc, char** argv);
static void parse_solve_hap_options( int argc, char** argv );
//static void parse_region_string( std::string samtools_region );

#endif  // RUN_SOLVE_HAP_H
