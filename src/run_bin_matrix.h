#ifndef RUN_BIN_MATRIX_H
#define RUN_BIN_MATRIX_H

//////////////// c++ include ///////////////////
#include <string>
#include <vector>
#include <unordered_map>
#include <getopt.h>
#include <iostream>
using namespace std;
//using namespace opt;

//////////////// linker include ////////////////
#include "bin_assembly.h"
#include "bin_reference.h"
#include "read_bam.h"

/////////////// structures /////////////////////
typedef std::unordered_map<std::string,int> contig_dict;
typedef std::unordered_map<std::string,contig_node> contig_bx;
//std::unordered_map<std::string,hg_contig_map> hg_contig;

/////////////// functions //////////////////////
static void parse_bin_matrix_options( int argc, char** argv );
void run_bin_matrix( int argc, char** argv );

#endif  // RUN_BIN_MATRIX_H
