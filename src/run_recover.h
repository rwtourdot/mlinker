#ifndef RUN_RECOVER_H
#define RUN_RECOVER_H

//////////////// c++ include //////////////////
#include <getopt.h>
#include <string>
#include <vector>
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <algorithm>
using namespace std;

//////////////// linker include ////////////////
#include "read_linker_output.h"
#include "coord_dict.h"
#include "variant_site.h"
#include "build_hap_matrix.h"
#include "write_linker_output.h"
#include "recovered_site.h"

/////////////// functions /////////////////////
void run_recover( int argc, char** argv);
static void parse_recover_options( int argc, char** argv );
void recover_links( std::string chr_choice, coord_dictionary& pdict, variant_graph& vgraph, std::unordered_map<std::string,read_tree>& rgraph, std::string technology, coord_dictionary& rdict, std::string scaffoldsolutionFile, std::map<int,recovered_node>& recovered_map );

#endif  // RUN_RECOVER_H
