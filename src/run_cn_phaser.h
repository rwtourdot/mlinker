#ifndef RUN_CN_PHASER_H
#define RUN_CN_PHASER_H

//////////////// c++ include //////////////////
#include <getopt.h>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
using namespace std;
//using namespace opt;

//////////////// linker include ////////////////
#include "read_linker_output.h"
#include "variant_site.h"
#include "coord_dict.h"
#include "bin_reference.h"
#include "copy_number_phasing.h"
#include "write_linker_output.h"
#include "read_vcf.h"

/////////////// functions /////////////////////
bool is_file_exist( std::string filename );
static void parse_cn_phaser_options( int argc, char** argv );
void run_cn_phaser(int argc, char** argv);

#endif  // RUN_CN_PHASER_H
