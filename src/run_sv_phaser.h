#ifndef RUN_SV_CHECK_H
#define RUN_SV_CHECK_H

//////////////// c++ include //////////////////
#include <getopt.h>
#include <string>
#include <vector>
#include <iostream>
using namespace std;

//////////////// linker include ////////////////
#include "read_bam.h"
#include "read_vcf.h"
#include "read_linker_output.h"
#include "sv_junction.h"
#include "write_linker_output.h"

/////////////// functions /////////////////////
void run_sv_phaser( int argc, char** argv );
static void parse_sv_options( int argc, char** argv );

#endif  // RUN_SV_CHECK_H
