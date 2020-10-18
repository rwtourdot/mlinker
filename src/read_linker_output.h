#ifndef READ_LINKER_OUTPUT_H
#define READ_LINKER_OUTPUT_H

//////////////// c++ include //////////////////
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
using namespace std;

//////////////// linker include ////////////////
#include "coord_dict.h"
#include "variant_site.h"
#include "build_hap_matrix.h"
#include "sv_junction.h"
#include "read_tree.h"
#include "read_bam.h"
#include "hic_links.h"
#include "mbin.h"
#include "bin_assembly.h"

//////////////// definitions //////////////////
#define min_sv_reads 10

/////////////// functions /////////////////////
void read_het_coverage( std::string coverageFile, variant_graph& vgraph );
void read_hap_solution( std::string hapFile, coord_dictionary& pdict );
void read_hap_solution_initialize( std::string hapFile, coord_dictionary& pdict );
void read_hic_links_file( std::string hic_links_file, hic_link& hlink );
std::string split_string_first( std::string s, std::string delimiter, int choice);
bool is_file_exist( std::string filename );
void read_sv_file( std::string svFile, std::vector<sv_entry>& sv_list );
std::vector<vcf_entry> read_het_coverage_vvec( std::string coverageFile, int chromosome );
void read_variant_graph_file( std::string graphFile, std::string chromosome, variant_graph& vgraph, std::unordered_map<std::string,read_tree>& read_graph );
void read_scaffold( std::string scaffile, coord_dictionary& pdict );
void read_bin_haplotype_cn_data( std::string hap_cn_file, cn_map& chromosome_map, int binsize, std::vector<int>& bin_array );

#endif  // READ_LINKER_OUTPUT

