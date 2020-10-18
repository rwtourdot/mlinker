#ifndef READ_BAM_H
#define READ_BAM_H

//////////////// c++ include //////////////////
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <string>
#include <map>
#include <unordered_map>
#include <algorithm>
#include <cmath>
using namespace std;

//////////////// linker include ////////////////
#include "variant_site.h"
#include "read_tree.h"
#include "bin_reference.h"
#include "sv_junction.h"
#include "coord_dict.h"

//////////////// bamtools //////////////////////
#include "api/BamWriter.h"
#include "api/BamReader.h"
#include "api/BamAlignment.h"
#include "api/BamAlgorithms.h"
#include "api/BamMultiReader.h"
#include "api/BamAux.h"

/////////////// definitions ///////////////////
#define samsum_lt 2000
#define minimum_mapq 20                 // read map quality greater than {val}
#define minimum_mapq_hic 20                 // read map quality greater than {val}
#define minimum_baseq 8                 // read base quality greater than {val}
#define minimum_baseq_hic 20                 // read base quality greater than {val}
#define minimum_baseq_tenx 20                 // read base quality greater than {val}

//#define binsize 10000                   // 10 kb
#define minimum_contig 1000000          // 1 Mb
#define sv_buffer 1000                   // 100
#define sv_bx_buffer 100000

//////////////// structures //////////////////

/////////////// functions /////////////////////
std::vector<std::string> split_string( std::string teststring );
std::string return_index_string( std::string );
std::string return_read_id( std::string );
std::string technology_hash( std::string tech, BamTools::BamAlignment al );
std::pair<char,int> get_base_readpos( std::vector<BamTools::CigarOp> &cigar_str, std::string aligned_bases, int variant_pos, int start_pos, std::string qual );
bool alignment_boolean_check( BamTools::BamAlignment al );
bool alignment_boolean_check_hic( BamTools::BamAlignment al );
bool alignment_boolean_check_lenient( BamTools::BamAlignment al );
bool get_tag_value( const std::string &tag, variant_graph &variant_cnx );
//void connect_up_variants( vcf_vector vvec, std::unordered_map<std::string,BamTools::BamAlignment*> bam_map );
void connect_up_variants_bam_pileup( vcf_vector vvec, std::string inputFilename,int chromosome, variant_graph &variant_cnx, read_graph &rgraph, std::string technology );
void connect_up_variants_hic_bam_pileup( coord_dictionary& pdict, std::string inputFilename,int chromosome, variant_graph &variant_cnx, read_graph &rgraph, std::string technology );
void connect_assembly_bam_pileup( std::string inputFilename, bx_map &gen_map, read_graph &bx_graph, int binsize, full_map& chr_map, std::unordered_map<std::string,int>& contig_dict );
//void connect_sv_bam_pileup( std::string inputFilename, read_graph &bx_graph, std::vector<sv_entry>& sv_list, std::map<std::string,int> chr_str_map );
void connect_sv_read_pileup( std::string inputFilename, read_graph &bx_graph, std::vector<sv_entry>& sv_list, std::map<std::string,int> chr_str_map );
void connect_sv_hets( vcf_vector vvec, std::string inputFilename, read_graph &bx_graph, std::vector<sv_entry>& sv_list, std::map<std::string,int> chr_str_map );
void contig_name_map( std::string inputFilename, std::map<std::string,int> &chr_str_map);
void contig_name_length( std::string inputFilename, std::unordered_map<std::string,int> &contig_dict );
static void open_bam_file( BamTools::BamReader &reader, std::string inputFilename );
void sv_phase_het( BamTools::BamReader& reader, BamTools::BamAlignment& al, sv_entry& sv_object, int start, int end, int chromosome, vcf_entry vnode, int& possible_hets, bool side );

#endif  // READ_BAM_H
