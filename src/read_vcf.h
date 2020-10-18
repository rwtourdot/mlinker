#ifndef READ_VCF_H
#define READ_VCF_H

//////////////// c++ include //////////////////
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <map>
#include <list>
#include <algorithm>
using namespace std;

///////////// het filtering cutoffs ///
#define minimum_cov 50
#define maximum_cov 200
#define minimum_frac 0.25

/////////////// definitions ///////////////////
#define DEBUG 0
//#define DEBUG 1
#define het_cutoff 1000

//////////////// linker include ////////////////
#include "variant_site.h"

//////////////// htslib ////////////////////////
#include <htslib/vcf.h>
#include <htslib/vcfutils.h>
#include <htslib/hts.h>
#include <htslib/kstring.h>
#include <htslib/kseq.h>

extern "C" {
  #include "htslib/vcf.h"
  #include "htslib/vcfutils.h"
  #include "htslib/synced_bcf_reader.h"
}

//////////////////////////////////////////////
//std::list<std::string> allowed_base_list = { "A","T","C","G"};

/////////////// functions /////////////////////
std::vector<vcf_entry> load_and_filter_vcf_file( std::string input_vcf_file, int chromosome, std::string chr_string, std::string hetfilterFile );
std::vector<vcf_entry> load_vcf_file( std::string input_vcf_file, int chromosome );
std::vector<vcf_entry> load_vcf_file_coverage( std::string input_vcf_file, int chromosome );
std::vector<vcf_entry> load_vcf_file_total( std::string input_vcf_file );
std::map<std::string,int> load_vcf_file_header( std::string input_vcf_file, std::map<std::string,int>& chr_str_map );
void subset_het_sites( vcf_vector& vvec, int start_bound, int end_bound );

#endif  // READ_VCF_H
