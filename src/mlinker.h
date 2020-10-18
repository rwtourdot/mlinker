#ifndef MLINKER_H
#define MLINKER_H

//////////////// c++ include //////////////////
#include <limits>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <unordered_map>
#include <unordered_set>
#include <typeinfo>
#include <algorithm>
#include <iterator>
#include <set>
#include <getopt.h>
using namespace std;

//////////////// bamtools //////////////////////
#include "api/BamWriter.h"
#include "api/BamReader.h"
#include "api/BamAlignment.h"
#include "api/BamAlgorithms.h"
#include "api/BamMultiReader.h"
#include "api/BamAux.h"

//////////////// htslib ////////////////////////
#include <htslib/vcf.h>
#include <htslib/hts.h>

extern "C" {
  #include "htslib/synced_bcf_reader.h"
}

//////////////// linker include ////////////////
//#include "run_link_phaser.h"
#include "run_extract_hash.h"
#include "run_solve_hap.h"
#include "run_het_coverage.h"
#include "run_cn_phaser.h"
#include "run_bin_matrix.h"
#include "run_sv_phaser.h"
#include "run_variant_filtering.h"
#include "run_hic_phaser.h"
#include "run_scaffold.h"
#include "run_recover.h"

//#include "run_bx_bin.h"
/////////////// definitions ///////////////////
////#define DEBUG 1
//#define DEBUG 0
//#define het_cutoff 1000

////////////// read cutoffs ///////
//#define minimum_mapq 20                 // read map quality greater than {val}
//#define minimum_baseq 8                 // read base quality greater than {val}

///////////// solver_cutoffs //////
//#define solver_loops 10                 // {val} the number of spin flip block flip loops
//#define switch_cut 150                        // -{val} is the minimum switchE required to flip block
//#define pos_diff_cutoff 100000          // {val} if the maximum delta genome distance - band width

///////////// com line input /////
#define no_argument 0
#define required_argument 1
#define optional_argument 2

///////////// link matrix cut ////
//#define minimum_link_number 2           // {val} is the number of link hashes (reads,bx tags) which need to span two hets

///////////// het filtering cutoffs ///
//#define minimum_cov 50
//#define maximum_cov 200
//#define minimum_frac 0.25

//////////////// definitions //////////////////
typedef std::map<std::string,int> map_str_int;
typedef std::unordered_map<std::string,int> umap_str_int;

/////////////// functions /////////////////////
int main(int argc, char** argv);

#endif  // MLINKER_H


/////////////// namespaces ///////////////////
//namespace opt {
//        static std::string chr_choice = "chr20";
//        static std::string technology = "pacbio";
//        static std::string input_vcf_file = "/czlab/Data/HCC1954_HaplotypeCalling/v.0/BL1954_PCRFree.hets.recalibrated.vcf";
        //static std::string input_vcf_file = "/czlab/Data/HCC1954_HaplotypeCalling/BL1954_PCRFree.hets.recalibrated.vcf";
        //static std::string input_vcf_file = "/czlab/Data/HCC1954_HaplotypeCalling/v.2/HCC1954BL_hets/HCC1954BL.hets.recalibrated.vcf.gz";
//        static std::string input_bam_file = "/czlab/Data/HCC1954_PacBio/BL1954.pacbio_rawReads.GRCh38.bam";
//        static bool multiple_bams = false;
//        static std::string second_input_bam_file = "/czlab/Data/HCC1954_10X_New/GRCh38/BL1954_10xG.bam";
//        static std::string id_string="default";
//};

/////////////// structures ///////////////////
//static const char* shortopts = "o:i:v:c:e:s:p:n:";
//static const struct option longopts[] = {
//        { "vcf-file",    no_argument, NULL, 'v' },
//        { "bam-file",    no_argument, NULL, 'i' },
//        { "technology",  no_argument, NULL, 'e' },
//        { "chr-choice",  no_argument, NULL, 'c' },
//        { "bam-file2",   no_argument, NULL, 's' },
//        { "paired_tm",   no_argument, NULL, 'p' },
//        { "id_string",   no_argument, NULL, 'n' }
//}
