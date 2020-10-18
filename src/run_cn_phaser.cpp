#include "run_cn_phaser.h"

///////////////////////////////////////////////
namespace opt {
        static std::string input_hap_file = "./output/hap_solution_default_tenx_chr20.dat";
        static std::string input_cov_file = "./output/het_coverage_default_tenx_chr20.dat";
        static std::string id_string = "default";
	static int binsize = 10000;  // 10kb bins have 8-10 het sites on average
};

static const char* shortopts = "ho:l:m:n:b:";
static const struct option longopts[] = {
        { "hap-file",    no_argument, NULL, 'l' },
        { "cov-file",    no_argument, NULL, 'm' },
        { "id_string",   no_argument, NULL, 'n' },
        { "binsize",     no_argument, NULL, 'b' }
};

///////////////////////////////////////////////
static const char *CN_PHASE_USAGE_MESSAGE =
"Usage: mlinker cn_phase [OPTION] -l /path/to/hap_solution.dat -m /path/to/het_coverage.dat \n\n"
"\n"
"  Options\n"
"  -l,      input haplotype solution path \n"
"  -m,      input het coverage path  \n"
"  -n,      id string for output files \n"
"  -b,      binsize - raw base number (default 10000 - 10kb) \n"
"\n";

////////////////////////////////////////////////////////////////////////////////////////////////////////
static void parse_cn_phaser_options( int argc, char** argv ) {
        bool die = false;
        if (argc < 2) {
          std::cerr << "\n" << CN_PHASE_USAGE_MESSAGE;
          exit(1);
        }
        if (string(argv[1]) == "help" || string(argv[1]) == "--help") {
                std::cerr << "\n" << CN_PHASE_USAGE_MESSAGE;
                exit(1);
        }
        for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
                std::istringstream arg(optarg != NULL ? optarg : "");
                switch (c) {
                case 'h': die = true; break;
                case 'l': arg >> opt::input_hap_file; break;
                case 'm': arg >> opt::input_cov_file; break;
                case 'n': arg >> opt::id_string; break;
                case 'b': arg >> opt::binsize; break;
                }
        }
        if (die) {
          std::cerr << "\n" << CN_PHASE_USAGE_MESSAGE;
          exit(1);
        }
        cout << endl;
        cout << "############### running linker copy number phase ### " << endl;
        cout << "== input hap  === " << opt::input_hap_file << endl;
        cout << "== input cov  === " << opt::input_cov_file << endl;
        cout << "== id string  === " << opt::id_string << endl;
        cout << "== binsize    === " << opt::binsize << endl;
        cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx " << endl;
        cout << endl;
};

////////////////////////////////////////////////////
void run_cn_phaser(int argc, char** argv) {
        parse_cn_phaser_options(argc, argv);
        std::string outputFile = "./output/cn_phased_" + opt::id_string + ".dat";
        variant_graph vgraph;
        coord_dictionary pdict;
	//read_het_coverage( opt::input_cov_file, vgraph );
	//bool paired = true;
	//initialize_pdict( vgraph, pdict, paired );
	/////////////////////////////////////////
	//pdict.hap_zero_initialization();
	//read_hap_solution( opt::input_hap_file, pdict );
	/////////////////////////////////////////
	cn_map chromosome_map;
	std::vector<int> bin_array,good_bins,merged_bins,rec_bins;
	read_bin_haplotype_cn_data( opt::input_cov_file, chromosome_map, opt::binsize, bin_array );
	//initialize_copy_num_map( chromosome_map, pdict, vgraph, opt::binsize, bin_array );
	cout << " cn phasing section " << endl;
	cn_phasing( chromosome_map, bin_array, good_bins, rec_bins, merged_bins );
	write_cn_phased_bins( chromosome_map, outputFile, good_bins, merged_bins );
	//write_cn_phased( chromosome_map, outputFile, rec_bins, merged_bins );
        //cout << " cn phasing coming soon " << endl;
        return;
};

        //if ( opt::input_cov_file.substr(opt::input_cov_file.length() - 4) == ".vcf" ) {
        //   vcf_vector vvec;
        //   vvec = load_vcf_file_coverage(opt::input_cov_file,chromosome);
        //}


        //int chromosome = 0;  //2
        //if ( opt::input_cov_file.substr(opt::input_cov_file.length() - 4) == ".vcf" ) {
        //   vcf_vector vvec;
        //   vvec = load_vcf_file_coverage(opt::input_cov_file,chromosome);
        //}
