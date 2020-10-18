#include "run_het_coverage.h"

///////////////////////////////////////////////
namespace opt {
        static std::string chr_choice = "chr20";
        static std::string technology = "tenx";
        static std::string input_vcf_file = "./sample_data/BL1954_PCRFree.hets.recalibrated.vcf";
        static std::string input_bam_file = "./sample_data/BL1954_10xG_chr20.bam";
        static bool multiple_bams = false;
        static std::string second_input_bam_file = "./sample_data/BL1954_10xG_chr20.bam";
        static std::string id_string = "default";
};

static const char* shortopts = "ho:i:v:c:e:s:n:";
static const struct option longopts[] = {
        { "vcf-file",    no_argument, NULL, 'v' },
        { "bam-file",    no_argument, NULL, 'i' },
        { "technology",  no_argument, NULL, 'e' },
        { "chr-choice",  no_argument, NULL, 'c' },
        { "bam-file2",   no_argument, NULL, 's' },
        { "id_string",   no_argument, NULL, 'n' }
};

///////////////////////////////////////////////
static const char *COVERAGE_USAGE_MESSAGE =
"Usage: mlinker coverage [OPTION] -v /path/to/vcffile.vcf -i /path/to/bamfile.bam \n\n"
"\n"
"  Options\n"
"  -i,      input bamfile path \n"
"  -v,      input vcffile path  \n"
"  -e,      long read tech (tenx,pacbio,illumina,nanopore,hic) \n"
"  -c,      chromosome name ( chr4 ) or contig name depending on bam \n"
"  -s,      second bam file path \n"
"  -n,      id string for output files \n"
"\n";

////////////////////////////////////////////////////////////////////////////////////////////////////////
static void parse_het_coverage_options( int argc, char** argv ) {
        bool die = false;
        if (argc < 2) {
          std::cerr << "\n" << COVERAGE_USAGE_MESSAGE;
          exit(1);
        }
        if (string(argv[1]) == "help" || string(argv[1]) == "--help") {
                std::cerr << "\n" << COVERAGE_USAGE_MESSAGE;
                exit(1);
        }
        for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
                std::istringstream arg(optarg != NULL ? optarg : "");
                switch (c) {
                case 'h': die = true; break;
                case 'i': arg >> opt::input_bam_file; break;
                case 'v': arg >> opt::input_vcf_file; break;
                case 'e': arg >> opt::technology; break;
                case 'c': arg >> opt::chr_choice; break;
                case 's': opt::multiple_bams = true; arg >> opt::second_input_bam_file; break;
                case 'n': arg >> opt::id_string; break;
                }
        }
        if (die) {
          std::cerr << "\n" << COVERAGE_USAGE_MESSAGE;
          exit(1);
        }
        cout << endl;
        cout << "############### running linker het_coverage ############### " << endl;
        cout << "== chromosome === " << opt::chr_choice << endl;
        cout << "== technology === " << opt::technology << endl;
        cout << "== input bam  === " << opt::input_bam_file << endl;
        cout << "== input vcf  === " << opt::input_vcf_file << endl;
        cout << "== id string  === " << opt::id_string << endl;
        if(opt::multiple_bams) { cout << "== second bam === " << opt::second_input_bam_file << endl; }
        cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx " << endl;
        cout << endl;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
void run_het_coverage(int argc, char** argv) {
        parse_het_coverage_options(argc, argv);
        std::string coverageFile = "./output/het_coverage_" + opt::id_string + "_" + opt::technology + "_" + opt::chr_choice + ".dat";
        cout << "- output cov file: " << coverageFile << endl;
        cout << endl;
        //#################### some global structures ################
	std::map<std::string,int> chr_str_map;   //coord_dictionary pdict;
	vcf_vector vvec;
	read_graph rgraph;
        variant_graph vgraph;
        contig_name_map(opt::input_bam_file,chr_str_map);
        int chromosome = chr_str_map[opt::chr_choice];
        //#################### start of code ##########################
        vvec = load_vcf_file(opt::input_vcf_file,chromosome);
        if (opt::technology == "pacbio" || opt::technology == "tenx" || opt::technology == "illumina" || opt::technology == "nanopore" || opt::technology == "hic") {
                connect_up_variants_bam_pileup(vvec,opt::input_bam_file,chromosome,vgraph,rgraph,opt::technology);
                if(opt::multiple_bams) { connect_up_variants_bam_pileup(vvec,opt::second_input_bam_file,chromosome,vgraph,rgraph,opt::technology); }
        }
        else { cout << "error: not a valid technology choice [pacbio,tenx,illumina,nanopore,hic] " << endl; return; }
	if (opt::technology == "tenx") {
		calc_coverage_unique_hash(vgraph);
		write_het_bx_coverage(vgraph,coverageFile,opt::chr_choice);
	}
	else { write_het_coverage(vgraph,coverageFile,opt::chr_choice); }
        return;
};
