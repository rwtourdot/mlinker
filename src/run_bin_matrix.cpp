#include "run_bin_matrix.h"

/////////////// namespaces ///////////////////
namespace opt {
        static std::string chr_choice = "all";   //"chr20"
        static std::string technology = "tenx";
        static std::string input_bam_file = "./sample_data/BL1954_10xG_chr20.bam";
        static std::string id_string="default";
	static int binsize = 10000;
};

/////////////// structures ///////////////////
static const char* shortopts = "ho:b:i:c:e:n:";
static const struct option longopts[] = {
        { "bam-file",    no_argument, NULL, 'i' },
        { "technology",  no_argument, NULL, 'e' },
        { "chr-choice",  no_argument, NULL, 'c' },
        { "id_string",   no_argument, NULL, 'n' },
        { "binsize",   no_argument, NULL, 'b' }
};

static const char *MATRIX_USAGE_MESSAGE =
"Usage: mlinker matrix [OPTION] -i /path/to/bamfile.bam -b 10000 \n\n"
"\n"
"  Options\n"
"  -i,      input bamfile path \n"
"  -e,      long read tech (tenx, pacbio, nanopore, illumina, hic) \n"
"  -c,      chromosome name ( chr4 ) or contig name depending on bam \n"
"  -n,      id string for output files \n"
"  -b,      binsize - raw base number (default 10000 - 10kb) \n"
"\n";

////////////////////////////////////////////////////////////////////////////////////////////////////////
static void parse_bin_matrix_options( int argc, char** argv ) {
        bool die = false;            //cout << "there are " << argc << " arguments:\n";
        if (argc < 2) {
          std::cerr << "\n" << MATRIX_USAGE_MESSAGE;
          exit(1);
        }
        if (string(argv[1]) == "help" || string(argv[1]) == "--help") {
                std::cerr << "\n" << MATRIX_USAGE_MESSAGE;
                exit(1);
        }
        for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
                std::istringstream arg(optarg != NULL ? optarg : "");
                switch (c) {
                case 'h': die = true; break;
                case 'i': arg >> opt::input_bam_file; break;
                case 'c': arg >> opt::chr_choice; break;
		case 'e': arg >> opt::technology; break;
		case 'n': arg >> opt::id_string; break;
		case 'b': arg >> opt::binsize; break;
                }
        }
        if (die) {
          std::cerr << "\n" << MATRIX_USAGE_MESSAGE;
          exit(1);
        }
        cout << endl;
        cout << "############### running link matrix ##################### " << endl;
        cout << "== chromosome === " << opt::chr_choice << endl;
        cout << "== input bam  === " << opt::input_bam_file << endl;
        cout << "== technology === " << opt::technology << endl;
        cout << "== id string  === " << opt::id_string << endl;
        cout << "== binsize    === " << opt::binsize << endl;
        cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx " << endl;
        cout << endl;
};

////////////////////////////////////////////////////
void run_bin_matrix( int argc, char** argv ) {
        parse_bin_matrix_options(argc, argv);
        std::string bsize_string = std::to_string(opt::binsize);
        //#################### start of code ##########################
	contig_dict cdict;
	contig_bx cbx;
        bx_map gen_map;
	full_map chr_map;
        contig_name_length(opt::input_bam_file,cdict);
        std::string output_file = "./output/bxmatrix_binsize_" + std::to_string(opt::binsize) + "_" + opt::chr_choice + "_contacts.dat";
        create_contig_dict(opt::input_bam_file,cdict,cbx,opt::chr_choice,opt::binsize,output_file,gen_map,chr_map); //hg_contig
        cout << " loaded contigs " << endl;
        write_bin_matrix_genome(chr_map,output_file,opt::binsize);
        return;
};



        //write_bin_matrix(gen_map,output_file,opt::chr_choice,opt::binsize);
        //check_bx_rearrangements(contig_bx,overlap_file,rearrangement_file,hg_contig,bx_link_file,submatrix_link_file);
