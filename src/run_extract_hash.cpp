#include "run_extract_hash.h"

/////////////// namespaces ///////////////////
namespace opt {
        static std::string chr_choice = "chr20";
        static std::string technology = "tenx";
        static std::string input_vcf_file = "./sample_data/BL1954_PCRFree.hets.recalibrated.vcf";
	static std::string input_cov_file = "./output/het_coverage_default_tenx_chr20.dat";
	static bool cov_load = false;
	static bool vcf_load = false;
        static std::string input_bam_file = "./sample_data/BL1954_10xG_chr20.bam";
        static bool multiple_bams = false;
        static std::string second_input_bam_file = "./sample_data/BL1954_10xG_chr20.bam";
        static std::string id_string = "default";
	static int start_bound = 0;
	static int end_bound = 300000000;
	static bool region_defined = false;
  static bool output_graph_hash = false;
};

/////////////// structures ///////////////////
static const char* shortopts = "hgo:i:v:m:c:e:s:n:";
static const struct option longopts[] = {
	{ "help",        no_argument, NULL, 'h' },
        { "vcf-file",    no_argument, NULL, 'v' },
        { "cov-file",    no_argument, NULL, 'm' },
        { "bam-file",    no_argument, NULL, 'i' },
        { "technology",  no_argument, NULL, 'e' },
        { "chr-choice",  no_argument, NULL, 'c' },
        { "bam-file2",   no_argument, NULL, 's' },
        { "id_string",   no_argument, NULL, 'n' },
        { "output_graph_hash",   no_argument, NULL, 'g' }
};

///////////////////////////////////////////////////////////////////////////////////////////////////////
static const char *EXTRACT_USAGE_MESSAGE =
"Usage: mlinker extract [OPTION] -v /path/to/vcffile.vcf -i /path/to/bamfile.bam \n\n"
"\n"
"  Options\n"
"  -i,      input bamfile path \n"
"  -v,      input vcffile path  \n"
"  -m,      input coverage file path ( diploid ) \n"
"  -e,      long read tech ( tenx, pacbio, nanopore, illumina, hic ) \n"
"  -c,      chromosome name ( chr4 ) or contig name depending on bam \n"
"  -s,      second bam file path \n"
"  -n,      id string for output files \n"
"  -g,      flag to output graph_hash file (default no) \n"
"\n";

///////////////////////////////////////////////////////////////////////////////////////////////////////
static void parse_region_string( std::string samtools_region ) {
	std::string chr_delimiter=":";
	std::string pos_delimiter="-";
	std::size_t index = samtools_region.find(chr_delimiter);
	if ( index != std::string::npos ) {
		std::string cid = samtools_region.substr(0,index);
		std::string gpos = samtools_region.substr(index+1);
		gpos.erase(std::remove(gpos.begin(),gpos.end(),','),gpos.end());
		std::size_t sub_index = gpos.find(pos_delimiter);
        	if ( sub_index != std::string::npos ) {
			opt::chr_choice = cid;
                	std::string start_string = gpos.substr(0,sub_index);
                	std::string end_string = gpos.substr(sub_index+1);
			int start_int = std::stoi(start_string);
			int end_int = std::stoi(end_string);
			if ( end_int > start_int ) {
				opt::chr_choice = cid;
				opt::start_bound = start_int;
				opt::end_bound = end_int;
				opt::region_defined = true;
			}
			else { std::cerr << " region string problem  \n" << EXTRACT_USAGE_MESSAGE; exit(1); }
		}
		else { std::cerr << " region string problem  \n" << EXTRACT_USAGE_MESSAGE; exit(1); }
	}
};

///////////////////////////////////////////////////////////////////////////////////////////////////////
static void parse_extract_hash_options( int argc, char** argv ) {
	bool die = false; //bool vcf_load = false; //bool cov_load = false;
	//if(argc <= 2) { die = true; }
  if (argc < 2) {
    std::cerr << "\n" << EXTRACT_USAGE_MESSAGE;
    exit(1);
  }
	if (string(argv[1]) == "help" || string(argv[1]) == "--help") {
		std::cerr << "\n" << EXTRACT_USAGE_MESSAGE;
		exit(1);
	}
        for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
                std::istringstream arg(optarg != NULL ? optarg : "");
                switch (c) {
		case 'h': die = true; break;
                case 'i': arg >> opt::input_bam_file; break;
                case 'v': opt::vcf_load = true; arg >> opt::input_vcf_file; break;
                case 'm': opt::cov_load = true; arg >> opt::input_cov_file; break;
                case 'e': arg >> opt::technology; break;
                case 'c': arg >> opt::chr_choice; break;
                case 's': opt::multiple_bams = true; arg >> opt::second_input_bam_file; break;
                case 'n': arg >> opt::id_string; break;
                case 'g': opt::output_graph_hash = true; break;
                }
        }
	//if (opt::input_bam_file.length() == 0) { die = true; }
	//if (opt::input_vcf_file.length() == 0) { die = true; }
	if (!opt::vcf_load && !opt::cov_load) { die = true; }
	if (die) {
	  std::cerr << "\n" << EXTRACT_USAGE_MESSAGE;
	  exit(1);
	}
	parse_region_string( opt::chr_choice );
        cout << endl;
        cout << "############### running link phaser ############### " << endl;
        cout << "== chromosome === " << opt::chr_choice << endl;
        cout << "== technology === " << opt::technology << endl;
        cout << "== input bam  === " << opt::input_bam_file << endl;
        cout << "== output graph hash  === " << opt::output_graph_hash << endl;
        if(opt::vcf_load) { cout << "== input vcf  === " << opt::input_vcf_file << endl; }
        if(opt::cov_load) { cout << "== input cov  === " << opt::input_cov_file << endl; }
        cout << "== id string  === " << opt::id_string << endl;
        if(opt::multiple_bams) { cout << "== second bam === " << opt::second_input_bam_file << endl; }
	if ( opt::region_defined ) {
		cout << "== start pos  === " << opt::start_bound << endl;
		cout << "== end pos    === " << opt::end_bound << endl;
	}
        cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx " << endl;
        cout << endl;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////
void run_extract_hash( int argc, char** argv ) {
        parse_extract_hash_options(argc, argv);
        std::string output_link_file = "./output/graph_variant_" + opt::id_string + "_" + opt::technology + "_" + opt::chr_choice + ".dat";
        std::string output_hash_file = "./output/graph_hash_" + opt::id_string + "_" + opt::technology + "_" + opt::chr_choice + ".dat";
        //std::string hapsolutionFile = "./output/hap_solution_" + opt::id_string  + "_" + opt::technology + "_" + opt::chr_choice + ".dat";
        cout << "- output variant link file: " << output_link_file << endl;
        if (opt::output_graph_hash) { cout << "- output hash link file: " << output_hash_file << endl; }
        cout << endl;
        std::map<std::string,int> chr_str_map;
        coord_dictionary pdict;
        vcf_vector vvec;
        read_graph rgraph;
        variant_graph vgraph;
        contig_name_map(opt::input_bam_file,chr_str_map);
        int chromosome = chr_str_map[opt::chr_choice];
        //#################### start of code ##########################
	//if (opt::filter) { };  // vvec = read_het_coverage( opt::input_cov_file, vgraph );
        //vvec = load_and_filter_vcf_file(opt::input_vcf_file,chromosome,opt::chr_choice,hetfilterFile);
	if (opt::cov_load) { vvec = read_het_coverage_vvec( opt::input_cov_file, chromosome ); }
        else { vvec = load_vcf_file( opt::input_vcf_file, chromosome ); }
	if (opt::region_defined) {
		subset_het_sites( vvec, opt::start_bound, opt::end_bound );
	}
        if (opt::technology == "pacbio" || opt::technology == "tenx" || opt::technology == "illumina" || opt::technology == "nanopore" || opt::technology == "hic") {
                connect_up_variants_bam_pileup(vvec,opt::input_bam_file,chromosome,vgraph,rgraph,opt::technology);
                if(opt::multiple_bams) { connect_up_variants_bam_pileup(vvec,opt::second_input_bam_file,chromosome,vgraph,rgraph,opt::technology); }
        }
        else { cout << "error: not a valid technology choice [pacbio,tenx,illumina,nanopore,hic] " << endl; return; }
        write_variant_link_list(vgraph,rgraph,output_link_file,opt::chr_choice);
        if (opt::output_graph_hash) { write_hash_link_list(vgraph,rgraph,output_hash_file,opt::chr_choice); }
	////////
        return;
};
