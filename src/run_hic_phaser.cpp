#include "run_hic_phaser.h"

/////////////// namespaces ///////////////////
namespace opt {
        static std::string chr_choice = "chr8";
        //static std::string input_vcf_file = "./sample_data/BL1954_PCRFree.hets.recalibrated.vcf";
        static std::string input_vcf_file = "/czlab/tourdot/K-562/K_562_variants_SECOND.g.vcf";
        static std::string input_het_cov = "./output/filtered_coverage_june26_K562_chr8.dat";
        //static std::string input_bam_file = "./sample_data/BL1954_10xG_chr20.bam";
        static std::string input_bam_file = "/czlab/Data/K-562/HiC/Bulk/HiC069.bam";
        static std::string id_string = "default";
        static std::string technology = "hic";
	static int start_bound = 0;
	static int end_bound = 300000000;
	static bool region_defined = false;
};

/////////////// structures ///////////////////
static const char* shortopts = "ho:i:v:c:e:s:n:m:";
static const struct option longopts[] = {
	{ "help",        no_argument, NULL, 'h' },
        { "vcf-file",    no_argument, NULL, 'v' },
        { "bam-file",    no_argument, NULL, 'i' },
        { "chr-choice",  no_argument, NULL, 'c' },
        { "id_string",   no_argument, NULL, 'n' },
	{ "cov-file",    no_argument, NULL, 'm' }
};

static const char *HIC_PHASER_USAGE_MESSAGE =
"Usage: mlinker hic_phase [OPTION] -v /path/to/vcffile.vcf -i /path/to/hic_bamfile.bam \n\n"
"\n"
"  Options\n"
"  -i,      input bamfile path \n"
"  -v,      input vcffile path  \n"
"  -m,      input het cov path  \n"
"  -c,      chromosome name ( chr4 ) or contig name depending on bam \n"
"  -n,      id string for output files \n"
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
			else { std::cerr << " region string problem  \n" << HIC_PHASER_USAGE_MESSAGE; exit(1); }
		}
		else { std::cerr << " region string problem  \n" << HIC_PHASER_USAGE_MESSAGE; exit(1); }
	}
};

///////////////////////////////////////////////////////////////////////////////////////////////////////
static void parse_hic_phaser_options( int argc, char** argv ) {
	bool die = false;
	//if(argc <= 2) { die = true; }
        if (string(argv[1]) == "help" || string(argv[1]) == "--help") {
                std::cerr << "\n" << HIC_PHASER_USAGE_MESSAGE;
                exit(1);
        }
        for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
                std::istringstream arg(optarg != NULL ? optarg : "");
                switch (c) {
		case 'h': die = true; break;
                case 'i': arg >> opt::input_bam_file; break;
                case 'v': arg >> opt::input_vcf_file; break;
                case 'm': arg >> opt::input_het_cov; break;
                case 'c': arg >> opt::chr_choice; break;
                case 'n': arg >> opt::id_string; break;
                }
        }
	//if (opt::input_bam_file.length() == 0) { die = true; }
	//if (opt::input_vcf_file.length() == 0) { die = true; }
	if (die) {
	  std::cerr << "\n" << HIC_PHASER_USAGE_MESSAGE;
	  exit(1);
	}
	parse_region_string( opt::chr_choice );
        cout << endl;
        cout << "############### running hic phaser ############### " << endl;
        cout << "== chromosome === " << opt::chr_choice << endl;
        cout << "== input bam  === " << opt::input_bam_file << endl;
        cout << "== input vcf  === " << opt::input_vcf_file << endl;
        cout << "== input cov  === " << opt::input_het_cov << endl;
        cout << "== id string  === " << opt::id_string << endl;
	if ( opt::region_defined ) {
		cout << "== start pos  === " << opt::start_bound << endl;
		cout << "== end pos    === " << opt::end_bound << endl;
	}
        cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx " << endl;
        cout << endl;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////
void run_hic_phaser( int argc, char** argv ) {
        parse_hic_phaser_options(argc, argv);
        std::string output_hic_file = "./output/hic_links_" + opt::id_string + "_" + opt::chr_choice + ".dat";
        cout << "- output hic file: " << output_hic_file << endl;
        cout << endl;
        std::map<std::string,int> chr_str_map;
        coord_dictionary pdict;
        vcf_vector vvec;
        read_graph rgraph;
        variant_graph vgraph;
        contig_name_map(opt::input_bam_file,chr_str_map);
        int chromosome = chr_str_map[opt::chr_choice];
	//##################
        read_het_coverage( opt::input_het_cov, vgraph );
        //vvec = load_vcf_file(opt::input_vcf_file,chromosome);
        initialize_pdict( vgraph, pdict, false );
        //#################### start of code ##########################
	//if (opt::region_defined) { subset_het_sites(vvec,opt::start_bound,opt::end_bound); }
        connect_up_variants_hic_bam_pileup(pdict,opt::input_bam_file,chromosome,vgraph,rgraph,opt::technology);
        link_hashes(vgraph,rgraph);
	calc_coverage_unique_hash(vgraph);
        //prune_graph(vgraph);
	write_hic_links(vgraph,output_hic_file,opt::chr_choice);
	///////////// create global datastructures for haplotype solver
	///////////// write haplotype output
        return;
}








//std::string hapsolutionFile = "./output/hic_hap_solution_" + opt::id_string + "_" + opt::chr_choice + ".dat";
//cout << "- output hap file: " << hapsolutionFile << endl;
