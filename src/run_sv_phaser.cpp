#include "run_sv_phaser.h"

///////////////////////////////////////////////
namespace opt {
        static std::string technology = "tenx";
        static std::string input_vcf_file = "./sample_data/BL1954_PCRFree.hets.recalibrated.vcf";
        static std::string input_bam_file = "./sample_data/BL1954_10xG_chr20.bam";
        static std::string input_sv_file = "./sample_data/Interchromosomal.txt";
        static std::string id_string = "default";
};

static const char* shortopts = "ho:i:v:e:u:n:";
static const struct option longopts[] = {
        { "help",        no_argument, NULL, 'h' },
        { "vcf-file",    no_argument, NULL, 'v' },
        { "bam-file",    no_argument, NULL, 'i' },
        { "technology",  no_argument, NULL, 'e' },
        { "sv_file",     no_argument, NULL, 'u' },
        { "id_string",   no_argument, NULL, 'n' }
};

///////////////////////////////////////////////
static const char *SV_USAGE_MESSAGE =
"Usage: mlinker sv_phase [OPTION] -v /path/to/vcffile.vcf -i /path/to/bamfile.bam -s /path/to/svfile.txt \n\n"
"\n"
"  Options\n"
"  -i,      input bamfile path \n"
"  -v,      input vcffile path  \n"
"  -e,      long read tech (tenx,pacbio,nanopore) \n"
"  -u,      sv file path \n"
"  -n,      id string for output files \n"
"\n";

////////////////////////////////////////////////////////////////////////////////////////////////////////
static void parse_sv_options( int argc, char** argv ) {
        bool die = false;
        if (argc < 2) {
          std::cerr << "\n" << SV_USAGE_MESSAGE;
          exit(1);
        }
        if (string(argv[1]) == "help" || string(argv[1]) == "--help") {
                std::cerr << "\n" << SV_USAGE_MESSAGE;
                exit(1);
        }
        for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
                std::istringstream arg(optarg != NULL ? optarg : "");
                switch (c) {
                case 'h': die = true; break;
                case 'i': arg >> opt::input_bam_file; break;
                case 'v': arg >> opt::input_vcf_file; break;
                case 'e': arg >> opt::technology; break;
                case 'u': arg >> opt::input_sv_file; break;
                case 'n': arg >> opt::id_string; break;
                }
        }
        if (die) {
          std::cerr << "\n" << SV_USAGE_MESSAGE;
          exit(1);
        }
        cout << endl;
        cout << "############### running sv check ################### " << endl;
        cout << "== technology === " << opt::technology << endl;
        cout << "== input bam  === " << opt::input_bam_file << endl;
        cout << "== input vcf  === " << opt::input_vcf_file << endl;
        cout << "== input sv   === " << opt::input_sv_file << endl;
        cout << "== id string  === " << opt::id_string << endl;
        cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx " << endl;
        cout << endl;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
void run_sv_phaser(int argc, char** argv) {
	parse_sv_options(argc, argv);
	std::string output_phased_sv_file = "./output/sv_phased_" + opt::id_string + ".dat";
	std::vector<sv_entry> sv_list;
	read_graph rgraph;
	read_sv_file( opt::input_sv_file, sv_list );
	std::map<std::string,int> chr_str_map;
	cout << "loaded sv file" << endl;
	contig_name_map(opt::input_bam_file,chr_str_map);         //connect_sv_bam_pileup( opt::input_bam_file, rgraph, sv_list, chr_str_map );
	connect_sv_read_pileup( opt::input_bam_file, rgraph, sv_list, chr_str_map );
	cout << "loaded sv reads" << endl;
        vcf_vector vvec;
	vvec = load_vcf_file_total(opt::input_vcf_file);
	connect_sv_hets( vvec, opt::input_bam_file, rgraph, sv_list, chr_str_map );
	write_phased_sv( sv_list, output_phased_sv_file );
        return;
}








        //for (int i = 0; i < sv_list.size(); i++) {
        //      cout << sv_list[i].chr_one << ":" << sv_list[i].pos_one << "  " << sv_list[i].chr_two << ":" << sv_list[i].pos_two << endl;
        //      for (auto& it : sv_list[i].hap_end1) {
        //              cout << "0 " << sv_list[i].chr_one << " " << it.first << "  " << it.second << "  " << sv_list[i].num_end1[it.first] << endl;
        //      }
        //      for (auto& it : sv_list[i].hap_end2) {
        //              cout << "1 " << sv_list[i].chr_two << " " << it.first << "  " << it.second << "  " << sv_list[i].num_end2[it.first] << endl;
        //      }
        //      cout << "######" << endl;
        //}
