#include "run_variant_filtering.h"

///////////////////////////////////////////////
namespace opt {
	static std::string tumor_cov_file = "./output/het_coverage_default_tenx_chr20.dat";
	static std::string normal_cov_file = "./output/het_coverage_default_tenx_chr20.dat";
        static std::string id_string = "default";
        static bool filter_tumor = false;
};

static const char* shortopts = "ho:t:m:n:";
static const struct option longopts[] = {
        { "normal-cov-file",    no_argument, NULL, 'm' },
        { "tumor-cov-file",   no_argument, NULL, 't' },
        { "id_string",   no_argument, NULL, 'n' }
};

///////////////////////////////////////////////
static const char *FILTER_USAGE_MESSAGE =
"Usage: mlinker filter [OPTION] -t /path/to/tumor_coverage_file.dat -m /path/to/normal_coverage_file.dat \n\n"
"\n"
"  Options\n"
"  -t,      tumor het coverage path \n"
"  -m,      normal het coverage path \n"
"  -n,      id string for output files \n"
"\n";

////////////////////////////////////////////////////////////////////////////////////////////////////////
static void parse_filtering_options( int argc, char** argv ) {
        bool die = false;
				if (argc < 2) {
          std::cerr << "\n" << FILTER_USAGE_MESSAGE;
          exit(1);
        }
        if (string(argv[1]) == "help" || string(argv[1]) == "--help") {
                std::cerr << "\n" << FILTER_USAGE_MESSAGE;
                exit(1);
        }
        for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
                std::istringstream arg(optarg != NULL ? optarg : "");
                switch (c) {
                case 'h': die = true; break;
                case 't': opt::filter_tumor = true; arg >> opt::tumor_cov_file; break;
                case 'm': arg >> opt::normal_cov_file; break;
                case 'n': arg >> opt::id_string; break;
                }
        }
        if (die) {
        	std::cerr << "\n" << FILTER_USAGE_MESSAGE;
        	exit(1);
        }
        cout << endl;
        cout << "############### running linker filter ############### " << endl;
        if(opt::filter_tumor) { cout << "== tumor  input cov  === " << opt::tumor_cov_file << endl; }
        cout << "== normal input cov  === " << opt::normal_cov_file << endl;
        cout << "== id string         === " << opt::id_string << endl;
        cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx " << endl;
        cout << endl;
};

////////////////////////////////////////////////////
void filter_het_sites(coord_dictionary& pdict, coord_dictionary& pdict2, variant_graph& vgraph, variant_graph& vgraph2) {
	int number_filter = 0;
        for (int j = 0; j < pdict.num_paired; j++) {
                int loop_pos = pdict.sorted_paired_positions[j];
                std::string het_string = pdict.paired_dict[loop_pos][0];
		std::string het_string2 = pdict.paired_dict[loop_pos][1];
		std::string het_string_variant;
		cout << loop_pos << " " << het_string << " " << het_string2 << endl;
		//if ( vgraph.find(het_string) != vgraph.end() && vgraph2.find(het_string) != vgraph2.end() ) {
		//if ( vgraph.find(het_string_variant) != vgraph.end() && vgraph2.find(het_string_variant) != vgraph2.end() ) {
		if ( vgraph2[het_string].var ) { het_string_variant = het_string; }
		else {                           het_string_variant = het_string2; }
                std::string rbase = vgraph[het_string].ref_base;        std::string rbase2 = vgraph2[het_string].ref_base;
                std::string vbase = vgraph[het_string].var_base;        std::string vbase2 = vgraph2[het_string].var_base;
                char rbase_char = rbase[0];                             char rbase_char2 = rbase2[0];
                char vbase_char = vbase[0];                             char vbase_char2 = vbase2[0];
		// het_string
                int rnum = vgraph[het_string].base_dict[rbase_char];    int rnum2 = vgraph2[het_string].base_dict[rbase_char2];
                int vnum = vgraph[het_string].base_dict[vbase_char];    int vnum2 = vgraph2[het_string].base_dict[vbase_char2];
                int total_normal = rnum + vnum;
		double rfrac = 0.0;
		double vfrac = 0.0;
		double min_fraction = 0.0;
		bool isvariant_normal = vgraph[het_string].var;
		bool isvariant_tumor = vgraph2[het_string].var;
		//cout << het_string << " " << het_string2 << " " << loop_pos << " " << rbase_char << " " << vbase_char << " " << rbase_char2 << " " << vbase_char2 << " " << rnum << " " << vnum << " " << rnum2 << " " << vnum2 << endl;
		if ( total_normal > 0 ) {
			rfrac = (double) rnum / (double) total_normal;
                	vfrac = (double) vnum / (double) total_normal;
			min_fraction = std::min(rfrac,vfrac);
		}
		if ( min_fraction > filter_fraction ) {
			if ( total_normal > min_total_cov && total_normal < max_total_cov ) {
				pdict.within_filter[j] = true;     pdict2.within_filter[j] = true;
				vgraph[het_string_variant].filter = true;  vgraph2[het_string_variant].filter = true;
				number_filter += 1;
			}
		}
                cout << het_string << " " << het_string_variant << " " << isvariant_normal << " " << isvariant_tumor << " " << loop_pos << " " << rbase << " " << vbase << " " << rnum << " " << vnum << "\t" <<rbase2 << " " << vbase2 << " " << rnum2 << " " << vnum2 << "  " << min_fraction << endl;
		//} //}
        }
	cout << "number filtered " << number_filter << endl;
	cout << "total number " << pdict.num_paired << endl;
}

////////////////////////////////////////////////////
void run_variant_filtering( int argc, char** argv ) {
        parse_filtering_options(argc, argv);
        std::string coverageFile = "./output/filtered_coverage_" + opt::id_string + ".dat";
	/////////// load normal coverage
        variant_graph vgraph;
        coord_dictionary pdict;
        read_het_coverage( opt::normal_cov_file, vgraph );
	initialize_pdict( vgraph, pdict, false );
	/////////// load tumor coverage
	variant_graph vgraph2;
	coord_dictionary pdict2;
	if(opt::filter_tumor) { read_het_coverage( opt::tumor_cov_file, vgraph2 ); }
	else { read_het_coverage( opt::normal_cov_file, vgraph2 ); }
	initialize_pdict( vgraph2, pdict2, false );
	///////////
	filter_het_sites( pdict, pdict2, vgraph, vgraph2 );
        //##############################################
	//cout << coverageFile << endl;
	write_het_coverage_filter( vgraph2, coverageFile );
        return;
};


//double vfrac = (double) vnum / (double) total_normal;
