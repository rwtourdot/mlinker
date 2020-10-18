#include "run_solve_hap.h"

/////////////// namespaces ///////////////////
namespace opt {
  static std::string chr_choice = "chr20";
	static std::string input_graph_file = "./output/graph_variant_jun10_BL1954_tenx_chr20.dat";
  static std::string id_string = "default";
	static int start_bound = 0;
	static int end_bound = 300000000;
  static int block_cutoff = -700;
  static bool loh_mode = false;
  static bool output_net = false;
};

/////////////// structures ///////////////////
static const char* shortopts = "hvlo:i:c:n:b:";
static const struct option longopts[] = {
	{ "help",        no_argument, NULL, 'h' },
  { "graph-file",    no_argument, NULL, 'i' },
  { "chr-choice",  no_argument, NULL, 'c' },
  { "id_string",   no_argument, NULL, 'n' },
  { "block_cutoff",  no_argument, NULL, 'b' },
  { "loh_mode",   no_argument, NULL, 'l' }, // these were 1's
  { "output_net",   no_argument, NULL, 'v' }
};

///////////////////////////////////////////////////////////////////////////////////////////////////////
static const char *PHASER_USAGE_MESSAGE =
"Usage: mlinker solve [OPTION] -i /output/graph_variant_*.dat -c chr20 -n trial \n\n"
"\n"
"  Options\n"
"  -i,      input graph_variant file \n"
"  -c,      chromosome name ( chr4 ) or contig name depending on bam \n"
"  -n,      id string for output files \n"
"  -b,      blockE cutoff (enter positive, default -700) \n"
"  -l,      flag to run in LOH mode (variant sites need not be heterozygous - default no) \n"
"  -v,      flag to output variant network files (default no) \n"
"\n";

///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
static void parse_solve_hap_options( int argc, char** argv ) {
	bool die = false; //bool vcf_load = false; //bool cov_load = false; //if(argc <= 2) { die = true; }
  if (argc < 2) {
    std::cerr << "\n" << PHASER_USAGE_MESSAGE;
    exit(1);
  }
  if (string(argv[1]) == "help" || string(argv[1]) == "--help") {
    std::cerr << "\n" << PHASER_USAGE_MESSAGE;
    exit(1);
  }
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
      case 'h': die = true; break;
      case 'i': arg >> opt::input_graph_file; break;
      case 'c': arg >> opt::chr_choice; break;
      case 'n': arg >> opt::id_string; break;
      case 'b': arg >> opt::block_cutoff; break;
      case 'l': opt::loh_mode = true; break;
      case 'v': opt::output_net = true; break;
    }
  }
	if (die) {
	  std::cerr << "\n" << PHASER_USAGE_MESSAGE;
	  exit(1);
	}
  if (opt::block_cutoff > 0) { int temp_bc = opt::block_cutoff; opt::block_cutoff = -temp_bc; }
	//parse_region_string( opt::chr_choice );
  cout << endl;
  cout << "############### running link phaser ############### " << endl;
  cout << "== chromosome === " << opt::chr_choice << endl;
  cout << "== input graph  === " << opt::input_graph_file << endl;
  cout << "== id string  === " << opt::id_string << endl;
  cout << "== block_cutoff  === " << opt::block_cutoff << endl;
  cout << "== loh mode  === " << opt::loh_mode << endl;
  cout << "== output network  === " << opt::output_net << endl;
  cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx " << endl;
  cout << endl;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////
void run_solve_hap( int argc, char** argv ) {
  parse_solve_hap_options(argc, argv);
  std::string output_network_file = "./output/variant_network_" + opt::id_string + "_" + opt::chr_choice + ".dat";
  std::string hapsolutionFile = "./output/hap_solution_" + opt::id_string + "_" + opt::chr_choice + ".dat";
  if (opt::output_net) {  cout << "- output net file: " << output_network_file << endl; }
  cout << "- output hap file: " << hapsolutionFile << endl;
  cout << endl;
  coord_dictionary pdict;
  read_graph rgraph;
  variant_graph vgraph;
  float t = clock();
  //#################### start of code ##########################
	cout << "loading graph " << endl;
	t = clock();
	read_variant_graph_file(opt::input_graph_file,opt::chr_choice,vgraph,rgraph);
	t = float(clock() - t)/ CLOCKS_PER_SEC;
	cout << "read graph variant file == time: " << t << endl;
	cout << "=============== " << endl;
	cout << "linking graph " << endl;
	//#############################################################
       t = clock();
  link_hashes(vgraph,rgraph);
        t = float(clock() - t)/ CLOCKS_PER_SEC;
        cout << "link hashes == time: " << t << endl;
       t = clock();
	cout << "=============== " << endl;
	cout << "pruning graph " << endl;
  prune_graph(vgraph);
        t = float(clock() - t)/ CLOCKS_PER_SEC;
        cout << "prune graph == time: " << t << endl;
	cout << "=============== " << endl;
	cout << "initialize matrix coordinates " << endl;
       t = clock();
  initialize_pdict(vgraph,pdict,opt::loh_mode);
        t = float(clock() - t)/ CLOCKS_PER_SEC;
        cout << "init coord dict == time: " << t << endl;
	cout << "=============== " << endl;
	//#############################################################
  if (opt::output_net) { write_link_network(vgraph,output_network_file,opt::chr_choice); }
  static const std::size_t length = pdict.num_paired;
	///////////// create global datastructures for haplotype solver
  map_matrix<int> num_matrix_second(length);
	map_matrix<double> diff_matrix(length);
       t = clock();
  if (opt::loh_mode) { initialize_solver_loh(vgraph,pdict,diff_matrix,num_matrix_second); }
  else { initialize_solver(vgraph,pdict,diff_matrix,num_matrix_second); }
        t = float(clock() - t)/ CLOCKS_PER_SEC;
        cout << "init solver == time: " << t << endl;
	cout << "=============== " << endl;
	solver_recursive(vgraph,pdict,diff_matrix,num_matrix_second);
	//solver(vgraph,pdict,diff_matrix,num_matrix_second);
	///////////// write haplotype output
	cout << " solver finished " << endl;
	call_blocks(pdict,opt::block_cutoff);
	cout << " called haplotype blocks " << endl;
  if (opt::loh_mode) { write_hap_solution_loh(vgraph,hapsolutionFile,pdict,opt::chr_choice); }
  else { write_hap_solution(vgraph,hapsolutionFile,pdict,opt::chr_choice); }
  return;
};
