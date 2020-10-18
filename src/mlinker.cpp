//####################################
#include "mlinker.h"

static const char *LINKER_USAGE_MESSAGE =
"############### running linker ############### \n"
"Program: mLinker \n"
"Version: 0.42 \n"
"Usage: mlinker <command> [options]\n\n"
"Commands:\n"
"\textract\t\textract links \n"
"\tmatrix\t\tbin genome and create linked read heatmap \n"
"\tsolve\t\tsolve germline phase from extracted links \n"
"\tscaffold\tcombine blocks and hic to phase whole chromosomes \n"
"\trecover\t\trecover heterozygous sites lost in scaffold \n"
"\tcoverage\textract coverage of each het site from bam file \n"
"\tcn_phase\tfurther phase genome due to LOH or aneuplody in genome \n"
"\tsv_phase\tphase sv call on long reads \n"
"\tfilter\t\tfilter normal variants by copy fraction \n"
"\n";
//"	phase		(extract + solve) phase germline haplotype from linked read sample \n"
//"	hic_phase	phase split reads in hic \n"


int main(int argc, char** argv) {
  if (argc <= 1) {
    std::cerr << LINKER_USAGE_MESSAGE;
    return 0;
  } else {
    std::string command(argv[1]);
    if (command == "help" || command == "--help") {
      std::cerr << LINKER_USAGE_MESSAGE;
      return 0;
    //} else if (command == "phase") { 
	//run_link_phaser(argc-1, argv+1);
    } else if (command == "extract") { 
	run_extract_hash(argc-1, argv+1);
    } else if (command == "solve") {
        run_solve_hap(argc-1, argv+1);
    } else if (command == "cn_phase") { 
	run_cn_phaser(argc-1, argv+1);
    } else if (command == "coverage") { 
	run_het_coverage(argc-1, argv+1);
    } else if (command == "matrix") {  	
	run_bin_matrix(argc-1, argv+1);
    } else if (command == "sv_phase") { 
	run_sv_phaser(argc-1, argv+1); 	
    } else if (command == "filter") {
        run_variant_filtering(argc-1, argv+1);
    //} else if (command == "hic_phase") {
    //    run_hic_phaser(argc-1, argv+1);
    } else if (command == "scaffold") {
        run_scaffold(argc-1, argv+1);
    } else if (command == "recover") {
        run_recover(argc-1, argv+1);
    }
    else {
      std::cerr << LINKER_USAGE_MESSAGE;
      return 0;
    }
  }

  std::cerr << "done with linker" << std::endl;
  return 0;
}

//"	bx_bin          get unique bx coverage from 10X \n"
    //} else if (command == "bx_bin") {
    //    run_bin_bxtag(argc-1, argv+1);

