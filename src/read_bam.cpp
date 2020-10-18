#include "read_bam.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////
bool alignment_boolean_check( BamTools::BamAlignment al ) {
        bool return_bool = false;
        if (al.IsPrimaryAlignment()) {
                if (al.AlignmentFlag < samsum_lt) {
                        if (!al.IsDuplicate()) {
                                if (al.MapQuality > minimum_mapq) {
                                        return_bool = true;
                                }
                        }
                }
        }
        return return_bool;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
bool alignment_boolean_check_hic( BamTools::BamAlignment al ) {
        bool return_bool = false;
        if (!al.IsDuplicate()) { 
        	if (al.MapQuality > minimum_mapq_hic) {
			return_bool = true;
		} 
	}
        return return_bool;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
bool alignment_boolean_check_lenient( BamTools::BamAlignment al ) {
        bool return_bool = false;
        if (!al.IsDuplicate()) { return_bool = true; }
        return return_bool;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
std::string return_read_id( std::string alignment_name ) {
        std::string delimiter = "/";
        std::string func_read_string;
        unsigned last = alignment_name.find_last_of(delimiter);
        func_read_string = alignment_name.substr(last+1,alignment_name.length()-last);
        return func_read_string;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
std::pair<char,int> get_base_readpos( std::vector<BamTools::CigarOp> &cigar_str, std::string aligned_bases, int variant_pos, int start_pos, std::string qual ) {
        std::pair<char,int> base_qual;
        std::vector<int> seq_map(aligned_bases.length());
        std::vector<int> qmap_backward(aligned_bases.length()); //std::vector<int> qmap(qual.size()); //for (int k=0; k<qmap.size(); k++) { qmap[k] = -1; }
        std::vector<int> qvector; for (char& ch : qual) { int xord = ch; int qscore = (xord-33); qvector.push_back(qscore); }
        char base_to_return;
        int quality_to_return;  int loop_until = 0;  int s_index = 0;  int num_bases_on_read=0;
        for (int k=0; k < (int)cigar_str.size(); k++) { //cout << cigar_str[k].Length << cigar_str[k].Type << ' ';
                if ((cigar_str[k].Type != 'S') && (cigar_str[k].Type != 'H')) {
                        int qindex = num_bases_on_read;
                        for (int l=0; l < (int)cigar_str[k].Length; l++) { //cout << aligned_bases.at(s_index); //cout << k << "  " << "  " << s_index << " " << endl;
                                //qmap[qindex] = s_index;
                                qmap_backward[s_index] = qindex;
                                seq_map[s_index] = start_pos;
                                if (start_pos == variant_pos) { loop_until = s_index; }
                                if ((cigar_str[k].Type == 'M') || (cigar_str[k].Type == 'D')) { start_pos++; }
                                qindex++;
                                s_index++;
                        }
                }
                if (cigar_str[k].Type != 'D') { num_bases_on_read += cigar_str[k].Length; }
        }  //cout << endl;
        for (int l=0; l < (loop_until+1); l++) {
                if (seq_map[l] == variant_pos) {    //cout << l << "   " << aligned_bases[l] << endl;
                        base_to_return = aligned_bases[l];
                        if (base_to_return == '-') { base_to_return = 'D'; quality_to_return = minimum_baseq; }   //check if site is a deletion
                        else { quality_to_return = qvector[qmap_backward[l]]; }
                        if (std::count(seq_map.begin(),seq_map.end(),variant_pos) > 1) { base_to_return = 'I';  }  //check if site is a indel
                }
        }    //cout << base_to_return << "  " << quality_to_return << endl;
        base_qual.first = base_to_return;
        base_qual.second = quality_to_return;   //cin.get();
        return base_qual;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
std::string return_index_string( std::string filename ) {
        std::string index_filename;
        index_filename = filename.substr(0,filename.size()-1);
        index_filename += "i";
        return index_filename;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
std::string technology_hash( std::string tech, BamTools::BamAlignment al ) {
        std::string readhash_long;
        if ( tech == "pacbio" || tech == "nanopore" ) { //|| tech == "illumina" || tech == "hic"
                std::string readname,start_pos,end_pos;
                std::vector<std::string> splitted_string;   splitted_string.reserve(10);
                splitted_string = split_string(al.Name);
                readname = splitted_string.back();          //readname should
                start_pos = std::to_string(al.Position);
                end_pos = std::to_string(al.GetEndPosition());
                readhash_long = readname + "_" + start_pos + "_" + end_pos; // readname should be bxn
        }
	else if ( tech == "illumina" || tech == "hic") {
		std::string readname;
		std::vector<std::string> splitted_string;   splitted_string.reserve(10);
                splitted_string = split_string(al.Name);
                readname = splitted_string.back();
		readhash_long = readname;
	}
        else if ( tech == "tenx" ) {
                std::string bxtag = "BX";
                std::string return_bx_string;
                al.GetTag(bxtag,return_bx_string);
                readhash_long = return_bx_string;
        }
        return readhash_long;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
bool get_tag_value( const std::string &tag, variant_graph &variant_cnx ) {
        auto check = variant_cnx.find(tag);
        if (check == variant_cnx.end()) return false;
        return true;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
static void open_bam_file( BamTools::BamReader &reader, std::string inputFilename ) {
        std::string indexFilename;
        indexFilename = return_index_string( inputFilename );
        if (!reader.Open(inputFilename)) { cerr << "Could not open input BAM file." << endl; return; }
        if (!reader.OpenIndex(indexFilename)) {
                indexFilename = inputFilename + ".bai";
                cout << "looking for index file " << indexFilename << endl;
                if (!reader.OpenIndex(indexFilename)) { cerr << "Could not find BAM index file." << endl; return; }
        }
        reader.Open(inputFilename);
        reader.OpenIndex(indexFilename);        //BamTools::BamRegion BamRegion;
	return;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
void connect_up_variants_bam_pileup( vcf_vector vvec, std::string inputFilename, int chromosome, variant_graph &variant_cnx, read_graph &rgraph, std::string technology ) {   ////// bamtools prep ///////
        BamTools::BamReader reader;
	open_bam_file(reader,inputFilename);
        BamTools::BamAlignment al;
        ///////////////////////////////////
        std::string variant_node_name,readhash_long;
        char base_at_pos;
        int total_variants = 0;
	int used_variants = 0;
        clock_t t;
        t = clock();
        for (int i=0; i < vvec.size(); i++) {
                if ((vvec[i].ref_base.length() == 1 ) && (vvec[i].var_base.length() == 1 )) {
		if (vvec[i].bounded) {
                //cout << vcf_vector[i].pos << "    " << vcf_vector[i].ref_base << "   " << vcf_vector[i].var_base << endl;
		int loop_pos = vvec[i].pos;
                std::string rbase = vvec[i].ref_base;
                std::string vbase = vvec[i].var_base;
                std::string variant_id = std::to_string(loop_pos) + "_" + rbase + "_" + vbase;
                std::string reference_id = std::to_string(loop_pos) + "_" + rbase + "_" + rbase;
                bool present1 = get_tag_value(variant_id,variant_cnx);
                bool present2 = get_tag_value(reference_id,variant_cnx);
                if(!present1) {
                        variant_node v_node;  
			variant_cnx[variant_id] = v_node;
                        variant_cnx[variant_id].set_values(loop_pos,true,vbase,rbase);
                }
                if(!present2) {
                        variant_node r_node;  
			variant_cnx[reference_id] = r_node;
                        variant_cnx[reference_id].set_values(loop_pos,false,vbase,rbase);
                }
                BamTools::BamRegion region(chromosome,loop_pos-2,chromosome,loop_pos+2);
                reader.SetRegion(region);
                while(reader.GetNextAlignmentCore(al)) {
                        if ( (al.Position < loop_pos) && (al.GetEndPosition() > loop_pos) ) {
			bool al_check = false;
			if ( technology != "hic" ) { al_check = alignment_boolean_check(al); }
			else { al_check = alignment_boolean_check_hic(al); }	
                        if ( al_check == true ) {  // added in 0.5
                                al.BuildCharData();
                                std::string qual = al.Qualities;
				std::string readname = al.Name;
                                std::pair<char,int> base_qual = get_base_readpos(al.CigarData,al.AlignedBases,loop_pos,al.Position,qual);
                                base_at_pos = base_qual.first;
                                int bqual = base_qual.second;
				//if (bqual >= minimum_baseq_hic) {
				bool bq_check = false;
				if ( technology == "hic" ) { if ( bqual >= minimum_baseq_hic ) { bq_check = true; } }
				else if ( technology == "tenx" ) { if ( bqual >= minimum_baseq_tenx ) { bq_check = true; } }
				else { if ( bqual >= minimum_baseq ) { bq_check = true; } }
                                if ( bq_check == true ) {
                                std::string readhash_long = technology_hash(technology,al);
                                std::stringstream ss;   std::string base_to_compare;
                                ss << base_at_pos;      ss >> base_to_compare;  // convert character to string
                                //cout << al.RefID << "  " << al.MapQuality << "  " << start_pos << "  " << end_pos << "   " << base_at_pos << endl;
				//cout << technology << " " << readhash_long << endl;
				if (readhash_long.size() > 1 ) {
                                        variant_cnx[variant_id].add_base_to_dict(base_at_pos,readhash_long);
                                        variant_cnx[reference_id].add_base_to_dict(base_at_pos,readhash_long);
				} else {
                                        variant_cnx[variant_id].add_base_to_dict(base_at_pos,"nohash");
                                        variant_cnx[reference_id].add_base_to_dict(base_at_pos,"nohash");
				}
                                if (readhash_long.size() > 1 ) {           // tenx if loop - sometimes reads dont have bx tags  readname long is empty
                                if ((base_to_compare == vbase) || (base_to_compare == rbase)) {
                                        std::string het_hash = std::to_string(loop_pos) + "_" + rbase + "_" + base_to_compare;
                                        if (rgraph.find(readhash_long) != rgraph.end()) { rgraph[readhash_long].add_connection(het_hash,readname); }
                                        else { rgraph[readhash_long].set_name(readhash_long); rgraph[readhash_long].add_connection(het_hash,readname); }
                                        if (base_to_compare == vbase) {      variant_cnx[variant_id].add_connected_read(readhash_long,readname); }
                                        else if (base_to_compare == rbase) { variant_cnx[reference_id].add_connected_read(readhash_long,readname); }
                                } } }
                        } }
                }
		used_variants += 1;
                } }
                total_variants += 1;
                //cout << iterate << "    " << vcf_vector.size() << endl;
        }
        t = clock() - t;
        cout << "loaded bamreads ---------------- time: " << t << endl;
        cout << "total variants (indels and unbounded): " << total_variants << endl;
        cout << "useful phasing variants -------------: " << used_variants << endl;
        reader.Close();
        return;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
void connect_up_variants_hic_bam_pileup( coord_dictionary& pdict, std::string inputFilename, int chromosome, variant_graph &variant_cnx, read_graph &rgraph, std::string technology ) {   ////// bamtools prep ///////
        BamTools::BamReader reader;
        open_bam_file(reader,inputFilename);
        BamTools::BamAlignment al;
        ///////////////////////////////////
        std::string readhash_long;
        char base_at_pos;  //int total_variants = 0;
        int used_variants = 0;
        clock_t t;
        t = clock();
	for (int i = 0; i < pdict.num_paired; i++) {
                int loop_pos = pdict.sorted_paired_positions[i];
                std::string het_string = pdict.paired_dict[loop_pos][0];
                std::string het_string2 = pdict.paired_dict[loop_pos][1];
		std::string rbase = variant_cnx[het_string].ref_base;
		std::string vbase = variant_cnx[het_string].var_base;
		std::string variant_id,reference_id;
                if ( variant_cnx[het_string].var ) { variant_id = het_string;  reference_id = het_string2; } 
		else {              		     variant_id = het_string2; reference_id = het_string;  }
		//cout << loop_pos << " " << reference_id << " " << variant_id << endl;
                BamTools::BamRegion region(chromosome,loop_pos-2,chromosome,loop_pos+2);
                reader.SetRegion(region);
                while(reader.GetNextAlignmentCore(al)) {
                	if ((al.Position < loop_pos) && (al.GetEndPosition() > loop_pos)) {
                        if (alignment_boolean_check_hic(al)) {  // added in 0.5   // added a hic specific check 
                                al.BuildCharData();
                                std::string qual = al.Qualities;
				std::string readname = al.Name;
                                std::pair<char,int> base_qual = get_base_readpos(al.CigarData,al.AlignedBases,loop_pos,al.Position,qual);
                                base_at_pos = base_qual.first;
                                int bqual = base_qual.second;
                                if (bqual >= minimum_baseq_hic) {
                                std::string readhash_long = technology_hash(technology,al);
                                std::stringstream ss;   std::string base_to_compare;
                                ss << base_at_pos;      ss >> base_to_compare;  // convert character to string
                                //cout << al.RefID << "  " << al.MapQuality << "  " << start_pos << "  " << end_pos << "   " << base_at_pos << endl;
                                //cout << technology << " " << readhash_long << endl;
                                //if (readhash_long == "H7CMTADXX:1:1108:9261:79691") { cout << readhash_long << endl; }
                                if (readhash_long.size() > 1 ) {
                                        variant_cnx[variant_id].add_base_to_dict(base_at_pos,readhash_long);
                                        variant_cnx[reference_id].add_base_to_dict(base_at_pos,readhash_long);
                                } else {
                                        variant_cnx[variant_id].add_base_to_dict(base_at_pos,"nohash");
                                        variant_cnx[reference_id].add_base_to_dict(base_at_pos,"nohash");
                                }
                                if (readhash_long.size() > 1 ) {           // tenx if loop - sometimes reads dont have bx tags  readname long is empty
                                if ((base_to_compare == vbase) || (base_to_compare == rbase)) {
                                        std::string het_hash = std::to_string(loop_pos) + "_" + rbase + "_" + base_to_compare;
                                        if (rgraph.find(readhash_long) != rgraph.end()) { rgraph[readhash_long].add_connection(het_hash,readname); }
                                        else { rgraph[readhash_long].set_name(readhash_long); rgraph[readhash_long].add_connection(het_hash,readname); }
                                        if (base_to_compare == vbase) {      variant_cnx[variant_id].add_connected_read(readhash_long,readname); }
                                        else if (base_to_compare == rbase) { variant_cnx[reference_id].add_connected_read(readhash_long,readname); }
                                } } }
                        } }
                }
                used_variants += 1;
        }
        t = clock() - t;
        cout << "loaded bamreads ---------------- time: " << t << endl;
        cout << "useful phasing variants -------------: " << used_variants << endl;
        reader.Close();
        return;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
void connect_assembly_bam_pileup( std::string inputFilename, bx_map &gen_map, read_graph &bx_graph, int binsize, full_map& chr_map, std::unordered_map<std::string,int>& contig_dict ) {   /////// bamtools prep /////////////
        BamTools::BamReader reader;
	open_bam_file(reader,inputFilename);
        BamTools::BamAlignment al;
        std::string start_pos,end_pos;
        int iterate = 0;
        int contig_intref,num_ref_count;
        clock_t t;   t = clock();
	for (auto& it : contig_dict) {
		std::string contig = it.first;
		int contig_length = it.second;
        	contig_intref = reader.GetReferenceID(contig);
        	num_ref_count = reader.GetReferenceCount();     //cout << contig_intref << endl;
        	BamTools::BamRegion region(contig_intref);
        	reader.SetRegion(region.LeftRefID,0,region.LeftRefID,contig_length);
        	while(reader.GetNextAlignment(al)) {   //al.BuildCharData();
                	if (al.RefID == contig_intref) {
                	if (alignment_boolean_check(al)) {
                        	std::string bxtag = "BX";
                        	std::string return_bx_string;
                        	al.GetTag(bxtag,return_bx_string);
                        	start_pos = std::to_string(al.Position);
                        	int index = std::floor(al.Position/binsize);  // this int needs to be converted to a string
                        	if (return_bx_string.size() > 1) {
                                	gen_map[index].add_connected_read(return_bx_string);
                                	chr_map[contig][index].add_connected_read(return_bx_string);
                                	//cout << al.RefID << "   " << index << "  " << start_pos << "  " << return_bx_string << endl;
                                	if (bx_graph.find(return_bx_string) != bx_graph.end()) { bx_graph[return_bx_string].add_connection_int(index,contig); }
                                	else { bx_graph[return_bx_string].set_name(return_bx_string); bx_graph[return_bx_string].add_connection_int(index,contig); }
                        	}
                        	iterate += 1;
                	} }
        	}
	}
        t = clock() - t;
        cout << "loaded bx tags - time: " << t << endl;
        cout << "number of reads: " << iterate << endl;
        reader.Close();
        return;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
void connect_sv_read_pileup( std::string inputFilename, read_graph &bx_graph, std::vector<sv_entry>& sv_list, std::map<std::string,int> chr_str_map ) {   /////// bamtools prep /////////////
        BamTools::BamReader reader;
        open_bam_file(reader,inputFilename);
        BamTools::BamAlignment al;
        std::string start_pos,end_pos;
        clock_t t;   t = clock();
	cout << "reading bamfile - sv parsing" << endl;
        for (int i = 0; i < sv_list.size(); i++) {  // 1 instead of 0
		//cout << sv_list[i].chr_one << "  " << sv_list[i].pos_one << "  " << sv_list[i].strand.first << "  " << sv_list[i].chr_two << "  " << sv_list[i].pos_two << "  " << sv_list[i].strand.second << endl;
                int chromosome1 = chr_str_map[sv_list[i].chr_one];
                int chromosome2 = chr_str_map[sv_list[i].chr_two];
		int start1,start2,end1,end2; 
		if ( sv_list[i].strand.first == 1 )        { start1 = sv_list[i].pos_one - sv_buffer; end1 = sv_list[i].pos_one + 2; }
	 	else if ( sv_list[i].strand.first == -1 )  { start1 = sv_list[i].pos_one - 2;         end1 = sv_list[i].pos_one + sv_buffer; }  
		if ( sv_list[i].strand.second == 1 ) 	   { start2 = sv_list[i].pos_two - sv_buffer; end2 = sv_list[i].pos_two + 2; }
		else if ( sv_list[i].strand.second == -1 ) { start2 = sv_list[i].pos_two - 2;         end2 = sv_list[i].pos_two + sv_buffer; }  
                BamTools::BamRegion region1( chromosome1, start1, chromosome1, end1 );
                BamTools::BamRegion region2( chromosome2, start2, chromosome2, end2 );
                reader.SetRegion(region1);
                while(reader.GetNextAlignment(al)) { 
                        if (alignment_boolean_check_lenient(al)) {
                                std::string return_bx_string;
                                al.GetTag("BX",return_bx_string);
                                std::string readname = al.Name;    
                                if  (return_bx_string.size() > 1) { sv_list[i].add_tag_read1(readname,return_bx_string); } 
                        }
                }
                reader.SetRegion(region2);
                while(reader.GetNextAlignment(al)) { 
                        if (alignment_boolean_check_lenient(al)) {
                                std::string return_bx_string;
                                al.GetTag("BX",return_bx_string);
				std::string readname = al.Name;
                                if  (return_bx_string.size() > 1) { sv_list[i].add_tag_read2(readname,return_bx_string); }
                        }
                }
        }
        reader.Close();
        for (int i = 0; i < sv_list.size(); i++) { sv_list[i].read_id_intersection(); }
        return;
                //int contig_intref = reader.GetReferenceID(sv_list[i].chr_one);
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
void sv_phase_het( BamTools::BamReader& reader, BamTools::BamAlignment& al, sv_entry& sv_object, int start, int end, int chromosome, vcf_entry vnode, int& possible_hets, bool side ) {
	int vpos = vnode.pos;
	int vchr = vnode.chromosome_id;
	if ( vchr == chromosome && vpos > start && vpos < end ) {
		possible_hets += 1;
		BamTools::BamRegion het_region( chromosome, vpos-2, chromosome, vpos+2 );
		reader.SetRegion(het_region);
		while(reader.GetNextAlignment(al)) {
			std::string return_bx_string;
			al.GetTag("BX",return_bx_string);
			std::vector<std::string> bxvector = sv_object.bx_intersection;
			if (std::find(bxvector.begin(), bxvector.end(), return_bx_string) != bxvector.end()) {
				if ((al.Position < vpos) && (al.GetEndPosition() > vpos)) {
					std::string qual = al.Qualities;
					std::pair<char,int> base_qual = get_base_readpos(al.CigarData,al.AlignedBases,vpos,al.Position,qual);
					char base_at_pos = base_qual.first;
					std::stringstream ss;
					std::string base_at_pos_str;
					ss << base_at_pos;
					ss >> base_at_pos_str;
					if ( base_at_pos_str == vnode.ref_base ) { sv_object.het_map(vpos,1,side); }
					else if ( base_at_pos_str == vnode.var_base ) { sv_object.het_map(vpos,-1,side); }
					cout << " matched bx read " << vchr << "  " << vpos << "  " << vnode.ref_base << "_" << vnode.var_base << "  " << base_at_pos << " - " << side << endl;
				}
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
void connect_sv_hets( vcf_vector vvec, std::string inputFilename, read_graph &bx_graph, std::vector<sv_entry>& sv_list, std::map<std::string,int> chr_str_map ) {   /////// bamtools prep ///////////
        BamTools::BamReader reader;
        open_bam_file(reader,inputFilename);
        BamTools::BamAlignment al;
        for (int i = 0; i < sv_list.size(); i++) {
		int chromosome1 = chr_str_map[sv_list[i].chr_one];
                int chromosome2 = chr_str_map[sv_list[i].chr_two];
                int start1,start2,end1,end2;
                if ( sv_list[i].strand.first == 1 )        { start1 = sv_list[i].pos_one - sv_bx_buffer; end1 = sv_list[i].pos_one + 2; }
                else if ( sv_list[i].strand.first == -1 )  { start1 = sv_list[i].pos_one - 2;            end1 = sv_list[i].pos_one + sv_bx_buffer; }
                if ( sv_list[i].strand.second == 1 )       { start2 = sv_list[i].pos_two - sv_bx_buffer; end2 = sv_list[i].pos_two + 2; }
                else if ( sv_list[i].strand.second == -1 ) { start2 = sv_list[i].pos_two - 2;            end2 = sv_list[i].pos_two + sv_bx_buffer; }
		int possible_hets1 = 0;
		int possible_hets2 = 0;
		for (int j=0; j < vvec.size(); j++) {
			vcf_entry vnode = vvec[j];
			if ((vnode.ref_base.length() == 1 ) && (vnode.var_base.length() == 1 )) {
				sv_phase_het(reader,al,sv_list[i],start1,end1,chromosome1,vnode,possible_hets1,true);
				sv_phase_het(reader,al,sv_list[i],start2,end2,chromosome2,vnode,possible_hets2,false);
			} 
		}	
		cout << sv_list[i].chr_one << "  " << sv_list[i].pos_one << "  " << sv_list[i].strand.first << "  " << sv_list[i].chr_two << "  " << sv_list[i].pos_two << "  " << sv_list[i].strand.second << "  " << sv_list[i].bx_intersection.size() << " phets " << possible_hets1 << endl;
	}
	reader.Close();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
void contig_name_map( std::string inputFilename, std::map<std::string,int> &chr_str_map )  {     /////// bamtools prep /////////////
        BamTools::BamReader reader;
	open_bam_file(reader,inputFilename);
        BamTools::RefVector refnames = reader.GetReferenceData();
        int bi = 0;  for (auto& it : refnames) { chr_str_map[it.RefName] = bi; bi++; }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
void contig_name_length( std::string inputFilename, std::unordered_map<std::string,int> &contig_dict )  {    /////// bamtools prep /////////////
        BamTools::BamReader reader;
	open_bam_file(reader,inputFilename);
        BamTools::RefVector refnames = reader.GetReferenceData();
        for (auto& it : refnames) { contig_dict[it.RefName] = it.RefLength; }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<std::string> split_string( std::string teststring ) {
        std::stringstream test(teststring);
        std::string segment;
        std::vector<std::string> seglist;
        while(std::getline(test,segment,'/')) { seglist.push_back(segment); }
        return seglist;
};













                //if ((vvec[i].ref_base.length() == 1 ) && (vvec[i].var_base.length() == 1 )) {
                //if (vvec[i].bounded) {
                        //cout << vcf_vector[i].pos << "    " << vcf_vector[i].ref_base << "   " << vcf_vector[i].var_base << endl;
                //std::string variant_id = std::to_string(vvec[i].pos) + "_" + vvec[i].ref_base + "_" + vvec[i].var_base;
                //std::string reference_id = std::to_string(vvec[i].pos) + "_" + vvec[i].ref_base + "_" + vvec[i].ref_base;
                //bool present1 = get_tag_value(variant_id,variant_cnx);
                //bool present2 = get_tag_value(reference_id,variant_cnx);
                //if(!present1) {
                //      variant_node v_node;
                //      variant_cnx[variant_id] = v_node;
                //      variant_cnx[variant_id].set_values(vvec[i].pos,true,vvec[i].var_base,vvec[i].ref_base);
                //}
                //if(!present2) {
                //      variant_node r_node;
                //      variant_cnx[reference_id] = r_node;
                //      variant_cnx[reference_id].set_values(vvec[i].pos,false,vvec[i].var_base,vvec[i].ref_base);
                //}









                //BamTools::BamRegion region1( chromosome1, start1, chromosome1, end1 );
                //BamTools::BamRegion region2( chromosome2, start2, chromosome2, end2 );


//for (int k = 0; k < sv_list[i].bx_intersection.size(); k++) {
//      cout << sv_list[i].bx_intersection[k] << endl;
//}


////////////////////////////////////////////////////////////////////////////////////////////////////////
//void connect_sv_bam_pileup( std::string inputFilename, read_graph &bx_graph, std::vector<sv_entry>& sv_list, std::map<std::string,int> chr_str_map ) {   /////// bamtools prep /////////////
//        BamTools::BamReader reader;
//        open_bam_file(reader,inputFilename);
//        BamTools::BamAlignment al;
//        std::string start_pos,end_pos;
//        clock_t t;   t = clock();
//      for (int i = 0; i < sv_list.size(); i++) {
//              cout << sv_list[i].chr_one << "  " << sv_list[i].pos_one << "  " << sv_list[i].chr_two << "  " << sv_list[i].pos_two << endl;
//              int chromosome1 = chr_str_map[sv_list[i].chr_one];
//              int chromosome2 = chr_str_map[sv_list[i].chr_two];
//              int start1 = sv_list[i].pos_one - 2; int start2 = sv_list[i].pos_two - 2;
//              int end1 = sv_list[i].pos_one + 2;   int end2 = sv_list[i].pos_two + 2;
//                BamTools::BamRegion region1( chromosome1, start1, chromosome1, end1 );
//                BamTools::BamRegion region2( chromosome2, start2, chromosome2, end2 );
//                reader.SetRegion(region1);
//              while(reader.GetNextAlignment(al)) {   //if (al.RefID == contig_intref) {
//                      if (alignment_boolean_check(al)) {
//                              std::string bxtag = "BX";
//                              std::string return_bx_string;
//                              al.GetTag(bxtag,return_bx_string);
//                              if (return_bx_string.size() > 1) { sv_list[i].add_tag_sv1(return_bx_string); }
//                      }
//              }
//                reader.SetRegion(region2);
//                while(reader.GetNextAlignment(al)) {   //if (al.RefID == contig_intref) {
//                        if (alignment_boolean_check(al)) {
//                                std::string bxtag = "BX";
//                                std::string return_bx_string;
//                                al.GetTag(bxtag,return_bx_string);
//                                if (return_bx_string.size() > 1) { sv_list[i].add_tag_sv2(return_bx_string); }
//                        }
//                }
//      }
//      reader.Close();
//        return;
                //int contig_intref = reader.GetReferenceID(sv_list[i].chr_one);
//};




