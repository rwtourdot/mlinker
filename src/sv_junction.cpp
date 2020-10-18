#include "sv_junction.h"

void sv_entry::set_locations( std::string chr1, int pos1, std::string chr2, int pos2, int str1, int str2, int tcov ) {
	chr_one = chr1; chr_two = chr2;
	pos_one = pos1; pos_two = pos2;
	strand.first = str1;
	strand.second = str2;
	tot_coverage = tcov;
}

void sv_entry::add_tag_read1( std::string tag, std::string bx ) { 
	readid_tags1.push_back(tag); 
	bxmap1[tag] = bx;
}

void sv_entry::add_tag_read2( std::string tag, std::string bx ) { 
	readid_tags2.push_back(tag);
	bxmap2[tag] = bx; 
}

void sv_entry::read_id_intersection() {
        std::set<std::string> set_read1(readid_tags1.begin(),readid_tags1.end());
        std::set<std::string> set_read2(readid_tags2.begin(),readid_tags2.end());
        std::set_intersection(set_read1.begin(),set_read1.end(),set_read2.begin(),set_read2.end(),std::back_inserter(read_intersection));
        std::set_union(set_read1.begin(),set_read1.end(),set_read2.begin(),set_read2.end(),std::back_inserter(read_union));
        intersection_length = read_intersection.size();
        union_length = read_union.size();
	for (int j = 0; j < intersection_length; j++) { bx_intersection.push_back(bxmap1[read_intersection[j]]); }
}

//void sv_entry::sv_bx_intersection() {
//        std::set<std::string> set_sv1(sv1_bxtags.begin(),sv1_bxtags.end());
//        std::set<std::string> set_sv2(sv2_bxtags.begin(),sv2_bxtags.end());
//        std::set_intersection(set_sv1.begin(),set_sv1.end(),set_sv2.begin(),set_sv2.end(),std::back_inserter(tag_intersection));
//        std::set_union(set_sv1.begin(),set_sv1.end(),set_sv2.begin(),set_sv2.end(),std::back_inserter(tag_union));
//        intersection_length = tag_intersection.size();
//        union_length = tag_union.size();
//}

void sv_entry::het_map( int pos, int haplotype, bool side ) {
	if (side) { 
		hap_end1[pos] = haplotype;
		if ( num_end1.find(pos) != num_end1.end() ) { num_end1[pos] += 1; }
		else { num_end1[pos] = 1; }
	}
	else { 
		hap_end2[pos] = haplotype; 
		if ( num_end2.find(pos) != num_end2.end() ) { num_end2[pos] += 1; }
		else { num_end2[pos] = 1; }
	}
}



