#include "hic_links.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////
std::string split_string_first_hic_link( std::string s, std::string delimiter, int choice ) {
        std::string return_string;
        if ( choice == 0 ) { return_string = s.substr(0, s.find(delimiter)); }
        if ( choice == 1 ) { return_string = s.substr(s.find(delimiter)+1, s.length()); }
        return return_string;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////
void hic_link::add_link( int pos1, int pos2, std::string var1, std::string var2, int depth ) {
        anchor_pos1.push_back(pos1);  //hlink
        anchor_pos2.push_back(pos2);
	anchor_var1.push_back(var1);
	anchor_var2.push_back(var2);
	link_depth.push_back(depth);
	std::string var1_temp_string = split_string_first_hic_link(var1,"_",1);
	std::string var2_temp_string = split_string_first_hic_link(var2,"_",1);
	std::string var1_ref_base = split_string_first_hic_link(var1_temp_string,"_",0);
	std::string var1_read_base = split_string_first_hic_link(var1_temp_string,"_",1);
	std::string var2_ref_base = split_string_first_hic_link(var2_temp_string,"_",0);
	std::string var2_read_base = split_string_first_hic_link(var2_temp_string,"_",1);
	anchor_ref_base1.push_back(var1_ref_base);
	anchor_read_base1.push_back(var1_read_base);
	anchor_ref_base2.push_back(var2_ref_base);
	anchor_read_base2.push_back(var2_read_base);
	hic_pos_map[pos1].push_back(pos2);
	hic_link_map[var1].push_back(var2);
	hic_link_depth_map[var1].push_back(depth);
	if ( var1_ref_base != var1_read_base ) { anchor_var_bool1.push_back(true); }
	else { anchor_var_bool1.push_back(false); }
	if ( var2_ref_base != var2_read_base ) { anchor_var_bool2.push_back(true); }
	else { anchor_var_bool2.push_back(false); }
};

