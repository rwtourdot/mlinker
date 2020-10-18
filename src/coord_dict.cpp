#include "coord_dict.h"

void coord_dictionary::initialize( std::unordered_map<std::string,variant_node> &var_dict , bool loh_mode ) {
    for (auto it : var_dict) {
        if (it.second.paired) {
            paired_dict[it.second.pos].push_back(it.first);
            unpaired_dict[it.second.pos] = it.first;
            double_positions.push_back(it.second.pos);  //char rbase = it.second.ref_base[0]; //char vbase = it.second.var_base[0];
	    paired[it.second.pos] = 1;
            if (it.second.var) { char vbase = it.second.var_base[0]; ref_index[it.second.pos] = 1; base_number[it.second.pos].push_back(it.second.base_dict[vbase]); }
            else {               char rbase = it.second.ref_base[0]; ref_index[it.second.pos] = 0; base_number[it.second.pos].push_back(it.second.base_dict[rbase]); }
        } else {  //cout << it.second.pos << " " << it.second.var << endl;
	    paired_dict[it.second.pos].push_back(it.first);
	    paired[it.second.pos] = 0;
            unpaired_dict[it.second.pos] = it.first;
            if (it.second.var) { char vbase = it.second.var_base[0]; ref_index[it.second.pos] = 1; base_number[it.second.pos].push_back(it.second.base_dict[vbase]); }
            else {               char rbase = it.second.ref_base[0]; ref_index[it.second.pos] = 0; base_number[it.second.pos].push_back(it.second.base_dict[rbase]); }
	}
        all_positions.push_back(it.second.pos);
    }
    cout << "double positions: " << double_positions.size() << endl;
    cout << "all positions: " << all_positions.size() << endl;
    std::sort( double_positions.begin(), double_positions.end() );
    /////////////////////////
    std::sort( all_positions.begin(), all_positions.end() );
    all_positions.erase( unique( all_positions.begin(), all_positions.end() ), all_positions.end() );
    //////////////////////////////////
    if (loh_mode) { for (std::vector<int>::iterator it= all_positions.begin(); it!= all_positions.end(); ++it) { sorted_paired_positions.push_back(*it); } }
    else { for (int i = 0; i < double_positions.size(); i+=2) { sorted_paired_positions.push_back(double_positions[i]); } }
    //////////////////////////////////
    for (int i = 0; i < all_positions.size(); i+=1) { sorted_all_positions.push_back(all_positions[i]); }
    num_paired = sorted_paired_positions.size();
    num_total = sorted_all_positions.size();
    for (int i = 0; i < num_paired; i++) {
      deltaE.push_back(0.0);
      switchE.push_back(0.0);
      span_bound.push_back(0);
      block.push_back(0);
      within_filter.push_back(false);
    };
}

void coord_dictionary::get_submatrix_bounds( map_matrix<int> &nmatrix ) {
    for (int i = 0; i < nmatrix.length; i++) { up_bound_submatrix.push_back(i); low_bound_submatrix.push_back(i); }
    for(auto const &ent1 : nmatrix.mat) {
        auto const &i = ent1.first;
        int lowest = i; int highest = i;
        for(auto const &ent2 : ent1.second) {
                auto const &j = ent2.first;
                if (j > highest) { highest = j; }
                if (j < lowest) { lowest = j; }
        }
        up_bound_submatrix[i] = highest;
        low_bound_submatrix[i] = lowest;
    }
    flip_up_bound = up_bound_submatrix;
    flip_low_bound = low_bound_submatrix;
    for (int i = 1; i < nmatrix.length; i++) {
        if (flip_up_bound[i-1] > flip_up_bound[i]) { flip_up_bound[i] = flip_up_bound[i-1]; }
    }
    for (int i = (nmatrix.length-2); i >= 0; i--) {
        if (flip_low_bound[i+1] < flip_low_bound[i]) { flip_low_bound[i] = flip_low_bound[i+1]; }
    }
    for (int i = 0; i < nmatrix.length; i++) { span_bound[i] = flip_up_bound[i] - flip_low_bound[i]; }
}

void coord_dictionary::hap_random_initialization() {
    for (int i = 0; i < num_paired; i++) { haplotype.push_back((int)(rand()%2)*2-1); }
}

void coord_dictionary::hap_zero_initialization() {
    for (int i = 0; i < num_paired; i++) {
        haplotype.push_back(0);
	reload_bool.push_back(false);
    }
}
