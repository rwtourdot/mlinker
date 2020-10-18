#include "block_dict.h"

void block_dictionary::add_bl_het( int block_num, int pos, int hap, std::string ref_string, std::string alt_string ) {
	bl_ref[block_num].push_back(ref_string);
	bl_alt[block_num].push_back(alt_string);
	bl_pos[block_num].push_back(pos);
	bl_hap[block_num].push_back(hap);
	num_hets[block_num] += 1;
	if (num_hets[block_num] == 1) {
		bl_min[block_num] = pos;
		bl_max[block_num] = pos;
	} else {
		if (pos > bl_max[block_num]) { bl_max[block_num] = pos; }
		if (pos < bl_max[block_num]) { bl_min[block_num] = pos; }
	}
};

void block_dictionary::subset_initialize(int num_blocks, int min_hets, int min_links, map_matrix<int> block_matrix) {
	int k = 0;
        for (int i = 0; i < num_blocks; i++) {
		if (num_hets[i] > min_hets) {
			if (bl_nlinks[i] >= min_links) {
				subset_blocks.push_back(i);
				block_map[i] = k;
				block_map_inverted[k] = i;
				k++;
			}
		}
	}
        length_subset = subset_blocks.size();
        map_matrix<int> temp_matrix(length_subset);
        int i = 0;
        for (auto it : subset_blocks) {
                int j = 0;
                for (auto it2 : subset_blocks) {
			if (it != it2) { //temp_matrix.add_to(i,j,block_matrix(it,it2));
				if (block_matrix(it,it2) != 0) {
					cout << it << "\t" << bl_min[it] << "\t" << bl_max[it] << "\t" << it2 << "\t" << bl_min[it2] << "\t" << bl_max[it2] << "\t" << i << "\t" << j << "\t" << block_matrix(it,it2) << "\t" << block_matrix(it2,it) << "\t" << bl_nlinks[it] << "\t" << bl_nlinks[it2] << "\t"<< endl;
				}
				temp_matrix.set_val(i,j,block_matrix(it,it2));
				//temp_matrix.set_val(j,i,block_matrix(it,it2));
			}
                        j++;
                }
                i++;
        }
	for (int i = 0; i < length_subset; i++) {
        	subset_deltaE.push_back(0);
        	subset_switchE.push_back(0);
	}
	subset_matrix = temp_matrix;
};

void block_dictionary::scale_matrix(int scale_limit, map_matrix<int> block_matrix) {
	map_matrix<int> temp_matrix(length_subset);
	int i = 0;
	int num_entries = 0;
	int total_energy = 0;
	cout << "scale limit: " << scale_limit << endl;
        for (auto it : subset_blocks) {
                int j = 0;
                for (auto it2 : subset_blocks) {
			if (it != it2) {
			int diff1 = bl_min[it] - bl_max[it2];
			int diff2 = bl_min[it2] - bl_max[it];
			int diffp = std::max(diff1,diff2);
			//cout << diff1 << "\t" << diff2 << "\t" << diffp << "\t" << block_matrix(it,it2) << "\t" << it << "\t" << it2 << endl ;
			if (diffp < scale_limit) {
				if (block_matrix(it,it2) != 0) {
					num_entries += 1;
					total_energy += std::abs(block_matrix(it,it2));
				}
                        	//temp_matrix.add_to(i,j,block_matrix(it,it2)); //cout << it << "\t" << it2 << "\t" << endl;
				//cout << it << "\t" << it2 << "\t" << diff1 << "\t" << diff2 << "\t" << diffp << "\t" << block_matrix(it,it2) << endl;
                        	temp_matrix.set_val(i,j,block_matrix(it,it2));
                        	//temp_matrix.set_val(j,i,block_matrix(it,it2));
			} }
                        j++;
                }
                i++;
        }
	cout << "num entries: " << num_entries << endl;
	cout << "total energy: " << total_energy << endl;
	distance_matrix = temp_matrix;
}

void block_dictionary::subset_hap_random_initialization() {
	//for (int i = 0; i < length_subset; i++) { subset_haplotype.push_back((int)(rand()%2)*2-1); }
	for (int i = 0; i < length_subset; i++) { subset_haplotype.push_back((int)(1)); }
};
