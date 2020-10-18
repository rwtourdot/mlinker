#include "mc_solver.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////
void solver( std::unordered_map<std::string,variant_node>& var_dict, coord_dictionary& pdict, map_matrix<double> diff_matrix, map_matrix<int> num_matrix_second ) {
	float t = clock();
	static const std::size_t length = pdict.num_paired;
	map_matrix<int> num_matrix_third(length);
	//cout << "========" << endl;
	cout << "monte carlo solver: " << endl;
	length_cutoff_nmatrix( pdict, diff_matrix, num_matrix_second, num_matrix_third );
	vector<bool> flip_maxima(length);
	for (int j = 0; j < pdict.num_paired; j++) { flip_maxima[j] = true; }
	for (int i = 0; i < solver_loops; i++) {
		cout << "******************* solver loop     " << i << endl;
		t = clock();
		single_spin_flip_map( pdict, diff_matrix, num_matrix_third );  // num_matrix_second
		if (i > 0) { get_flip_positions( pdict, flip_maxima ); }
		t = float(clock() - t)/ CLOCKS_PER_SEC;
		cout << "spin flip --------- time: " << t << endl;
		t = clock();
		cout << "-- finished flip spins     " << endl;
		block_flip_map( pdict, diff_matrix, num_matrix_third, flip_maxima );  // num_matrix_second
		t = float(clock() - t)/ CLOCKS_PER_SEC;
		cout << "block flip -------- time: " << t << endl;
		cout << "-- finished block flip     " << endl;
		energy_sum_min_max( pdict );
	}
	set_switch_energy( pdict, flip_maxima );   /// create a new haplotype block class and find the segments
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
void solver_recursive( std::unordered_map<std::string,variant_node>& var_dict, coord_dictionary& pdict, map_matrix<double> diff_matrix, map_matrix<int> num_matrix_second ) {
	float t = clock();
	static const std::size_t length = pdict.num_paired;
	map_matrix<int> num_matrix_third(length);
	//cout << "========" << endl;
	cout << "monte carlo solver: " << endl;
	length_cutoff_nmatrix( pdict, diff_matrix, num_matrix_second, num_matrix_third );
	for (int i = 0; i < solver_loops; i++) {
		cout << "******************* solver loop     " << i << endl;
		t = clock();
		single_spin_flip_recursive( pdict, diff_matrix, num_matrix_third );
		t = float(clock() - t)/ CLOCKS_PER_SEC;
		cout << "spin flip --------- time: " << t << endl;
		t = clock();
		block_flip_recursive( pdict, diff_matrix, num_matrix_third );
		t = float(clock() - t)/ CLOCKS_PER_SEC;
		cout << "block flip -------- time: " << t << endl;
		energy_sum_min_max( pdict );
	}
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
void solver_recursive_hic( block_dictionary& bdict, std::vector<int> hic_limit_loop, map_matrix<int> block_matrix ) {
	float t = clock();
	//cout << "========" << endl;
	cout << "monte carlo solver:" << endl;
	for( int k = 0; k < hic_limit_loop.size(); k++ ) {
		int hlimit = hic_limit_loop[k];
		bdict.scale_matrix( hlimit, block_matrix );
		for (int i = 0; i < solver_loops_hic; i++) {
			cout << "******************* solver_loop" << "\t" << i << "\tscale_loop\t" << hlimit << endl;
			t = clock();
			single_spin_flip_recursive_hic( bdict );
			//block_phasing_quick_flip( bdict, graph_flip_list );
			t = float(clock() - t)/ CLOCKS_PER_SEC;
			cout << "spin flip --------- time: " << t << endl;
			t = clock();
			block_flip_brute_force_hic( bdict );  //block_flip_recursive_hic( bdict );
			t = float( clock() - t )/ CLOCKS_PER_SEC;
			cout << "block flip -------- time: " << t << endl;
			energy_sum_min_max_hic( bdict );
		}
	}
};

static void block_phasing_quick_flip( block_dictionary& bdict, vector<int>& graph_flip_list ) {
	//
	vector<int> loop_haplotype = bdict.subset_haplotype;
	std::unordered_map<int,std::unordered_map<int,int>> correlation_graph;
	std::unordered_map<int,std::unordered_map<int,int>> high_correlation_graph;
	for (int l = 0; l < bdict.length_subset; l++) {
		for ( auto const &ent1 : bdict.subset_matrix.mat[l] ) {
			auto const &m = ent1.first;
			if ( bdict.subset_matrix(l,m) != 0 & l != m ) {
			int lm_energy = loop_haplotype[l]*loop_haplotype[m]*bdict.subset_matrix(l,m);
			for ( auto const &ent2 : bdict.subset_matrix.mat[m] ) {
				auto const &p = ent2.first;
				if ( bdict.subset_matrix(l,p) != 0 & l != p & m != p ) {
				int lp_energy = loop_haplotype[l]*loop_haplotype[p]*bdict.subset_matrix(l,p);
				int correlation = lm_energy*lp_energy;
				//cout << correlation << endl;
				if ( correlation > 0 ) {
					//cout << l << "\t" << m << "\t" << p << "\t" << correlation << endl;
					correlation_graph[m][p] += 1;
					correlation_graph[p][m] += 1; //.push_back(m);
					if (correlation_graph[m][p] > 30) {  //2-
						high_correlation_graph[m][p] += 1;
						high_correlation_graph[p][m] += 1;
					}
				} }
			} }
                }
	}

	vector<int> covered_list = {};
	vector<int> traverse_list;
	for (int l = 0; l < bdict.length_subset; l++) {
		if (std::find(covered_list.begin(), covered_list.end(), l) != covered_list.end()) {
			cout << "covered" << endl;
		} else {
        		traverse_graph( high_correlation_graph, traverse_list, l, bdict.length_subset );
			cout << "start " << l << endl;
			for (int l = 0; l < traverse_list.size(); l++) {
				cout << traverse_list[l] << " ";
				covered_list.push_back(traverse_list[l]);
			}
			cout << endl;
		}
	}

	for (int l = 0; l < traverse_list.size(); l++) {
		cout << traverse_list[l] << " ";
		bdict.subset_haplotype[l] = loop_haplotype[l]*-1;
	}
	cout << endl;
	graph_flip_list = traverse_list;
}

static void traverse_graph( std::unordered_map<int,std::unordered_map<int,int>> high_correlation_graph, vector<int>& traverse_list, int start_node, int subset_length ) {
	traverse_list.clear();
        bool *visited = new bool[subset_length];
        for(int i = 0; i < subset_length; i++) { visited[i] = false; }
        // Create a queue for BFS
        int s = start_node;
	visited[s] = true;
        std::list<int> queue;
        queue.push_back(s);
	traverse_list.push_back(s);

        /////////////////////////////////////////////////////////
        while(!queue.empty()) {
                s = queue.front();  //std::cout << s << " ";
                queue.pop_front();
                for (auto& it : high_correlation_graph[s]) {
                        int j = it.first;
                        if( !visited[j] ) {
				//cout << high_correlation_graph[s][j] << endl;
                                visited[j] = true;
                                queue.push_back(j);
                                traverse_list.push_back(j);
                        }
                }
        }  //cout << endl;
        std::sort( traverse_list.begin(),traverse_list.end() );
}



////////////////////////////////////////////////////////////////////////////////////////////////////////
static void length_cutoff_nmatrix( coord_dictionary& pdict, map_matrix<double>& diff_matrix, map_matrix<int>& nmatrix, map_matrix<int>& nmatrix2 ) {
        for (int i = 0; i < pdict.num_paired; i++) {
                for (auto const &ent1 : nmatrix.mat[i]) {
                        auto const &m = ent1.first;      //cout << i << " " << m << " " << diff_pos << " " << nmatrix(i,m) << endl;
                        int diff_pos = abs(pdict.sorted_paired_positions[m] - pdict.sorted_paired_positions[i]);
                        if (diff_pos < pos_diff_cutoff) { nmatrix2.add_to(i,m,nmatrix(i,m)); }
                }
        }
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
void call_blocks(coord_dictionary& pdict,int switch_cutoff ) {
	int block_num = 0;
        for (int i = 0; i < pdict.num_paired; i++) {
		if (pdict.switchE[i] > switch_cutoff) { block_num += 1; }
		cout << i << " " << pdict.switchE[i] << " " << block_num << endl;
		pdict.block[i] = block_num;
	}
	//cin.get();
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
static void set_switch_energy( coord_dictionary& pdict, vector<bool>& flip_maxima ) {
        for (int i = 0; i < pdict.num_paired; i++) {
                if (!flip_maxima[i]) { pdict.switchE[i] = -10000.0; }
        }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
static void get_flip_positions( coord_dictionary& pdict, vector<bool>& flip_maxima ) {
        for (int i = 0; i < pdict.num_paired; i++) {
                vector<int> neighbors{i-1,i+1};  //cout << "  " << i << "  ";
                double energy = pdict.deltaE[i];  bool check = true;
		if ( energy < -200 ) {
                for (int j=0; j < neighbors.size(); j++) {   //cout << neighbors[j] << "  ";
                        if ( neighbors[j] >= 0 && neighbors[j] < pdict.num_paired ) {
				if ( pdict.deltaE[neighbors[j]] > energy ) { check = false; }
			}
                } }
                flip_maxima[i] = check;
        }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
static void energy_sum_min_max( coord_dictionary& pdict ) {
        double sum = 0.0; double max = 0.0; double min = 0.0; int az_num = 0;
        double sum_blockE = 0.0; double max_blockE = 0.0; double min_blockE = 0.0; int az_num_blockE = 0;
        for (int i = 0; i < pdict.num_paired; i++) {
                sum += pdict.deltaE[i];
                sum_blockE += pdict.switchE[i];
                if (pdict.deltaE[i] > max) { max = pdict.deltaE[i]; }
                if (pdict.deltaE[i] < min) { min = pdict.deltaE[i]; }
                if (pdict.deltaE[i] > 0) { az_num += 1; }
                if (pdict.switchE[i] > max_blockE) { max_blockE = pdict.switchE[i]; }
                if (pdict.switchE[i] < min_blockE) { min_blockE = pdict.switchE[i]; }
                if (pdict.switchE[i] > 0) { az_num_blockE += 1; }
		//cout << i << " " << pdict.haplotype[i] << " " << pdict.deltaE[i] << " " << pdict.switchE[i] << endl;
        }
	cout << "-- num paired: " << pdict.num_paired << endl;
        cout << "-- spin flip energy sum: " << sum << " minimum: " << min << " maximum: " << max << " #_above_zero: " << az_num << endl;
        cout << "-- block flip energy sum: " << sum_blockE << " minimum: " << min_blockE << " maximum: " << max_blockE << " #_above_zero: " << az_num_blockE << endl;
};

// added negative signs
////////////////////////////////////////////////////////////////////////////////////////////////////////
static void single_spin_flip_recursive( coord_dictionary& pdict, map_matrix<double>& diff_matrix, map_matrix<int>& nmatrix ) {
        int nsf = 0;
        for (int i = 0; i < pdict.num_paired; i++) {
                int init_i = pdict.haplotype[i]; //int final_i = -1*pdict.haplotype[i];
                double initial_energy = 0.0;     //double final_energy = 0.0;
                for (auto const &ent1 : nmatrix.mat[i]) {
                        auto const &m = ent1.first;
                        int diff_pos = abs(pdict.sorted_paired_positions[m] - pdict.sorted_paired_positions[i]);
                        initial_energy += -1.0*init_i*pdict.haplotype[m]*diff_matrix(i,m);  //final_energy += final_i*pdict.haplotype[m]*diff_matrix(i,m);
                        //cout << i << " " << m << " " << pdict.haplotype[i] << " " << diff_pos << " " << nmatrix(i,m) << " " << diff_matrix(i,m) << endl;
                }
                double final_energy = -1*initial_energy;
                double diff_energy = final_energy - initial_energy;
                //cout << "init e " << initial_energy << " final e " << final_energy << " e diff " << diff_energy << endl;
                if (diff_energy < 0.0) { pdict.haplotype[i] = pdict.haplotype[i]*(-1); nsf++; }
                pdict.deltaE[i] = (-1.0)*diff_energy;
        }
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
static void block_flip_recursive( coord_dictionary& pdict, map_matrix<double>& diff_matrix, map_matrix<int>& nmatrix ) {
        int nbf = 0;
        double initial_energy = 0.0;
        for (int i = 0; i < pdict.num_paired; i++) {
                for (auto const &ent1 : nmatrix.mat[i]) {
                        auto const &m = ent1.first;
                        initial_energy += -1.0*pdict.haplotype[i]*pdict.haplotype[m]*diff_matrix(i,m);
                }
        }
        double energyA = initial_energy;
        double energyB = initial_energy;
        vector<int> loop_haplotype = pdict.haplotype;
        for (int i = 0; i < pdict.num_paired; i++) {
                double save_energyB = energyB;
                vector<int> switch_haplotype = loop_haplotype;
                for (int j = 0; j < i+1; j++) { switch_haplotype[j] = loop_haplotype[j]*(-1); }
                int hap_spin  = loop_haplotype[i];
                int flip_spin = -1*loop_haplotype[i];
                double spin_i_energy = 0.0;
                double spin_f_energy = 0.0;
                for (auto const &ent1 : nmatrix.mat[i]) {
                        auto const &m = ent1.first;
                        spin_i_energy += -1.0*hap_spin*switch_haplotype[m]*diff_matrix(i,m);  /// may need to add negative sign
                        spin_f_energy += -1.0*flip_spin*switch_haplotype[m]*diff_matrix(i,m);
                }
                double switch_energy = spin_f_energy - spin_i_energy;
                energyA = energyA + switch_energy;
                double diff_energy = energyA - energyB;
                if ( diff_energy < 0.0 ) {
                        loop_haplotype = switch_haplotype;
                        energyB = energyA;
                        energyA = save_energyB;
                        nbf++;
                }
                else {
                        energyB = save_energyB;
                        energyA = energyA;
                }
                pdict.switchE[i] = (-1.0)*diff_energy;
        }
        pdict.haplotype = loop_haplotype;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
static void energy_sum_min_max_hic( block_dictionary& bdict ) {
        double sum = 0.0; double max = 0.0; double min = 0.0;
        for (int i = 0; i < bdict.length_subset; i++) {
                sum += bdict.subset_deltaE[i];
                if (bdict.subset_deltaE[i] > max) { max = bdict.subset_deltaE[i]; }
                if (bdict.subset_deltaE[i] < min) { min = bdict.subset_deltaE[i]; }
        }
        cout << "-- energy sum: " << sum << " minimum: " << min << " maximum: " << max << endl;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
static void single_spin_flip_recursive_hic( block_dictionary& bdict ) {
        int nsf = 0;
        for (int i = 0; i < bdict.length_subset; i++) {
                int init_i = bdict.subset_haplotype[i]; //int final_i = -1*pdict.haplotype[i];
                double initial_energy = 0.0;     //double final_energy = 0.0;
                for (auto const &ent1 : bdict.distance_matrix.mat[i]) {
                        auto const &m = ent1.first;
                        //int diff_pos = abs(pdict.sorted_paired_positions[m] - pdict.sorted_paired_positions[i]);
                        initial_energy += -1.0*init_i*bdict.subset_haplotype[m]*bdict.distance_matrix(i,m);  //final_energy += final_i*pdict.haplotype[m]*diff_matrix(i,m);
                }
                double final_energy = -1*initial_energy;
                double diff_energy = final_energy - initial_energy;
                //cout << "init e " << initial_energy << " final e " << final_energy << " e diff " << diff_energy << endl;
                if (diff_energy < 0.0) { bdict.subset_haplotype[i] = bdict.subset_haplotype[i]*(-1); nsf++; }
                bdict.subset_deltaE[i] = (-1.0)*diff_energy;
        }
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
static void block_flip_recursive_hic( block_dictionary& bdict ) {
        int nbf = 0;
        double initial_energy = 0.0;
        for (int i = 0; i < bdict.length_subset; i++) {
                for (auto const &ent1 : bdict.distance_matrix.mat[i]) {
                        auto const &m = ent1.first;
                        initial_energy += -1.0*bdict.subset_haplotype[i]*bdict.subset_haplotype[m]*bdict.distance_matrix(i,m);
                }
        }
        double energyA = initial_energy;
        double energyB = initial_energy;
        //vector<int> loop_haplotype = pdict.haplotype;
        vector<int> loop_haplotype = bdict.subset_haplotype;
        for (int i = 0; i < bdict.length_subset; i++) {
                double save_energyB = energyB;
                vector<int> switch_haplotype = loop_haplotype;
                for (int j = 0; j < i+1; j++) { switch_haplotype[j] = loop_haplotype[j]*(-1); }
                int hap_spin  = loop_haplotype[i];
                int flip_spin = -1*loop_haplotype[i];
                double spin_i_energy = 0.0;
                double spin_f_energy = 0.0;
                for (auto const &ent1 : bdict.distance_matrix.mat[i]) {
                        auto const &m = ent1.first;
                        spin_i_energy += -1.0*hap_spin*switch_haplotype[m]*bdict.distance_matrix(i,m);  /// may need to add negative sign
                        spin_f_energy += -1.0*flip_spin*switch_haplotype[m]*bdict.distance_matrix(i,m);
                }
                double switch_energy = spin_f_energy - spin_i_energy;
                energyA = energyA + switch_energy;
                double diff_energy = energyA - energyB;
                if ( diff_energy < 0.0 ) {
                        loop_haplotype = switch_haplotype;
                        energyB = energyA;
                        energyA = save_energyB;
                        nbf++;
                }
                else {
                        energyB = save_energyB;
                        energyA = energyA;
                }
                bdict.subset_switchE[i] = (-1.0)*diff_energy;
        }
        bdict.subset_haplotype = loop_haplotype;
};


////////////////////////////////////////////////////////////////////////////////////////////////////////
static void single_spin_flip_map( coord_dictionary& pdict, map_matrix<double>& diff_matrix, map_matrix<int>& nmatrix ) {
        int nsf = 0; clock_t t;
        for (int i = 0; i < pdict.num_paired; i++) {
                vector<int> sparse_map,sparse_hap,init_sub_haplotype,final_sub_haplotype;
                int l = 0; bool cross_i = true; int spin_pos;  double init_energy,final_energy,diff_energy;
                for(auto const &ent1 : nmatrix.mat[i]) {
                        auto const &m = ent1.first;
                        if( m>i && cross_i) { sparse_map.push_back(i); cross_i = false; spin_pos = l; l++; }
                        int diff_pos = abs(pdict.sorted_paired_positions[m] - pdict.sorted_paired_positions[i]);   ////// CHANGED THIS
                        //int diff_pos2 = abs(pdict.double_positions[m] - pdict.double_positions[i]);
                        //cout << "i " << i << "\tm " << m << "\t" << pdict.sorted_paired_positions[i] << "\t" << pdict.sorted_paired_positions[m] << "\t  sorted_paired " << diff_pos << "\t" << pdict.double_positions[i] << "\t" << pdict.double_positions[m] << "\t double_pos " << diff_pos2 << endl;      //////////////////////////
                        if (diff_pos < pos_diff_cutoff) { sparse_map.push_back(m); l++;}
                }
                if(cross_i) { sparse_map.push_back(i); cross_i = false; spin_pos = l; }
                for (int j = 0; j < sparse_map.size(); j++) { sparse_hap.push_back(pdict.haplotype[sparse_map[j]]); }
                submatrix_opt tempsub;  tempsub.initialize(sparse_map,diff_matrix);
                init_sub_haplotype = sparse_hap;
                final_sub_haplotype = init_sub_haplotype;
                //cout << "flip  position: " << i << " sparse_map_size " << sparse_map.size() << " spin_pos " << spin_pos << "  total " << l << endl;
                final_sub_haplotype[spin_pos] = final_sub_haplotype[spin_pos]*(-1);
                init_energy  = energy_function_opt(tempsub,init_sub_haplotype);
                final_energy = energy_function_opt(tempsub,final_sub_haplotype);
                diff_energy  = final_energy - init_energy;        //cout << "init e " << init_energy << " final e " << final_energy << " e diff " << diff_energy << endl;
                if (diff_energy < 0.0) { pdict.haplotype[i] = pdict.haplotype[i]*(-1); nsf++; }
                pdict.deltaE[i] = (-1.0)*diff_energy;
        }
        cout << "-- number of spin flips: " << nsf << endl;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
static void block_flip_map( coord_dictionary& pdict, map_matrix<double>& diff_matrix, map_matrix<int>& nmatrix, vector<bool>& flip_maxima ) {
        int nbf = 0;
        for (int i = 0; i < pdict.num_paired; i++) {
                if (flip_maxima[i]) {                   //if (pdict.switchE[i] > -switch_cut) {
                vector<int> sparse_map,sparse_hap,init_sub_haplotype,final_sub_haplotype;
                int flip_pos,l=0;   double init_energy,final_energy,diff_energy;
                for ( int j = pdict.flip_low_bound[i]; j <= pdict.flip_up_bound[i]; j++ ) {
                        if (j == i) { flip_pos = l; }
                        int diff_pos = abs(pdict.sorted_paired_positions[j] - pdict.sorted_paired_positions[i]);   ////// CHANGED THIS
                        //int diff_pos2 = abs(pdict.double_positions[j] - pdict.double_positions[i]);
                        //cout << "i " << i << "\tm " << j << "\t" << pdict.sorted_paired_positions[i] << "\t" << pdict.sorted_paired_positions[j] << "\t  sorted_paired " << diff_pos << "\t" << pdict.double_positions[i] << "\t" << pdict.double_positions[j] << "\t double_pos " << diff_pos2 << endl;
                        if (diff_pos < pos_diff_cutoff) { sparse_map.push_back(j); l++; }
                }
                for (int j = 0; j < sparse_map.size(); j++) { sparse_hap.push_back(pdict.haplotype[sparse_map[j]]); }
                //cout << "block position: " << i << " sparse_map_size " << sparse_map.size() << " flip_pos " << flip_pos << "  total " << l << endl;
                submatrix_opt tempsub; tempsub.initialize(sparse_map,diff_matrix);
                init_sub_haplotype = sparse_hap;
                final_sub_haplotype = init_sub_haplotype;
                for (int j = 0; j < flip_pos; j++) { final_sub_haplotype[j] = final_sub_haplotype[j]*(-1); }
                init_energy  = energy_function_opt(tempsub,init_sub_haplotype);
                final_energy = energy_function_opt(tempsub,final_sub_haplotype);
                diff_energy  = final_energy - init_energy;        //cout << "init e " << init_energy << " final e " << final_energy << " e diff " << diff_energy << endl;
                if (diff_energy < 0.0) { for (int l = 0; l < i; l++) { pdict.haplotype[l] = pdict.haplotype[l]*(-1); } nbf++; }
                pdict.switchE[i] = (-1.0)*diff_energy;
                }    //}
        }
        cout << "-- number of block flips: " << nbf << endl;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
static void block_flip_brute_force_hic( block_dictionary& bdict ) {
        int nbf = 0;
        vector<int> loop_haplotype = bdict.subset_haplotype;
        for (int i = 0; i < bdict.length_subset; i++) {
                vector<int> switch_haplotype = loop_haplotype;
                for (int j = 0; j < i+1; j++) { switch_haplotype[j] = loop_haplotype[j]*(-1); }
                double spin_i_energy = 0.0;
		double spin_f_energy = 0.0;
		for (int l = 0; l < bdict.length_subset; l++) {
                	for (auto const &ent1 : bdict.distance_matrix.mat[l]) {
                        	auto const &m = ent1.first;
                        	spin_i_energy += -1.0*loop_haplotype[l]*loop_haplotype[m]*bdict.distance_matrix(l,m);  /// may need to add negative sign
                        	spin_f_energy += -1.0*switch_haplotype[l]*switch_haplotype[m]*bdict.distance_matrix(l,m);
				//cout << l << "\t" << m << "\t"<<  bdict.distance_matrix(l,m) << "\t" << bdict.distance_matrix(m,l) << endl;
                	}
		}
                double switch_energy = spin_f_energy - spin_i_energy;
		//cout << "switch energy:\t" << i << "\t" << bdict.block_map_inverted[i] << "\t" << bdict.num_hets[bdict.block_map_inverted[i]] << "\t" << spin_i_energy << "\t" << spin_f_energy << "\t" << switch_energy << endl;
		bool flipped = false;
                if ( switch_energy < 0.0 ) {
                        loop_haplotype = switch_haplotype;
			flipped = true;
                        nbf++;
                }
		if ( flipped == true ) { bdict.subset_switchE[i] = switch_energy; }
		else { bdict.subset_switchE[i] = (-1.0)*switch_energy; }
        }
        bdict.subset_haplotype = loop_haplotype;
        cout << "-- number of block flips: " << nbf << endl;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
static double energy_function_opt( submatrix_opt tempsub, vector<int> haplotype ) {
        vector<double> temp_vector = dot_product_matrix(haplotype,tempsub.smat);
        double energy_return = (-0.5)*dot_product(haplotype,temp_vector);
        return energy_return;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
static double dot_product( vector<int> a, vector<double> b ) {
        double dot = 0.0; for (int i = 0; i < a.size(); i++) { dot += (double)(a[i])*(b[i]); };
        return dot;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
static vector<double> dot_product_matrix( vector<int> a, std::vector< std::vector<double> > b ) {
        vector<double> dot_vector;
        for (int i = 0; i < a.size(); i++) {
                double i_value = 0.0; for (int j = 0; j < b.size(); j++) { i_value += (double)a[j]*b[i][j]; }
                dot_vector.push_back(i_value);
        }
        return dot_vector;
};
