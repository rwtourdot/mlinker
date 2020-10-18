#include "build_hap_matrix.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////
void link_hashes( std::unordered_map<std::string,variant_node>& var_dict, std::unordered_map<std::string,read_tree>& read_graph ) {
        //cout << "========" << endl;
        int i = 0;   clock_t t;
        t = clock();
        for (auto& it : var_dict) {  //int j = 0;
		it.second.make_unique();
                for (int l=0; l < it.second.connected_reads_long_form.size(); l++) {
                        for (int m=0; m < read_graph[it.second.connected_reads_long_form[l]].connected_strings.size(); m++) {
                                std::string connected_het = read_graph[it.second.connected_reads_long_form[l]].connected_strings[m];
				//cout << i << "\t" << it.second.connected_reads_long_form[l] << "\t" << endl;
                                if (connected_het != it.first) { it.second.add_connection(connected_het); }
                        }
                }
                i++;
        }
        for (auto& it : var_dict) { it.second.count_connections(); }
        t = clock() - t;
        //cout << "connected hashes ----- time: " << t << endl;
        return;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
void calc_coverage_unique_hash( std::unordered_map<std::string,variant_node>& var_dict ) {
	for (auto& it : var_dict) { it.second.unique_hash(); }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
void prune_graph( std::unordered_map<std::string,variant_node>& var_dict ) {
        std::vector<int> het_pos;
        std::vector<int> hom_pos;
        std::vector<int> pos_intersection,pos_union;
        int j = 0;
        for (auto& it : var_dict) {
                if (it.second.num_hets > 0) {
                        if (it.second.var) { het_pos.push_back(it.second.pos); }
                        else { hom_pos.push_back(it.second.pos); }
                }
                j += 1;
        }  //int p = j/2;
        std::set<int> set_het(het_pos.begin(),het_pos.end());
        std::set<int> set_hom(hom_pos.begin(),hom_pos.end());
        std::set_intersection(set_het.begin(),set_het.end(),set_hom.begin(),set_hom.end(),std::back_inserter(pos_intersection));
        std::set_union(set_het.begin(),set_het.end(),set_hom.begin(),set_hom.end(),std::back_inserter(pos_union));
        cout << "number of variants: " << j << endl;    //cout << "paired_variants:    " << p << endl;
        cout << "var intersection:   " << pos_intersection.size() << endl;
        cout << "var union:          " << pos_union.size() << endl;
        for (auto& it : var_dict) {
                bool found = (std::find(pos_intersection.begin(), pos_intersection.end(), it.second.pos) != pos_intersection.end());
                if (found) { it.second.paired = true; }
                else { it.second.paired = false; }
        }
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
void calculate_epso( int i, int j, map_matrix_vector& expanded, map_matrix<int>& nmatrix, map_matrix<int>& span_matrix, map_matrix<double>& corr_matrix, map_matrix<double>& diff_matrix, int &ecount, double &esum ) {
        int vv=0,rv=0,total=0,span=0;
        double diff = 0.0,frac=0.0,corr=0.0;
        for (int m = 0; m < nmatrix(i,j); m++) {
                total++;
                if (expanded.return_val(i,j,m)) { vv++; } else { rv++; }
        }
        //if (i == j) { cout << " within the loop " << total << "  " << vv << "  " << rv << endl; }
        if (total > 0 ) {
                if (vv > rv) { span = rv; }
                if (rv >= vv) { span = vv; }
                diff = (double)(vv-rv)*abs(vv-rv)/(vv+rv);
                corr = (double)span/total;
				/*
                diff_matrix.set_val(i,j,corr);
                corr_matrix.set_val(i,j,corr);
                span_matrix.set_val(i,j,span);
				*/
				esum+=corr;
				ecount++;
        };

};

////////////////////////////////////////////////////////////////////////////////////////////////////////
//void calculate_link_fractions( int i, int j, map_matrix_vector& expanded, map_matrix<int>& nmatrix, map_matrix<int>& span_matrix, map_matrix<double>& corr_matrix, map_matrix<double>& diff_matrix ) {
//        int vv=0,rv=0,total=0,span=0;
//        double diff = 0.0,frac=0.0,corr=0.0;
//        for (int m = 0; m < nmatrix(i,j); m++) {
//                total++;
//                if (expanded.return_val(i,j,m)) { vv++; } else { rv++; }
//        }
        //if (i == j) { cout << " within the loop " << total << "  " << vv << "  " << rv << endl; }
//        if (total > 0 ) {
//                if (vv > rv) { span = vv; }
//                if (rv >= vv) { span = rv; }
//                diff = (double)(vv-rv)*abs(vv-rv)/(vv+rv);
//                corr = (double)(vv-rv)/total;
//                diff_matrix.set_val(i,j,diff);
//                corr_matrix.set_val(i,j,corr);
//                span_matrix.set_val(i,j,span);
//        };
//};

////////////////////////////////////////////////////////////////////////////////////////////////////////
//void link_matrix_calculations( map_matrix_vector& expanded, map_matrix<int>& nmatrix, map_matrix<int>& span_matrix, map_matrix<double>& corr_matrix, map_matrix<double>& diff_matrix ) {
//        for(auto const &ent1 : nmatrix.mat) {
//                auto const &i = ent1.first;
//                auto const &inner_map = ent1.second;
//                for(auto const &ent2 : ent1.second) {
//                        auto const &j = ent2.first;
//                        auto const &inner_value = ent2.second;   //cout << i << "  " << j << endl;
//                        calculate_link_fractions(i,j,expanded,nmatrix,span_matrix,corr_matrix,diff_matrix);
//                }
//        }
//};

void calculate_link_fractions( int i, int j, map_matrix_vector& expanded, map_matrix<int>& nmatrix, map_matrix<int>& span_matrix, map_matrix<double>& corr_matrix, map_matrix<double>& diff_matrix, double epso ) {
        int vv=0,rv=0,total=0,span=0;
        double diff = 0.0,frac=0.0,corr=0.0,eps=0.0;
        for (int m = 0; m < nmatrix(i,j); m++) {
                total++;
                if (expanded.return_val(i,j,m)) { vv++; } else { rv++; }
        }
        //if (i == j) { cout << " within the loop " << total << "  " << vv << "  " << rv << endl; }
        if (total > 0 ) {
                if (vv > rv) { span = rv; }
                if (rv >= vv) { span = vv; }
                eps = (double)span/total;
				if (epso > eps) { eps = epso; }
				//if (eps < 0 ) {cout << "eps under zero!" << endl; }
				corr = log((1-eps)/eps)*(vv-rv);

                diff_matrix.set_val(i,j,corr);
                corr_matrix.set_val(i,j,corr);
                span_matrix.set_val(i,j,span);
        };

};

////////////////////////////////////////////////////////////////////////////////////////////////////////
void link_matrix_calculations( map_matrix_vector& expanded, map_matrix<int>& nmatrix, map_matrix<int>& span_matrix, map_matrix<double>& corr_matrix, map_matrix<double>& diff_matrix ) {
	int ecount=0;
	double esum=0.0;
        for(auto const &ent1 : nmatrix.mat) {
                auto const &i = ent1.first;
                auto const &inner_map = ent1.second;
                for(auto const &ent2 : ent1.second) {
                        auto const &j = ent2.first;
                        auto const &inner_value = ent2.second;   //cout << i << "  " << j << endl;
						calculate_epso(i,j,expanded,nmatrix,span_matrix,corr_matrix,diff_matrix,ecount,esum);
                }
        }

		double epso = esum/ecount;
		cout << "epso is: " << epso << endl;
		cout << "ecount is: " << ecount << endl;
		cout << log((1-epso)/epso)  << endl;

		for(auto const &ent1 : nmatrix.mat) {
                auto const &i = ent1.first;
                auto const &inner_map = ent1.second;
                for(auto const &ent2 : ent1.second) {
                        auto const &j = ent2.first;
                        auto const &inner_value = ent2.second;   //cout << i << "  " << j << endl;
                        calculate_link_fractions(i,j,expanded,nmatrix,span_matrix,corr_matrix,diff_matrix,epso);
                }
        }
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
ptrdiff_t get_index_var( vector<int> paired_pos, int find_p ) { return std::find(paired_pos.begin(),paired_pos.end(),find_p) - paired_pos.begin(); };

////////////////////////////////////////////////////////////////////////////////////////////////////////
void initialize_pdict( std::unordered_map<std::string,variant_node>& var_dict, coord_dictionary& pdict, bool loh_mode ) { pdict.initialize(var_dict,loh_mode); }

////////////////////////////////////////////////////////////////////////////////////////////////////////
void initialize_solver( std::unordered_map<std::string,variant_node>& var_dict, coord_dictionary& pdict, map_matrix<double>& diff_matrix, map_matrix<int>& num_matrix_second ) {
        clock_t t;   t = clock();
        static const std::size_t length = pdict.num_paired;
        // added this loop to presearch pdict for index positions
        std::unordered_map<int,int> index_map;
        for (auto& it : var_dict) {
          int vpos = it.second.pos;
          ptrdiff_t i = get_index_var( pdict.sorted_paired_positions, vpos );
          int index = i;
          if (index < pdict.sorted_paired_positions.size() ) { index_map[vpos] = index; }
          else { index_map[vpos] = -1; }
        }
        // end of added code
        map_matrix<int> num_matrix(length),span_matrix(length);
        map_matrix<double> corr_matrix(length),binomial_matrix(length); //diff_matrix(length),
        map_matrix_vector link_matrix(length);
        //cout << "========" << endl;
        cout << "initializing solver: " << endl;
        int i = 0;
        for (auto& it : pdict.sorted_paired_positions) {
                std::string node1 = pdict.paired_dict[it].at(0);
                std::string node2 = pdict.paired_dict[it].at(1);  //number_of_connections = var_dict[node1].num_hets + var_dict[node2].num_hets
                for (auto it2 : var_dict[node1].connections) {
                        //ptrdiff_t j = get_index_var( pdict.sorted_paired_positions,var_dict[it2.first].pos );    //cout << i << "  " << j << endl;
                        int j = index_map[var_dict[it2.first].pos];
                        if (j < pdict.sorted_paired_positions.size() && j >= 0 && i != j) {
                                link_matrix.add(i,j,var_dict[node1].var,var_dict[it2.first].var,it2.second);
                                num_matrix.add_to(i,j,it2.second);
                        }
                }
                for (auto it2 : var_dict[node2].connections) {
                        //ptrdiff_t j = get_index_var( pdict.sorted_paired_positions,var_dict[it2.first].pos );
                        int j = index_map[var_dict[it2.first].pos];
                        if (j < pdict.sorted_paired_positions.size() && j >= 0 && i != j) {
                                link_matrix.add(i,j,var_dict[node2].var,var_dict[it2.first].var,it2.second);
                                num_matrix.add_to(i,j,it2.second);
                        }
                }
                i++;
        }
        //////////////// create a second matrix which only includes connections with at least a certain weight (number of links)
        for (int i = 0; i < pdict.num_paired; i++) {
                for(auto const &ent1 : num_matrix.mat[i]) {
                        if (ent1.second >= minimum_link_number) { num_matrix_second.set_val(i,ent1.first,ent1.second); }
                }
        }
        cout << "size of num matrix: " << diff_matrix.mat.size() << endl;
        //t = clock() - t;
        //cout << "created link matrix -- time: " << t << endl;
        link_matrix_calculations(link_matrix,num_matrix_second,span_matrix,corr_matrix,diff_matrix);
        pdict.get_submatrix_bounds(num_matrix_second);
        cout << "double positions size: " << pdict.double_positions.size() << endl;
        cout << "paired positions size: " << pdict.sorted_paired_positions.size() << endl;
        cout << "size of num matrix: " << diff_matrix.mat.size() << endl;
        pdict.hap_random_initialization();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
void initialize_solver_loh( std::unordered_map<std::string,variant_node>& var_dict, coord_dictionary& pdict, map_matrix<double>& diff_matrix, map_matrix<int>& num_matrix_second ) {
        clock_t t;   t = clock();
        static const std::size_t length = pdict.num_paired;
        // added this loop to presearch pdict for index positions
        std::unordered_map<int,int> index_map;
        for (auto& it : var_dict) {
          int vpos = it.second.pos;
          ptrdiff_t i = get_index_var( pdict.sorted_paired_positions, vpos );
          int index = i;
          if (index < pdict.sorted_paired_positions.size() ) { index_map[vpos] = index; }
          else { index_map[vpos] = -1; }
        }
        // end of added code
        cout << length << endl;
        map_matrix<int> num_matrix(length),span_matrix(length);
        map_matrix<double> corr_matrix(length),binomial_matrix(length); //diff_matrix(length),
        map_matrix_vector link_matrix(length);
        cout << "========" << endl;
        cout << "initializing solver: " << endl;
        int i = 0;
        for (auto& it : pdict.sorted_paired_positions) { //if(it == 10558150) {
          cout << it << "  " << pdict.paired_dict[it].size() << endl;
          for (int k = 0; k < pdict.paired_dict[it].size(); k++) {
                  std::string node_loop = pdict.paired_dict[it].at(k);
                  for (auto it2 : var_dict[node_loop].connections) {
                          //ptrdiff_t j = get_index_var( pdict.sorted_paired_positions,var_dict[it2.first].pos );
                          int j = index_map[var_dict[it2.first].pos]; //if(node_loop == "10557985_G_G" || node_loop == "10558795_C_T" || node_loop == "10558197_A_A") {
                          //cout << node_loop << " " << i << "  " << j << "  " << pdict.sorted_paired_positions.size() << endl;
                          if (j < pdict.sorted_paired_positions.size() && j >= 0 && i != j) {
                                  link_matrix.add(i,j,var_dict[node_loop].var,var_dict[it2.first].var,it2.second);
                                  num_matrix.add_to(i,j,it2.second);
                          }
                  }
          }
          i++;
        }
        //////////////// create a second matrix which only includes connections with at least a certain weight (number of links)
        for (int i = 0; i < pdict.num_paired; i++) {
          for(auto const &ent1 : num_matrix.mat[i]) {
            if (ent1.second >= minimum_link_number) { num_matrix_second.set_val(i,ent1.first,ent1.second); }
          }
        }
        cout << "size of num matrix: " << diff_matrix.mat.size() << endl;
        t = clock() - t;
        cout << "created link matrix -- time: " << t << endl;
        link_matrix_calculations(link_matrix,num_matrix_second,span_matrix,corr_matrix,diff_matrix);
        pdict.get_submatrix_bounds(num_matrix_second);
        cout << "double positions size: " << pdict.double_positions.size() << endl;
        cout << "paired positions size: " << pdict.sorted_paired_positions.size() << endl;
        cout << "size of num matrix: " << diff_matrix.mat.size() << endl;
        pdict.hap_random_initialization();
}
