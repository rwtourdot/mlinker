#include "sub_matrix.h"

void submatrix_opt::initialize( vector<int> &sparse_map, fmatrix &diff_matrix ) {
    length = sparse_map.size();
    smat.resize(length);
    for (int i = 0; i < length; i++) { smat[i].resize(length); }
    for (int i = 0; i < length; i++) {
        int idex = sparse_map[i];
        for (int j = 0; j < length; j++) {
                int jdex = sparse_map[j];
                smat[i][j] = diff_matrix(idex,jdex);
        }
    }
};

