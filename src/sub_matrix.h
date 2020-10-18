#ifndef SUB_MATRIX_H
#define SUB_MATRIX_H

//////////////// c++ include //////////////////
#include <vector>
using namespace std;

//////////////// linker include ////////////////
#include "map_matrix.h"

///////////////////////////////////////////////
class submatrix_opt {
public:
    int length;
    std::vector< std::vector<double> > smat;
    void initialize(vector<int>&,fmatrix&);
};

#endif  // SUB_MATRIX_H
