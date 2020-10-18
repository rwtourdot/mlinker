#ifndef MAP_MATRIX_H
#define MAP_MATRIX_H

//////////////// c++ include //////////////////
#include <vector>
#include <string>
#include <map>
using namespace std;

///////////////////////////////////////////////
template <class T>
class map_matrix
{
public:
    int length;
    typedef std::map<size_t, std::map<size_t , T> > mat_t;
    typedef typename mat_t::iterator row_iter;
    typedef std::map<size_t, T> col_t;
    typedef typename col_t::iterator col_iter;
    mat_t mat;
    map_matrix(){ }
    map_matrix(size_t i){ m=i; n=i; length = i; }
    map_matrix(size_t i, size_t j){ m=i; n=j; length = i; }
    T& operator()(size_t i, size_t j) { if(i>=m || j>=n) throw; return mat[i][j]; }
    T operator()(size_t i, size_t j) const { if(i>=m || j>=n) throw; return mat[i][j]; }
    void add_to(int i, int j, T number) { mat[i][j] = mat[i][j] + number; }
    void set_val(int i, int j, T value) { mat[i][j] = value; }
private:
    size_t m;
    size_t n;
};

///////////////////////////////////////////////
class map_matrix_vector
{
public:
    int length;
    typedef std::map<size_t, std::map<size_t , std::vector<bool> > > mat_t;
    typedef typename mat_t::iterator row_iter;
    typedef std::map<size_t, std::vector<bool> > col_t;
    typedef typename col_t::iterator col_iter;
    mat_t mat;
    map_matrix_vector(size_t i){ m=i; n=i; length = i; }
    map_matrix_vector(size_t i, size_t j){ m=i; n=j; length = i; }
    bool return_val(int i, int j, int m) { return mat[i][j][m]; }
    void add(int i, int j, bool b1 ,bool b2 , int number) {
        bool newbool = true;
        if (b1 && !b2) { newbool = false; }
        if (!b1 && b2) { newbool = false; }
        for  (int m = 0; m < number; m++) { mat[i][j].push_back(newbool); };
    }
private:
    size_t m;
    size_t n;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
typedef map_matrix<double> fmatrix;
typedef map_matrix<int> imatrix;

#endif  // MAP_MATRIX_H
