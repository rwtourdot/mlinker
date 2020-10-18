#ifndef READ_TREE_H
#define READ_TREE_H

//////////////// c++ include //////////////////
#include <vector>
#include <string>
#include <unordered_map>
#include <algorithm>
using namespace std;

//////////////// structures ///////////////////
class read_tree {
public:
    std::string readname;
    int num_cnx = 0;
    std::vector<std::string> connected_strings;
    std::vector<std::string> connected_readnames;
    std::vector<int> connected_ints;
    std::vector<std::string> connected_chr;
    void set_name(std::string);   //void add_connection(std::string);
    void add_connection(std::string,std::string);
    void add_connection_int(int,std::string);
};

//////////////// definitions //////////////////
typedef std::unordered_map<std::string,read_tree> read_graph;

#endif  // READ_TREE_H
