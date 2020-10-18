#ifndef SV_JUNCTION_H
#define SV_JUNCTION_H

//////////////// c++ include //////////////////
#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>
#include <set>
#include <algorithm>
using namespace std;

//////////////// structures ///////////////////
class sv_entry {
    int num_cnx = 0;
public:
    std::string chr_one,chr_two;
    int pos_one,pos_two,tot_coverage;
    int intersection_length,union_length;
    //int bx_intersection_length,bx_union_length;
    std::pair<int,int> strand;
    std::vector<std::string> bx_intersection;
    std::vector<std::string> readid_tags1,readid_tags2;
    std::vector<std::string> read_intersection,read_union;
    //std::vector<std::string> sv1_bxtags,sv2_bxtags;
    //std::vector<std::string> tag_intersection,tag_union;
    std::unordered_map<std::string,std::string> bxmap1,bxmap2;
    std::unordered_map<int,int> hap_end1,hap_end2;
    std::unordered_map<int,int> num_end1,num_end2;
    void set_locations(std::string,int,std::string,int,int,int,int);
    void add_tag_read1(std::string,std::string);
    void add_tag_read2(std::string,std::string);
    void read_id_intersection();
    //void sv_bx_intersection();
    void het_map(int,int,bool);
};

#endif // SV_JUNCTION_H
