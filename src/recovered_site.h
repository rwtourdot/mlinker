#ifndef RECOVERED_SITE_H
#define RECOVERED_SITE_H

//////////////// c++ include //////////////////
#include <string>
#include <vector>
#include <unordered_map>
#include <set>
#include <iostream>
using namespace std;

class recovered_node {
public:
  int pos;
  int hap;
  int ref_numA,ref_numB,alt_numA,alt_numB;
  int upos_ref_numA,upos_ref_numB,upos_alt_numA,upos_alt_numB;
  int bx_ref_numA,bx_ref_numB,bx_alt_numA,bx_alt_numB;
  int ref_numA_hashcall,ref_numB_hashcall,alt_numA_hashcall,alt_numB_hashcall;
};

#endif  // RECOVERED_SITE_H
