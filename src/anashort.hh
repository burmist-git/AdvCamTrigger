#ifndef ana_hh
#define ana_hh

//My
#include "anabaseshort.hh"

//root
#include <TROOT.h>

class TChain;
class TFile;
class TTree;
class TString;
class TBranch;

class anashort: public anabaseshort {
public:

  anashort(TString fileList) : anabaseshort(fileList)
  {
  }

  anashort(TString file, Int_t key) : anabaseshort(file, key)
  {
  }

  void Loop(TString histOut);
  
};

#endif
