#ifndef anashort_hh
#define anashort_hh

//My
#include "anabase.hh"

//root
#include <TROOT.h>

class TChain;
class TFile;
class TTree;
class TString;
class TBranch;

class anashort: public anabase {
public:

  anashort(TString fileList) : anabase(fileList, true)
  {
  }

  anashort(TString file, Int_t key) : anabase(file, key, true)
  {
  }

  void Loop(TString histOut);
  
};

#endif
