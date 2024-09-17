#ifndef anastereo_hh
#define anastereo_hh

//My
#include "anabasestereo.hh"

//root
#include <TROOT.h>

//C, C++
#include <vector>

class TChain;
class TFile;
class TTree;
class TString;
class TBranch;
class TGraph;

class anastereo: public anabasestereo {
public:
  
  anastereo(TString fileList) : anabasestereo(fileList)
  {
  }

  anastereo(TString file, Int_t key) : anabasestereo(file, key)
  {
  }

  void Loop(TString histOut);
  
};

#endif
