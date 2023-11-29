#ifndef anaTrg_hh
#define anaTrg_hh

//My
#include "anabase.hh"

//root
#include <TROOT.h>

class TChain;
class TFile;
class TTree;
class TString;
class TBranch;

class anaTrg: public anabase {
public:

  anaTrg(TString fileList) : anabase(fileList), _npe_max(500)
  {
  }

  anaTrg(TString file, Int_t key) : anabase(file, key), _npe_max(500)
  {
  }

  void Loop(TString histOut);
  static void set_n_pe_bins(TH1D *h1, Int_t &npe_max, Int_t &npe_min);
  static void dump_n_pe_bins(TH1D *h1);

private:
  
  Bool_t cut(TH1D *h1 = NULL);
  Int_t _npe_max;
  Int_t _npe_min;
  
};

#endif
