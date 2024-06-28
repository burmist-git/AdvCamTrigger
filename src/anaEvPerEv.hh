#ifndef anaEvPerEv_hh
#define anaEvPerEv_hh

//My
#include "anabase.hh"

//root
#include <TROOT.h>
#include <TVector3.h>

class TChain;
class TFile;
class TTree;
class TString;
class TBranch;

class anaEvPerEv: public anabase {
public:

  anaEvPerEv(TString fileList) : anabase(fileList)
  {
  }

  anaEvPerEv(TString file, Int_t key) : anabase(file, key)
  {
  }

  void save_event_to_bin_file(TString binFileOut, Int_t event_ID, Int_t rndseed);
  void get_events_map(TString txtFileOut);
  
private:
  //
};

#endif
