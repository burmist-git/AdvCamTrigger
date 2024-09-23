#ifndef anastereo_hh
#define anastereo_hh

//My
#include "anabasestereo.hh"
#include "anaConf.hh"

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
  
  anastereo(TString fileList);
  anastereo(TString file, Int_t key);

  void Loop(TString histOut);
  void get_event_acceptance_vs_dt( TH1D *h1_dtime, TH1D *h1_event);
  inline void set_ana_conf_file(TString val){_name_ana_conf_file = val;};

private:

  TString _name_ana_conf_file;
  anaConf _anaConf;
  
};

#endif
