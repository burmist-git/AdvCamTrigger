#ifndef anaPCAp_hh
#define anaPCAp_hh

//My
#include "anabase.hh"
#include "anashort.hh"
#include "anaConf.hh"
#include <vector>

//root
#include <TROOT.h>

class TString;
struct anaConf;

class anaPCAp: public anashort {
public:

  anaPCAp(TString fileList, TString anaConfFile);

  anaPCAp(TString fileList) : anashort(fileList)
  {
  }

  anaPCAp(TString file, Int_t key) : anashort(file, key)
  {
  }

  void Loop(TString histOut);
  void draw_principal(TString histOut);
  void draw_reco( TString histOut, Int_t nEv_max, TString recofile);
  bool cuts();

  anaConf _anaConf;

  Double_t _r_core;
  Double_t _theta_core;
  
  void load_S_Vh_data(TString name_S, TString name_Vh);
  void load_shower_data(TString name_shower, Int_t nEv_max, std::vector<std::vector<Double_t>> &shower_v);
  static const Int_t _dd_im = 887;
  Double_t _data_S[_dd_im];
  Double_t _data_Vh[_dd_im][_dd_im];

  Double_t _phi0_shift;
  
};

#endif
