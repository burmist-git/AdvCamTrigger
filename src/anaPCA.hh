#ifndef anaPCA_hh
#define anaPCA_hh

//My
#include "anabase.hh"
#include "anashort.hh"
#include "anaConf.hh"

//root
#include <TROOT.h>

class TString;
struct anaConf;

class anaPCA: public anashort {
public:

  anaPCA(TString fileList, TString anaConfFile);

  anaPCA(TString fileList) : anashort(fileList)
  {
  }

  anaPCA(TString file, Int_t key) : anashort(file, key)
  {
  }

  void Loop(TString histOut);
  bool cuts();

  Double_t _r_core;
  Double_t _theta_core;

  anaConf _anaConf;

  void load_S_Vh_data(TString name_S, TString name_Vh);
  static const Int_t _dd_im = 887;
  Double_t _data_S[_dd_im];
  Double_t _data_Vh[_dd_im][_dd_im];
  
};

#endif
