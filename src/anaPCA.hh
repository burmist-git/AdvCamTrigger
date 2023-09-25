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
  
};

#endif
