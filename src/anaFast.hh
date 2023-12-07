#ifndef anaFast_hh
#define anaFast_hh

//My
#include "anabase.hh"
#include "anashort.hh"
#include "anaConf.hh"

//root
#include <TROOT.h>
#include <TVector3.h>

class TString;
class TVector3;
struct anaConf;

class anaFast: public anashort {
public:

  anaFast(TString fileList, TString anaConfFile);

  anaFast(TString fileList) : anashort(fileList)
  {
  }

  anaFast(TString file, Int_t key) : anashort(file, key)
  {
  }

  void Loop(TString histOut);
  void Loop_scan(TString histOut, TString simtel_all_dat);
  bool cuts(Double_t theta_p_t_deg = 10000);
  
  Double_t _r_core;
  Double_t _theta_core;

  anaConf _anaConf;
  TVector3 _v_det;

  Double_t get_theta_p_t();
  
};

#endif
