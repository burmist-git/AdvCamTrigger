#ifndef anaSuperFast_hh
#define anaSuperFast_hh

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

class anaSuperFast: public anashort {
public:

  anaSuperFast(TString fileList, TString anaConfFile);

  anaSuperFast(TString fileList) : anashort(fileList)
  {
  }

  anaSuperFast(TString file, Int_t key) : anashort(file, key)
  {
  }

  void Loop(TString histOut);  

  anaConf _anaConf;
  TVector3 _v_det;

  Double_t get_theta_p_t();
  void get_cumulative( TH1D *h1, TH1D *h1_int);
  void get_cumulative( TH1D *h1, TH1D *h1_int, Double_t norm);
  void save_hist_to_csv(TString outfilename, TH1D *h1);
  
};

#endif
