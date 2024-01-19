#ifndef anaTrgA_hh
#define anaTrgA_hh

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

class anaTrgA: public anabase {
public:

  anaTrgA(TString fileList) : anabase(fileList)
  {
  }

  anaTrgA(TString file, Int_t key) : anabase(file, key)
  {
  }

  void Loop(TString histOut, Int_t binE, Int_t binTheta, Int_t binDist, Int_t npe_min, Int_t npe_max, Int_t nEv_max);

  void SiPM_dist(TString histOut);

  void transform_SiPM_distadd( TGraph *gr_tr, const TGraph *gr);

  void test_single_pe_amplitude_generator(TString histOut);
  
private:
  //
  Bool_t cut(Int_t nevsim, Double_t theta_p_t_deg, Double_t rcore);
  //
  Double_t _E_min;
  Double_t _E_max;
  Double_t _theta_deg_min;
  Double_t _theta_deg_max;
  Double_t _dist_min;
  Double_t _dist_max;
  Int_t _npe_min;
  Int_t _npe_max;
  Int_t _nEv_max;
  //
  const void print_cuts();
  //
  Double_t get_theta_p_t( const TVector3 &v_det, Double_t altitude_v, Double_t azimuth_v);
  
};

#endif
