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

  anaTrgA(TString fileList) : anabase(fileList), _disable_energy_theta_rcore_binwise_cuts(false), _k_dist_graph_flag(false)
  {
  }

  anaTrgA(TString file, Int_t key) : anabase(file, key), _disable_energy_theta_rcore_binwise_cuts(false), _k_dist_graph_flag(false)
  {
  }

  void Loop(TString histOut, Int_t binE, Int_t binTheta, Int_t binDist, Int_t npe_min, Int_t npe_max, Int_t nEv_max, Int_t rndseed, Int_t data_chunk_ID = -999, bool NGBsim = false);

  void SiPM_dist(TString histOut);

  void transform_SiPM_distadd( TGraph *gr_tr, const TGraph *gr);

  void test_single_pe_amplitude_generator(TString histOut);

  inline void set_disable_energy_theta_rcore_binwise_cuts(bool val){_disable_energy_theta_rcore_binwise_cuts = val;};
  inline void set_k_dist_graph_flag(bool val){_k_dist_graph_flag = val;};
  
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
  bool _disable_energy_theta_rcore_binwise_cuts; 
  //
  const void print_cuts();
  //
  Double_t get_theta_p_t( const TVector3 &v_det, Double_t altitude_v, Double_t azimuth_v);

  bool _k_dist_graph_flag;
  
};

#endif
