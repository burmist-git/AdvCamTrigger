#ifndef anamuon_hh
#define anamuon_hh

//My
#include "ana.hh"

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
class TRandom3;

class anamuon: public ana {
public:
  
  anamuon(TString fileList) : ana(fileList)
  {
  }
  
  anamuon(TString file, Int_t key) : ana(file, key)
  {
  }

  void Loop(TString histOut);
  Double_t get_theta_p_t();

  void gen_ring(TGraph *gr, Int_t np, Double_t x0, Double_t y0, Double_t R);
  void get_average_approximate_ring_parameters(TGraph *gr, TRandom3 *rnd, Int_t niterations, Double_t &x0, Double_t &y0, Double_t &R);
  void get_approximate_ring_parameters(TGraph *gr, TRandom3 *rnd, Double_t &x0, Double_t &y0, Double_t &R);
  void get_Rvs_theta_and_theta_dist(TGraph *gr, Double_t x0, Double_t y0, Double_t R, TGraph *gr_R, TH1D *h1_theta_deg);
  void fit_muon_ring(TGraph *gr, TGraph *gr_frame, TCanvas *c1, Int_t jentry_int, Double_t x0app_average, Double_t y0app_average, Double_t Rapp_average);
  
};

#endif
