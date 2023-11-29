#pragma once

//root
#include <TObject.h>
#include <TH2Poly.h>
#include <TGraph.h>
#include <TVector2.h>
#include <TCanvas.h>
#include <TH1D.h>

//c, c++
#include <string>
#include <vector>
#include <map>

class evstHist: public TH2Poly {
 public:
  
  evstHist();
  evstHist(const char* name, const char* title,
	   Double_t val_Emin = 1.0, Double_t val_Emax = 100000, Int_t val_N_bins_E = 25,
	   Double_t val_Thetamin = 0.0, Double_t val_Thetamax = 10.0, Int_t val_N_bins_t = 10);
  ~evstHist();
  void test();
  void test_get_bin(Double_t E, Double_t th, Double_t val);
  TCanvas* Draw_hist(TString fileName, TString frame_title="");

  inline TH1D* get_theta_hist() {return _h1_theta;}
  inline TH1D* get_E_hist() {return _h1_E;}
  inline std::vector<TH1D*> get_v_r() {return _v_r;};
  
  TString _hist_name;
  TString _hist_title;

  TString _title;
  
  void Divide(evstHist *evH_cut, evstHist *evH_all);
  void Multiply(evstHist *evH_eff, evstHist *evH_flux);
  void DumpBinContent(TString data_out);
  void LoadBinContent(TString data_in);

  Double_t GetTotIntegral();
  Double_t GetIntegral(Double_t e_min, Double_t e_max, Double_t theta_min, Double_t theta_max) const;
  
  static const void PrintBinsInfo(const TH1D *h1);

  static void set_r_core_bins(TH1D *h1, Double_t r_core_max = 2000);

  Int_t get_bin_ID( Double_t E, Double_t th);

  void Fill_rcore( Double_t th, Double_t E, Double_t r_core);
  
 private:

  Double_t _Emin;
  Double_t _Emax;
  Int_t _N_bins_E;

  Double_t *_Eminarr;
  Double_t *_Emaxarr;

  Double_t *_Thetaminarr;
  Double_t *_Thetamaxarr;
  
  Double_t _Thetamin;
  Double_t _Thetamax;
  Int_t _N_bins_t;
  
  Double_t _d;
  Double_t _le;
  Double_t _lt;

  TH1D* _h1_theta;
  TH1D* _h1_E;
  
  std::vector<TH1D*> _v_r;
  
};
