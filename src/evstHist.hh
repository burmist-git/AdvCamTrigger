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
  TCanvas* Draw_hist(TString fileName);

  inline TH1D* get_theta_hist() {return _h1_theta;}
  inline TH1D* get_E_hist() {return _h1_E;}

  TString _hist_name;
  TString _hist_title;

  TString _title;
  
  void Divide(evstHist *evH_cut, evstHist *evH_all);
  void Multiply(evstHist *evH_eff, evstHist *evH_flux);
  void DumpBinContent(TString data_out);
  void LoadBinContent(TString data_in);

  static const void PrintBinsInfo(const TH1D *h1);
  
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
  
};
