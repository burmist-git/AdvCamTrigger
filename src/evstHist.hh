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
  TCanvas* Draw_hist_core(Int_t e_bin_i, TString fileName, TString frame_title);
  
  inline TH1D* get_theta_hist() {return _h1_theta;}
  inline TH1D* get_E_hist() {return _h1_E;}
  inline std::vector<TH1D*> get_v_r() {return _v_r;};
  
  TString _hist_name;
  TString _hist_title;

  TString _title;
  
  void Divide(evstHist *evH_cut, evstHist *evH_all, bool with_r_core = false);
  void Multiply(evstHist *evH_eff, evstHist *evH_flux, bool with_r_core = false);
  void DumpBinContent(TString data_out, bool with_r_core = false);
  void LoadBinContent(TString data_in, bool with_r_core = false);

  Double_t GetTotIntegral();
  Double_t GetIntegral(Double_t e_min, Double_t e_max, Double_t theta_min, Double_t theta_max) const;
  
  static const void PrintBinsInfo(const TH1D *h1);

  static void set_r_core_bins(TH1D *h1, Double_t r_core_max = 2000);
  static void init_core_hist(TH1D *h1);
  
  Int_t get_bin_ID( Double_t E, Double_t th);

  void Fill_rcore( Double_t th, Double_t E, Double_t r_core);

  void Get_th_bin_ID_and_e_bin_ID( Int_t cellID, Int_t &th_bin_ID, Int_t &e_bin_ID);

  inline Double_t get_Emin() const {return _Emin;};
  inline Double_t get_Emax() const {return _Emax;};
  inline Int_t get_N_bins_E() const {return _N_bins_E;};
  inline Double_t get_Thetamin() const {return _Thetamin;};
  inline Double_t get_Thetamax() const {return _Thetamax;};
  inline Int_t get_N_bins_t() const {return _N_bins_t;};

  bool check_bin_compatibility(const evstHist *evH, bool with_r_core = false);
  bool check_r_core_bin_compatibility(const evstHist *evH);
  
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

  std::vector<Int_t> _v_cellID_th_bin_ID;
  std::vector<Int_t> _v_cellID_e_bin_ID;
  
};
