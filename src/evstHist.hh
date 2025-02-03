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
  void test_Get_th_bin_ID_and_e_bin_ID();
  void test_get_bin(Double_t E, Double_t th, Double_t val);
  TCanvas* Draw_hist(TString fileName, TString frame_title="");
  TCanvas* Draw_hist_core(Int_t e_bin_i, TString fileName, TString frame_title);
  
  inline TH1D* get_theta_hist() {return _h1_theta;}
  inline TH1D* get_E_hist() {return _h1_E;}
  inline std::vector<TH1D*>& get_v_r() {return _v_r;};
  inline std::vector<TH1D*>& get_h1_arb_v() {return _h1_arb_v;};
  
  TString _hist_name;
  TString _hist_title;

  TString _title;
  
  void selfNorm();
  void Divide(evstHist *evH_cut, evstHist *evH_all, bool with_r_core = false);
  void Divideh1(evstHist *evH_cut, evstHist *evH_all, Double_t normval = 1.0);
  void Multiply(evstHist *evH_eff, evstHist *evH_flux, bool with_r_core = false);
  void Multiply_core(evstHist *evH);
  void DumpBinContent(TString data_out, bool with_r_core = false);
  void LoadBinContent(TString data_in, bool with_r_core = false);

  Double_t GetTotIntegral();
  Double_t GetIntegral(Double_t e_min, Double_t e_max, Double_t theta_min, Double_t theta_max) const;
  
  static const void PrintBinsInfo(const TH1D *h1, Int_t binID=-1);
  const void Print_hist_BinsInfo(Int_t binE = -1, Int_t binTheta = -1, Int_t binDist =-1);
  const void get_Bin_Edge(const TH1D *h1, Int_t binID, Double_t &bin_l, Double_t &bin_r);
  
  static void set_r_core_bins(TH1D *h1, Double_t r_core_max = 2000);
  static void init_core_hist(TH1D *h1);

  Double_t Get_r_core_max_val();
  
  Int_t get_bin_ID( Double_t E, Double_t th);

  void Fill_rcore( Double_t th, Double_t E, Double_t r_core, Double_t weight_val = 1.0);

  bool Get_th_bin_ID_and_e_bin_ID( Int_t cellID, Int_t &th_bin_ID, Int_t &e_bin_ID);

  inline Double_t get_Emin() const {return _Emin;};
  inline Double_t get_Emax() const {return _Emax;};
  inline Int_t get_N_bins_E() const {return _N_bins_E;};
  inline Double_t get_Thetamin() const {return _Thetamin;};
  inline Double_t get_Thetamax() const {return _Thetamax;};
  inline Int_t get_N_bins_t() const {return _N_bins_t;};

  bool check_bin_compatibility(const evstHist *evH, bool with_r_core = false);
  bool check_r_core_bin_compatibility(const evstHist *evH);

  void Weight_TeV();
  Double_t Weight_integral_GeV( Double_t e_r_GeV, Double_t e_l_GeV);

  static Double_t get_Weight_ETeV(Double_t ETev);
  static Double_t get_Weight_ETeV_soft(Double_t ETev);
  
  bool get_arbitrary_hist_ID( Double_t E, Double_t th, Double_t r_core, unsigned int &arbitrary_hist_ID);

  void test_get_arbitrary_hist_ID();

  void init_h1_arb_v(const char* name, const char* title, Int_t nBins, Double_t valmin, Double_t valmax);
  void fill_h1_arb_v(Double_t E, Double_t th, Double_t r_core, Double_t val, Double_t valweight = -999.0);
  
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

  Double_t _r_core_max_val;
  Int_t _N_bins_r_core;
  
  Double_t _d;
  Double_t _le;
  Double_t _lt;

  TH1D* _h1_theta;
  TH1D* _h1_E;
  
  std::vector<TH1D*> _v_r;

  std::vector<Int_t> _v_cellID_th_bin_ID;
  std::vector<Int_t> _v_cellID_e_bin_ID;

  std::vector<TH1D*> _h1_arb_v;
  
};
