#ifndef anamuon_hh
#define anamuon_hh

//My
#include "ana.hh"
#include "dbscan.hh"

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
  
  anamuon(TString fileList) : ana(fileList), _dbs(new dbscan())
  {
  }
  
  anamuon(TString file, Int_t key) : ana(file, key), _dbs(new dbscan())
  {
  }

  void Loop(TString histOut);
  Double_t get_theta_p_t();

  void gen_ring(TGraph *gr, Int_t np, Double_t x0, Double_t y0, Double_t R);
  void get_average_approximate_ring_parameters(TGraph *gr, TRandom3 *rnd, Int_t niterations, Double_t &x0, Double_t &y0, Double_t &R);
  void get_approximate_ring_parameters(TGraph *gr, TRandom3 *rnd, Double_t &x0, Double_t &y0, Double_t &R);
  void get_Rvs_theta_and_theta_dist(TGraph *gr, Double_t x0, Double_t y0, Double_t R, TGraph *gr_R, TH1D *h1_theta_deg);
  void fit_muon_ring(TGraph *gr, TGraph *gr_frame, TGraph *gr_clean, TCanvas *c1, Int_t jentry_int,
		     Double_t x0app_average, Double_t y0app_average, Double_t Rapp_average,
		     Double_t x0out_fit, Double_t y0out_fit, Double_t Rout_fit,
		     Int_t event_id_val = -999,
		     Double_t energy_val = -999.0,
		     Double_t xcore_val = -999.0,
		     Double_t ycore_val = -999.0,
		     Double_t ev_time_val = -999.0,
		     Int_t nphotons_val = -999,
		     Int_t n_pe_val = -999,
		     Int_t n_pixels_val = -999,
		     Double_t azimuth_val = -999.0,
		     Double_t altitude_val = -999.0,
		     Double_t h_first_int_val = -999.0,
		     Double_t hmax_val = -999.0,
		     Double_t thetaDeg_val = -999.0,
		     TH1D *h1_pne_per_pix = NULL,
		     TString particle_type = "NONE",
		     TString pdf_out_file = "NONE");
  //
  static Double_t get_canonical_phi_from_V2phi(Double_t phiv2);
  void gen_azimuth_line(TGraph *gr, Int_t np, Double_t x0, Double_t y0, Double_t azimuth_deg);
  void gen_altitude_line(TGraph *gr, Int_t np, Double_t x0, Double_t y0, Double_t azimuth_val, Double_t altitude_val);

  void dbscan_cleaning( TGraph *gr, TGraph *gr_db, unsigned int minPts, Double_t eps);
  
  dbscan *_dbs;

};

#endif
