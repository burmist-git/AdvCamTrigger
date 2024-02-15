#pragma once

//root
#include <TString.h>

//c, c++
#include <string>
#include <vector>

class evstHist;
class TH1D;

class rateCalculator {
public:
  
  rateCalculator();
  rateCalculator( const char* name, const char* title, TString hist_file_prefix, TString output_hist_file, evstHist* evH_flux,
		  Int_t n_jobs = -999, bool disable_energy_theta_rcore_binwise_cuts = false);
  ~rateCalculator();
  
private:
  
  TString _name;
  TString _title;
  TString _hist_file_prefix;
  TString _output_hist_file;
  evstHist *_evH_flux;
  
  bool get_histogram_file_name(TString hist_file_prefix, Int_t i_binE, Int_t i_binTheta, Int_t i_binDist, TString &hist_file_name, Int_t jobID = -999);
  bool get_histogram_file_name(TString  hist_file_dir, TString file_prefix, TString &hist_file_name, Int_t jobID);
  void calculate_rate_binwise();
  void calculate_rate();
  void calculate_rate_tot();
  void calculate_rate_NSB(Int_t n_jobs_NSB = 100, TString hist_file_dir = "../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/nsb_1x_268MHz/trgA/");
  void append_histogram(TH1D *h1, TH1D *h1_to_add);
  void get_histogram(TString hist_file_name, TH1D *h1, TString hist_name_to_get = "h1_dbc_number_of_points_w");
  
  void copyHist(TH1D *h1, TH1D *h1_to_cp);
  void addHist(TH1D *h1, TH1D *h1_to_add);
  
  TH1D *_h1_rate;
  TH1D *_h1_N_dbc;
  TH1D *_h1_N_dbc_mean;
  TH1D *_h1_ev_per_job;
  TH1D *_h1_dbc_number_of_points_tot_w;
  TH1D *_h1_rate_vs_th;
  TH1D *_h1_rate_NSB_vs_th;
  TH1D *_h1_rate_tot_vs_th;
  TH1D *_h1_dbc_number_of_points_NSB;
  TH1D *_h1_N_dbc_mean_NSB;
  TH1D *_h1_dbc_mean_time_ii;
  TH1D *_h1_dbc_mean_time_ii_NSB;
  Double_t _n_ev_tot;
  Double_t _NSB_rate_lowth; 
  std::vector<TH1D*> _h1_dbc_number_of_points_w_v;

  bool _disable_energy_theta_rcore_binwise_cuts;
  Int_t _n_jobs;
  
  void save_output();
  
};
