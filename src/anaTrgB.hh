#ifndef anaTrgB_hh
#define anaTrgB_hh

//My
#include "anabase.hh"
#include "anaConf.hh"

//root
#include <TROOT.h>
#include <TVector3.h>

class TChain;
class TFile;
class TTree;
class TString;
class TBranch;

class anaTrgB: public anabase {
public:

  anaTrgB(TString fileList) : anabase(fileList), _trg_conf_file("trg_setup.conf"), _NGB_rate_in_MHz(268.0), _fadc_electronic_noise_RMS(3.8082498), _name_ana_conf_file("./anaTrgB_p.conf"), _n_tot_ev_after_cuts(0.0), _n_tot_ev_after_cuts_processed(0.0)
  {
  }

  anaTrgB(TString file, Int_t key) : anabase(file, key), _trg_conf_file("trg_setup.conf"), _NGB_rate_in_MHz(268.0), _fadc_electronic_noise_RMS(3.8082498), _name_ana_conf_file("./anaTrgB_p.conf"), _n_tot_ev_after_cuts(0.0), _n_tot_ev_after_cuts_processed(0.0)
  {
  }

  void Loop(TString histOut,
	    Float_t NGB_rate_in_MHz = 268.0, Float_t fadc_electronic_noise_RMS = 3.8082498,
	    Int_t npe_min = 0, Int_t npe_max = 100000,
	    Int_t nEvSim_max = -999, Int_t nEv_max = -999,
	    Int_t rndseed = 12345, bool NGBsim = false);  

  Double_t get_n_tot_ev_after_cuts();

  inline void set_trg_conf_file(TString trgf){_trg_conf_file = trgf;};
  inline void set_ana_conf_file(TString val){_name_ana_conf_file = val;};
  
  static void digital_sum_rate_calculations( TH1D *h1_digital_sum_rate, TH1D *h1_digital_sum, TH1D *h1_npe, Int_t nn_fadc_point, Float_t fadc_sample_in_ns);
  static void dbc_nsb_rate_calculations( TH1D *h1_dbc_nsb_rate, TH1D *h1_dbc_number_of_points, TH1D *h1_npe, Int_t nn_fadc_point, Float_t fadc_sample_in_ns);
  static void rate_calculations( TH1D *h1_rate_vs_th, TH1D *h1_w, Double_t tot_flux, Double_t tot_number_of_sim_events,
				 Double_t n_tot_ev_after_cuts = 1.0, Double_t n_tot_ev_after_cuts_processed = 1.0);
  static Double_t get_proton_flux_pers_persr_perm2(Double_t e_min_GeV, Double_t e_max_GeV);
  static Double_t get_solid_angle(Double_t theta_deg);
  static Double_t calculate_simulation_sampling_factor(Double_t e_min_GeV, Double_t e_max_GeV);
  
private:
  //
  Bool_t cut();
  //
  Int_t _npe_min;
  Int_t _npe_max;
  Int_t _nEv_max;
  Int_t _nEvSim_max;
  //
  Double_t _n_tot_ev_after_cuts;
  Double_t _n_tot_ev_after_cuts_processed;
  //
  Float_t _NGB_rate_in_MHz;
  Float_t _fadc_electronic_noise_RMS;  
  //
  const void print_cuts();
  
  TString _trg_conf_file;
  TString _name_ana_conf_file;
  
  anaConf _anaConf;
};

#endif
