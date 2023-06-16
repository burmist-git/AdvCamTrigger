#ifndef wfCamSim_hh
#define wfCamSim_hh

//C, C++
#include <iostream>
#include <vector>

//root
#include "TROOT.h"

//root
class TString;
class TRandom3;
class TGraph;
class TH1D;

class wfCamSim {

public :
  wfCamSim( TRandom3 *rnd, TString wf_tamplete, TString spe_dat);
  ~wfCamSim();

public :
  //
  void getWF_ampl(TString name, Double_t &Ampl_Prompt_max, Double_t &Prompt_max);
  void getWF_tmpl(TString name);
  Double_t generate_wf_ampl_from_file();
  const inline TGraph *getTemplate() {return _gr_wf_tmpl;}
  const inline TGraph *get_gr_wf_ampl() {return _gr_wf_ampl;}
  const inline TH1D *get_h1_wf_ampl_ADC() {return _h1_wf_ampl_ADC;}
  const inline TH1D *get_h1_wf_ampl() {return _h1_wf_ampl;}
  void setupe_ADC_ampl(Double_t bits_per_pe = 8.25, Int_t n_max_pe = 10);
  void test_single_pe_amplitude_generator( TH1D *h1, Int_t n_pe_to_sim);
  void test_single_pe_amplitude_from_hist_generator( TH1D *h1, Int_t n_pe_to_sim);
  //
  //
  void simulate_cam_event(const Int_t nn_fadc_point,
			  const Int_t nn_PMT_channels,
			  std::vector<std::vector<Int_t>> &wf,
			  const Float_t NGB_rate_in_MHz,
			  const Float_t ev_time,
			  const Float_t time_offset,
			  const Int_t n_pe,
			  const Int_t *pe_chID,
			  const Float_t *pe_time);
  void simulate_cam_event(const Int_t nn_fadc_point,
			  const Int_t nn_PMT_channels,
			  std::vector<std::vector<Int_t>> &wf,
			  const Float_t NGB_rate_in_MHz);
  //
  
private:
  //
  static Double_t integrate_spe( TGraph *gr, Double_t x0, Double_t Dx);
  //  
  TRandom3 *_rnd;
  TString _wf_tamplete;
  TString _spe_dat;
  //
  TGraph *_gr_wf_ampl;
  TH1D *_h1_wf_ampl;
  TH1D *_h1_wf_ampl_ADC;
  static const unsigned int _n_ADC_arr = 2000000000;
  unsigned char *_ampl_ADC_arr;
  TGraph *_gr_wf_tmpl;
  Double_t _n_ADC_max_for_generator;
  Int_t generate_single_pe_amplitude();
  Int_t generate_single_pe_amplitude_from_hist();
  //  
  Double_t _Ampl_Prompt_max;
  Double_t _Prompt_max;
  Double_t _t_max_ampl_wf_tmpl;
};

#endif
