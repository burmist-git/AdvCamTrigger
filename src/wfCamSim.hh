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
class anabase;

class wfCamSim {

public :
  //
  wfCamSim(TRandom3 *rnd, TString wf_tamplete, TString spe_dat,
	   const unsigned int nn_fadc_point,
	   const unsigned int nn_PMT_channels,
	   const Float_t fadc_offset,
	   const Float_t fadc_sample_in_ns,
	   const Float_t NGB_rate_in_MHz);
  wfCamSim(TRandom3 *rnd, TString wf_tamplete, TString spe_dat,
	   const unsigned int nn_fadc_point,
	   const unsigned int nn_PMT_channels,
	   const Float_t fadc_offset,
	   const Float_t fadc_sample_in_ns,
	   const Float_t NGB_rate_in_MHz,
	   Float_t fadc_electronic_noise_RMS);
  ~wfCamSim();
  //
  void getWF_ampl(TString name, Double_t &Ampl_Prompt_max, Double_t &Prompt_max);
  void getWF_tmpl(TString name);
  Double_t generate_wf_ampl_from_file();
  const inline TGraph *getTemplate() {return _gr_wf_tmpl;}
  const inline TGraph *get_gr_wf_ampl() {return _gr_wf_ampl;}
  const inline TH1D *get_h1_wf_ampl_ADC() {return _h1_wf_ampl_ADC;}
  const inline TH1D *get_h1_wf_ampl() {return _h1_wf_ampl;}
  const inline TH1D *get_h1_adc_NGB_pedestal(){return _h1_adc_NGB_pedestal;}
  const inline TH1D *get_h1_dadc_NGB_pedestal(){return _h1_dadc_NGB_pedestal;}
  //
  void setupe_ADC_ampl(Double_t bits_per_pe = 8.25, Int_t n_max_pe = 10);
  void test_single_pe_amplitude_generator( TH1D *h1, Int_t n_pe_to_sim);
  void test_single_pe_amplitude_from_hist_generator( TH1D *h1, Int_t n_pe_to_sim);
  //
  //
  void simulate_cam_event(const Int_t nn_fadc_point,
			  const Int_t nn_PMT_channels,
			  std::vector<std::vector<Int_t>> &wf,
			  const Float_t ev_time,
			  const Float_t time_offset,
			  const Int_t n_pe,
			  const Int_t *pe_chID,
			  const Float_t *pe_time);
  void simulate_cam_event(const Int_t nn_fadc_point,
			  const Int_t nn_PMT_channels,
			  std::vector<std::vector<Int_t>> &wf);
  //
  void print_wfCamSim_configure();
  void get_gr_WF_tmpl_array(TGraph *gr);
  //
  static void calculate_camera_ADC_mean_and_std(const std::vector<std::vector<Int_t>> &wf, Float_t &meanv, Float_t &stdv);
  static void calculate_camera_ADC_mean_and_std(const std::vector<std::vector<Int_t>> &wf, Float_t &meanv, Float_t &stdv, TH1D *h1_adc, TH1D *h1_dadc);
  static Int_t get_charge(const std::vector<int>& wf, Int_t offset);
  Int_t get_charge(const std::vector<int>& wf);
  
  void test_calculate_pedestal(TH1D *h1_adc, TH1D *h1_dadc);
  void generate_gif_for_event(TString pathPref, Int_t event_id,
			      const std::vector<std::vector<Int_t>> &wf);
  void generate_gif_for_event(TString pathPref, Int_t event_id,
			      const std::vector<std::vector<Int_t>> &wf,
			      const std::vector<std::vector<Int_t>> &wf_ref,
			      const std::vector<std::vector<unsigned int>> &trg_vector,
			      const anabase *ab = NULL);
  
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
  //
  TH1D *_h1_adc_NGB_pedestal;
  TH1D *_h1_dadc_NGB_pedestal;
  //
  static const unsigned int _n_ADC_arr = 2000000000;
  Double_t _n_ADC_max_for_generator;
  unsigned char *_ampl_ADC_arr;
  TGraph *_gr_wf_tmpl;
  Float_t *_arr_wf_tmpl;
  unsigned int _n_arr_wf_tmpl;
  Int_t generate_single_pe_amplitude();
  Int_t generate_single_pe_amplitude_from_hist();
  //
  void generateNGB(std::vector<int> &wf);
  void generate_zero_wf(std::vector<int> &wf);
  void generate_zero_wf(std::vector<int> &wf, Int_t pedestal);
  void generate_wf(std::vector<int> &wf, Float_t pe_time);
  void generate_wf_from_gr(std::vector<int> &wf, Float_t pe_time);
  void generate_electronic_noise(std::vector<int> &wf);
  void generate_electronic_noise_pedestal_removal(std::vector<int> &wf);
  //
  void calculate_pedestal();
  void calculate_pedestal(TH1D *h1_adc, TH1D *h1_dadc);
  //
  void generateWF_tmpl_array();
  //
  bool check_pe_chID(Int_t id_ch);
  //  
  Double_t _Ampl_Prompt_max;
  Double_t _Prompt_max;
  Double_t _t_max_ampl_wf_tmpl;
  //
  Float_t _fadc_offset;
  Float_t _fadc_sample_in_ns;
  Float_t _NGB_rate_in_MHz;
  //
  Float_t _NGB_pedestal_mean;
  Float_t _NGB_pedestal_std;
  //  
  Float_t _t_left_in_ns;
  Float_t _t_righ_in_ns;
  Float_t _n_NGB_pe_average_in_window;
  //
  Float_t _wf_tmpl_t_left;
  Float_t _wf_tmpl_t_right;
  Float_t _wf_tmpl_t_start;
  Float_t _wf_tmpl_t_stop;
  //
  Float_t _fadc_electronic_noise_RMS;
  Int_t _fadc_pedestal;
  //
  unsigned int _nn_fadc_point;
  unsigned int _nn_PMT_channels;
  //
  Float_t _dt_arr_wf_tmpl;
};

#endif
