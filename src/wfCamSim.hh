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
  //void gen_WF( TGraph *gr_wf, TGraph *gr_wf_sig, TGraph *gr_wf_sig_only, unsigned int n_signals, TH1D *h1_photon_time);
  //void gen_WF( TGraph *gr_wf, TGraph *gr_wf_sig, TGraph *gr_wf_sig_only, unsigned int n_signals);
  
private:
  //  
  TRandom3 *_rnd;
  TString _wf_tamplete;
  TString _spe_dat;
  //
  TGraph *_gr_wf_ampl;
  TGraph *_gr_wf_tmpl;
  TH1D *_h1_wf_ampl;
  Double_t _Ampl_Prompt_max;
  Double_t _Prompt_max;
  Double_t _t_max_ampl_wf_tmpl;
  //
  //
  //Double_t DCR_rate;
  //Int_t npe_tot;
  //Double_t *pe_time;
  //Int_t *pe_ch;
  //
};

#endif
