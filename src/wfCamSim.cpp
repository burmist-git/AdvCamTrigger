//my
#include "wfCamSim.hh"
#include "sipmCameraHist.hh"
#include "anabase.hh"

//root
#include "TRandom3.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TMath.h"
#include "TH1D.h"
#include "TColor.h"

//C, C++
#include <iostream>
#include <vector>
#include <fstream>
#include <assert.h>
#include <iomanip>
#include <stdlib.h>
#include <fstream>
#include <time.h>
#include <bits/stdc++.h>
#include <sys/stat.h>

wfCamSim::wfCamSim( TRandom3 *rnd, TString wf_tamplete, TString spe_dat,
		    const unsigned int nn_fadc_point,
		    const unsigned int nn_PMT_channels,
		    const Float_t fadc_offset,
		    const Float_t fadc_sample_in_ns,
		    const Float_t NGB_rate_in_MHz) : wfCamSim( rnd, wf_tamplete, spe_dat, nn_fadc_point,
							       nn_PMT_channels, fadc_offset,
							       fadc_sample_in_ns, NGB_rate_in_MHz, 3.94)
{ 
}

wfCamSim::wfCamSim( TRandom3 *rnd, TString wf_tamplete, TString spe_dat,
		    const unsigned int nn_fadc_point,
		    const unsigned int nn_PMT_channels,
		    const Float_t fadc_offset,
		    const Float_t fadc_sample_in_ns,
		    const Float_t NGB_rate_in_MHz,
		    Float_t fadc_electronic_noise_RMS) : _fadc_electronic_noise_RMS(fadc_electronic_noise_RMS), _electronic_noise_with_pedestal_removal(true)
{
  //
  _h1_wf_ampl_ADC = NULL;
  _ampl_ADC_arr = NULL;
  //
  _rnd = rnd;
  _wf_tamplete = wf_tamplete;
  _spe_dat = spe_dat;
  //
  //_fadc_electronic_noise_RMS = 3.94;
  //_fadc_electronic_noise_RMS = 0.01;
  _fadc_pedestal = 300;
  //_fadc_pedestal = 0;
  //
  _wf_tmpl_t_start =  0.0; //ns
  _wf_tmpl_t_stop  = 10.0; //ns
  //
  _fadc_offset = fadc_offset;
  _fadc_sample_in_ns = fadc_sample_in_ns;  
  _NGB_rate_in_MHz = NGB_rate_in_MHz;
  //
  _nn_fadc_point = nn_fadc_point;
  _nn_PMT_channels = nn_PMT_channels;
  //
  _t_left_in_ns = -fadc_sample_in_ns*30.0;
  _t_righ_in_ns =  fadc_sample_in_ns*_nn_fadc_point;
  //
  _n_NGB_pe_average_in_window = (_t_righ_in_ns - _t_left_in_ns)*1.0e-9*_NGB_rate_in_MHz*1.0e+6;
  //
  _dt_arr_wf_tmpl = _fadc_sample_in_ns/100.0; //ns
  //
  _gr_wf_tmpl = new TGraph();
  _gr_wf_ampl = new TGraph();
  _h1_wf_ampl = new TH1D();
  getWF_ampl(spe_dat, _Ampl_Prompt_max, _Prompt_max);
  getWF_tmpl(wf_tamplete);
  generateWF_tmpl_array();
  setupe_ADC_ampl();
  //
  _h1_adc_NGB_pedestal = new TH1D("_h1_adc_NGB_pedestal","_h1_adc_NGB_pedestal",200,-100,100);
  _h1_dadc_NGB_pedestal = new TH1D("_h1_dadc_NGB_pedestal","_h1_dadc_NGB_pedestal",200,-100,100);
  calculate_pedestal(_h1_adc_NGB_pedestal,_h1_dadc_NGB_pedestal);
  //calculate_pedestal();
  //
}

wfCamSim::~wfCamSim(){
}

void wfCamSim::generateWF_tmpl_array(){
  //std::cout<<"wfCamSim::generateWF_tmpl_array"<<std::endl;  
  Float_t safety_factor = 1.05; 
  _wf_tmpl_t_left = -(_t_righ_in_ns -_t_left_in_ns)*safety_factor;
  _wf_tmpl_t_right = (_t_righ_in_ns -_t_left_in_ns)*safety_factor;
  _wf_tmpl_t_left = ((Int_t)(_wf_tmpl_t_left/_dt_arr_wf_tmpl))*_dt_arr_wf_tmpl;
  _wf_tmpl_t_right = ((Int_t)(_wf_tmpl_t_right/_dt_arr_wf_tmpl))*_dt_arr_wf_tmpl;
  _n_arr_wf_tmpl = (_wf_tmpl_t_right - _wf_tmpl_t_left)/_dt_arr_wf_tmpl;
  _wf_tmpl_t_right = _n_arr_wf_tmpl*_dt_arr_wf_tmpl + _wf_tmpl_t_left;
  _arr_wf_tmpl = new Float_t[_n_arr_wf_tmpl];
  Float_t t_current;
  for(unsigned int i = 0;i<_n_arr_wf_tmpl;i++){
    t_current = _wf_tmpl_t_left + _dt_arr_wf_tmpl*i;
    if(t_current<=_wf_tmpl_t_start)
      _arr_wf_tmpl[i] = 0.0;
    else if (t_current>_wf_tmpl_t_start && t_current<=_wf_tmpl_t_stop)
      _arr_wf_tmpl[i] = _gr_wf_tmpl->Eval(t_current);
    else if (t_current>_wf_tmpl_t_stop)
      _arr_wf_tmpl[i] = 0.0;
  }
}

void wfCamSim::get_gr_WF_tmpl_array(TGraph *gr){
  for(unsigned int i = 0;i<_n_arr_wf_tmpl;i++)
    gr->SetPoint( i, (_wf_tmpl_t_left + _dt_arr_wf_tmpl*i), _arr_wf_tmpl[i]);
}

Double_t wfCamSim::integrate_spe( TGraph *gr, Double_t x0, Double_t Dx){
  Double_t gr_integral = 0.0;
  Double_t x, y;
  Int_t nn = 1000;
  Double_t d_x = Dx/nn;
  for(Int_t i = 0;i<nn;i++){
    x = x0 + d_x*i;
    y = gr->Eval(x);
    gr_integral += y;
  }
  return gr_integral*d_x;
}

Int_t wfCamSim::generate_single_pe_amplitude_invf(){
  Double_t val = _rnd->Uniform(0.0,1.0);
  if(val >= _h1_wf_ampl_ADC_int_arr_min[0] && val < _h1_wf_ampl_ADC_int_arr_max[0])
    return 1;
  else if(val >= _h1_wf_ampl_ADC_int_arr_min[1] && val < _h1_wf_ampl_ADC_int_arr_max[1])
    return 2;
  else if(val >= _h1_wf_ampl_ADC_int_arr_min[2] && val < _h1_wf_ampl_ADC_int_arr_max[2])
    return 3;
  else if(val >= _h1_wf_ampl_ADC_int_arr_min[3] && val < _h1_wf_ampl_ADC_int_arr_max[3])
    return 4;
  else if(val >= _h1_wf_ampl_ADC_int_arr_min[4] && val < _h1_wf_ampl_ADC_int_arr_max[4])
    return 5;
  else if(val >= _h1_wf_ampl_ADC_int_arr_min[5] && val < _h1_wf_ampl_ADC_int_arr_max[5])
    return 6;
  else if(val >= _h1_wf_ampl_ADC_int_arr_min[6] && val < _h1_wf_ampl_ADC_int_arr_max[6])
    return 7;
  else if(val >= _h1_wf_ampl_ADC_int_arr_min[7] && val < _h1_wf_ampl_ADC_int_arr_max[7])
    return 8;
  else if(val >= _h1_wf_ampl_ADC_int_arr_min[8] && val < _h1_wf_ampl_ADC_int_arr_max[8])
    return 9;
  else if(val >= _h1_wf_ampl_ADC_int_arr_min[9] && val < _h1_wf_ampl_ADC_int_arr_max[9])
    return 10;
  else if(val >= _h1_wf_ampl_ADC_int_arr_min[10] && val < _h1_wf_ampl_ADC_int_arr_max[10])
    return 11;
  //
  Int_t i = 11;
  while(i<__n_ADC_bins){
    if(val >= _h1_wf_ampl_ADC_int_arr_min[i] && val < _h1_wf_ampl_ADC_int_arr_max[i])
      return i+1;
    i++;
  }
  assert(0);
  return 0;   
}

Int_t wfCamSim::generate_single_pe_amplitude(){
  return (Int_t)_ampl_ADC_arr[((unsigned int)_rnd->Uniform(0.0,_n_ADC_max_for_generator))];
} 

Int_t wfCamSim::generate_single_pe_amplitude_from_hist(){
  Double_t r_y;
  Int_t    r_x;
  //if(_h1_wf_ampl_ADC->GetMaximum() <= 0){
  //std::cout<<" ERROR--> : _h1_wf_ampl_ADC->GetMaximum() <= 0 "<<std::endl
  //         <<"            _h1_wf_ampl_ADC->GetMaximum()  =   "<<_h1_wf_ampl_ADC->GetMaximum()<<std::endl;
  //assert(0);
  //}
  while(1){
    r_y = _rnd->Uniform(0.0,_h1_wf_ampl_ADC->GetMaximum());
    r_x = (Int_t)_rnd->Uniform(1,_h1_wf_ampl_ADC->GetNbinsX());
    if(r_y<=_h1_wf_ampl_ADC->GetBinContent(r_x))
      return _h1_wf_ampl_ADC->GetBinCenter(r_x);
  }
  return 0.0;
}

void wfCamSim::test_single_pe_amplitude_invf_generate(TH1D *h1, Int_t n_pe_to_sim){
  clock_t start, finish;
  h1->SetBins(_h1_wf_ampl_ADC->GetNbinsX(),
	      _h1_wf_ampl_ADC->GetBinLowEdge(1),
	      _h1_wf_ampl_ADC->GetBinLowEdge(_h1_wf_ampl_ADC->GetNbinsX()) + _h1_wf_ampl_ADC->GetBinWidth(_h1_wf_ampl_ADC->GetNbinsX()));
  start = clock();
  for(Int_t i = 0;i<n_pe_to_sim;i++)
    h1->Fill(((Double_t)generate_single_pe_amplitude_invf()-0.1),1.0/n_pe_to_sim);
  finish = clock();
  std::cout<<"test_single_pe_amplitude_invf_generate : "<<std::endl
	   <<"function time + fill histogram         : "<<((finish - start)/CLOCKS_PER_SEC)<<" (sec)"<<std::endl;
  start = clock();
  for(Int_t i = 0;i<n_pe_to_sim;i++)
    generate_single_pe_amplitude_invf();
  finish = clock();
  std::cout<<"function time                          : "<<((finish - start)/CLOCKS_PER_SEC)<<" (sec)"<<std::endl;
}

void wfCamSim::test_single_pe_amplitude_generator( TH1D *h1, Int_t n_pe_to_sim){
  clock_t start, finish;
  h1->SetBins(_h1_wf_ampl_ADC->GetNbinsX(),
	      _h1_wf_ampl_ADC->GetBinLowEdge(1),
	      _h1_wf_ampl_ADC->GetBinLowEdge(_h1_wf_ampl_ADC->GetNbinsX()) + _h1_wf_ampl_ADC->GetBinWidth(_h1_wf_ampl_ADC->GetNbinsX()));
  //for(unsigned int i = 0;i<((unsigned int)_n_ADC_max_for_generator);i++)
  //h1->Fill(((Double_t)_ampl_ADC_arr[i]-0.1),1.0/_n_ADC_max_for_generator);
  //_ampl_ADC_arr[((unsigned int)_rnd->Uniform(0.0,_n_ADC_max_for_generator))];
  start = clock();
  for(Int_t i = 0;i<n_pe_to_sim;i++)
    h1->Fill(((Double_t)generate_single_pe_amplitude()-0.1),1.0/n_pe_to_sim);
  finish = clock();
  std::cout<<"test_single_pe_amplitude_generator : "<<std::endl
	   <<"function time + fill histogram         : "<<((finish - start)/CLOCKS_PER_SEC)<<" (sec)"<<std::endl;
  start = clock();
  for(Int_t i = 0;i<n_pe_to_sim;i++)
    generate_single_pe_amplitude();
  finish = clock();
  std::cout<<"function time                          : "<<((finish - start)/CLOCKS_PER_SEC)<<" (sec)"<<std::endl;
}

void wfCamSim::test_single_pe_amplitude_from_hist_generator( TH1D *h1, Int_t n_pe_to_sim){
  clock_t start, finish;
  h1->SetBins(_h1_wf_ampl_ADC->GetNbinsX(),
	      _h1_wf_ampl_ADC->GetBinLowEdge(1),
	      _h1_wf_ampl_ADC->GetBinLowEdge(_h1_wf_ampl_ADC->GetNbinsX()) + _h1_wf_ampl_ADC->GetBinWidth(_h1_wf_ampl_ADC->GetNbinsX()));
  start = clock();
  for(Int_t i = 0;i<n_pe_to_sim;i++)
    h1->Fill(((Double_t)generate_single_pe_amplitude_from_hist()),1.0/n_pe_to_sim);
  finish = clock();
  std::cout<<"test_single_pe_amplitude_from_hist_generator : "<<std::endl
	   <<"function time + fill histogram         : "<<((finish - start)/CLOCKS_PER_SEC)<<" (sec)"<<std::endl;
  start = clock();
  for(Int_t i = 0;i<n_pe_to_sim;i++)
    generate_single_pe_amplitude_from_hist();
  finish = clock();
  std::cout<<"function time                          : "<<((finish - start)/CLOCKS_PER_SEC)<<" (sec)"<<std::endl;
}

void wfCamSim::generateNGB(std::vector<int> &wf){
  Int_t n_pe = _rnd->Poisson(_n_NGB_pe_average_in_window);
  for( Int_t j = 0; j < n_pe; j++)
    generate_wf(wf,(Float_t)_rnd->Uniform( _t_left_in_ns, _t_righ_in_ns));
}

void wfCamSim::generate_zero_wf(std::vector<int> &wf){
  generate_zero_wf(wf, 0.0);
}

void wfCamSim::generate_zero_wf(std::vector<int> &wf, Int_t pedestal){
  for( unsigned int i = 0; i < wf.size(); i++)
    wf.at(i) = pedestal;
}
  
void wfCamSim::generate_wf(std::vector<int> &wf, Float_t pe_time){
  unsigned int i0 = (-pe_time + 1.0 - _wf_tmpl_t_left)/_dt_arr_wf_tmpl;
  Int_t ampl_single_pe = generate_single_pe_amplitude_invf();
  for( unsigned int i = 0; i < wf.size(); i++){
    unsigned int jj = i0+i*100;
    if(jj>=0 &&jj<_n_arr_wf_tmpl){
      wf.at(i) += _arr_wf_tmpl[jj]*ampl_single_pe;
      if(wf.at(i) > 16384)
	wf.at(i) = 16384;
    }
  }
}

void wfCamSim::generate_wf_from_gr(std::vector<int> &wf, Float_t pe_time){
  Int_t ampl_single_pe = generate_single_pe_amplitude_invf();
  for( unsigned int i = 0; i < wf.size(); i++){
    wf.at(i) += (Int_t)(_gr_wf_tmpl->Eval(pe_time-i*_fadc_sample_in_ns+5)*ampl_single_pe);
    if(wf.at(i)>16384)
      wf.at(i) = 16384;
  }
}

void wfCamSim::calculate_pedestal(){
  calculate_pedestal(NULL, NULL);
}

void wfCamSim::calculate_pedestal(TH1D *h1_adc, TH1D *h1_dadc){
  std::vector<std::vector<Int_t>> wf_pedestal(_nn_PMT_channels, std::vector<Int_t>(_nn_fadc_point));
  for( unsigned int i = 0; i < wf_pedestal.size(); i++){
    generate_zero_wf(wf_pedestal.at(i));
    generateNGB(wf_pedestal.at(i));
  }
  calculate_camera_ADC_mean_and_std( wf_pedestal, _NGB_pedestal_mean, _NGB_pedestal_std, h1_adc, h1_dadc);
}
  
void wfCamSim::calculate_camera_ADC_mean_and_std(const std::vector<std::vector<Int_t>> &wf, Float_t &meanv, Float_t &stdv){
  calculate_camera_ADC_mean_and_std(wf, meanv, stdv, NULL, NULL);
}

void wfCamSim::calculate_camera_ADC_mean_and_std(const std::vector<std::vector<Int_t>> &wf, Float_t &meanv, Float_t &stdv, TH1D *h1_adc, TH1D *h1_dadc){
  meanv = 0.0;
  stdv = 0.0;
  for( unsigned int i = 0; i < wf.size(); i++){
    for( unsigned int j = 0; j < wf.at(i).size(); j++){
      meanv += wf.at(i).at(j);
      if(h1_adc != NULL)
	h1_adc->Fill(wf.at(i).at(j));
      if(h1_dadc != NULL && j>0)
	h1_dadc->Fill(wf.at(i).at(j)-wf.at(i).at(j-1));
    }
  }
  meanv /= (wf.size()*wf.at(0).size());
  for( unsigned int i = 0; i < wf.size(); i++)
    for( unsigned int j = 0; j < wf.at(i).size(); j++)
      stdv += TMath::Power((meanv - wf.at(i).at(j)),2);
  stdv = TMath::Sqrt(stdv/(wf.size()*wf.at(0).size()));
}
  
void wfCamSim::print_wfCamSim_configure(){
  std::cout<<"_Ampl_Prompt_max            "<<_Ampl_Prompt_max<<std::endl
	   <<"_Prompt_max                 "<<_Prompt_max<<std::endl
	   <<"_t_max_ampl_wf_tmpl         "<<_t_max_ampl_wf_tmpl<<std::endl
	   <<"_fadc_offset                "<<_fadc_offset<<std::endl
	   <<"_fadc_sample_in_ns          "<<_fadc_sample_in_ns<<std::endl
    	   <<"_fadc_electronic_noise_RMS  "<<_fadc_electronic_noise_RMS<<std::endl
    	   <<"_fadc_pedestal              "<<_fadc_pedestal<<std::endl
	   <<"_NGB_rate_in_MHz            "<<_NGB_rate_in_MHz<<std::endl
	   <<"_t_left_in_ns               "<<_t_left_in_ns<<std::endl
	   <<"_t_righ_in_ns               "<<_t_righ_in_ns<<std::endl
	   <<"_n_NGB_pe_average_in_window "<<_n_NGB_pe_average_in_window<<std::endl
	   <<"_nn_fadc_point              "<<_nn_fadc_point<<std::endl
	   <<"_nn_PMT_channels            "<<_nn_PMT_channels<<std::endl;
  std::cout<<"_n_arr_wf_tmpl              "<<_n_arr_wf_tmpl<<std::endl
	   <<"_dt_arr_wf_tmpl             "<<_dt_arr_wf_tmpl<<std::endl
    	   <<"_wf_tmpl_t_start            "<<_wf_tmpl_t_start<<std::endl
    	   <<"_wf_tmpl_t_stop             "<<_wf_tmpl_t_stop<<std::endl
    	   <<"_wf_tmpl_t_left             "<<_wf_tmpl_t_left<<std::endl
  	   <<"_NGB_pedestal_mean          "<<_NGB_pedestal_mean<<std::endl
    	   <<"_NGB_pedestal_std           "<<_NGB_pedestal_std<<std::endl;
}

Int_t wfCamSim::get_charge(const std::vector<int>& wf, Int_t offset){
  Int_t charge_wf = 0;
  for( unsigned int i = 0; i < wf.size(); i++ )
    charge_wf += wf.at(i) - offset;
  return charge_wf;
}

Int_t wfCamSim::get_charge(const std::vector<int>& wf){
  return get_charge(wf, _fadc_offset);
}

void wfCamSim::generate_electronic_noise(std::vector<int> &wf){
  for( unsigned int i = 0; i < wf.size(); i++)
    wf.at(i) += (Int_t)_rnd->Gaus(0.0,_fadc_electronic_noise_RMS);
}

void wfCamSim::generate_electronic_noise_pedestal_removal(std::vector<int> &wf){
  for( unsigned int i = 0; i < wf.size(); i++)
    wf.at(i) += ((Int_t)_rnd->Gaus(0.0,_fadc_electronic_noise_RMS) + (Int_t)_rnd->Uniform(0,3.5)-1);
}

void wfCamSim::simulate_cam_event(const Int_t nn_fadc_point,
				  const Int_t nn_PMT_channels,
				  std::vector<std::vector<int>>& wf){
  for( unsigned int i = 0; i < wf.size(); i++ ){
    generate_zero_wf(wf.at(i),(_fadc_offset-_NGB_pedestal_mean));
    //generate_zero_wf(wf.at(i),_fadc_offset);
    generateNGB(wf.at(i));
    if(_electronic_noise_with_pedestal_removal)
      generate_electronic_noise_pedestal_removal(wf.at(i));
    else
      generate_electronic_noise(wf.at(i));
  }
}

void wfCamSim::simulate_cam_event_NGB(std::vector<std::vector<Int_t>> &wf){
  for( unsigned int i = 0; i < wf.size(); i++ ){
    generate_zero_wf(wf.at(i),(_fadc_offset-_NGB_pedestal_mean));
    //generate_zero_wf(wf.at(i),_fadc_offset);
    generateNGB(wf.at(i));
    generate_electronic_noise_pedestal_removal(wf.at(i));
  }
}

void wfCamSim::simulate_cam_event(const Int_t nn_fadc_point,
				  const Int_t nn_PMT_channels,
				  std::vector<std::vector<Int_t>> &wf,
				  const Float_t ev_time,
				  const Float_t time_offset,
				  const Int_t n_pe,
				  const Int_t *pe_chID,
				  const Float_t *pe_time){
  simulate_cam_event( nn_fadc_point, nn_PMT_channels, wf); 
  for(Int_t i = 0;i<n_pe;i++)
    if(check_pe_chID(pe_chID[i]))
      generate_wf(wf.at((unsigned int)pe_chID[i]), pe_time[i] - ev_time + time_offset);
}

bool wfCamSim::check_pe_chID(Int_t id_ch){
  if( (id_ch<0) || ((unsigned int)id_ch>= _nn_PMT_channels)){
    std::cout<<" ERROR --> id_ch<0 || id_ch >= _nn_PMT_channels"<<std::endl
	     <<"                       id_ch = "<<id_ch<<std::endl
      	     <<"            _nn_PMT_channels = "<<_nn_PMT_channels<<std::endl;
    assert(0);
    return false;
  }
  return true;
}
  
void wfCamSim::setupe_ADC_ampl(Double_t bits_per_pe, Int_t n_max_pe){
  //
  if(_h1_wf_ampl_ADC == NULL){
    Int_t n_ADC_bins = (n_max_pe+(floor(n_max_pe*bits_per_pe+0.5) - n_max_pe*bits_per_pe)/bits_per_pe)*bits_per_pe;
    //
    std::cout<<"n_max_pe    "<<n_max_pe<<std::endl
	     <<"n_ADC_bins  "<<n_ADC_bins<<std::endl
	     <<"bits_per_pe "<<bits_per_pe<<std::endl;
    if(n_ADC_bins>255){
      std::cout<<"  ---> ERROR : n_ADC_bins>255 "<<n_ADC_bins<<std::endl;
      assert(0);
    }
    //std::cout<<"_h1_wf_ampl_ADC == NULL"<<std::endl;
    //
    _h1_wf_ampl_ADC = new TH1D("_h1_wf_ampl_ADC","_h1_wf_ampl_ADC",n_ADC_bins,0.0,n_ADC_bins);
    _h1_wf_ampl_ADC_integral = new TH1D("_h1_wf_ampl_ADC_integral","_h1_wf_ampl_ADC_integral",n_ADC_bins,0.0,n_ADC_bins);
    //
    _h1_wf_ampl_ADC_int_arr_max = new Double_t[n_ADC_bins];
    _h1_wf_ampl_ADC_int_arr_min = new Double_t[n_ADC_bins];
    __n_ADC_bins = n_ADC_bins;
    for(Int_t i = 1;i<=n_ADC_bins;i++){
      _h1_wf_ampl_ADC->SetBinContent(i,integrate_spe(_gr_wf_ampl,
						     _h1_wf_ampl_ADC->GetBinLowEdge(i)/bits_per_pe,
						     _h1_wf_ampl_ADC->GetBinWidth(i)/bits_per_pe));
      _h1_wf_ampl_ADC_integral->SetBinContent(i,_h1_wf_ampl_ADC_integral->GetBinContent(i-1)+_h1_wf_ampl_ADC->GetBinContent(i));
      if(i==1){
	_h1_wf_ampl_ADC_int_arr_max[i-1] = 0.0 + _h1_wf_ampl_ADC->GetBinContent(i);
      }
      else if (i == n_ADC_bins){
	_h1_wf_ampl_ADC_int_arr_max[i-1] = 1.0;
      }
      else{
	_h1_wf_ampl_ADC_int_arr_max[i-1] = _h1_wf_ampl_ADC_int_arr_max[i-1-1] + _h1_wf_ampl_ADC->GetBinContent(i);
      }
    }
    //
    for(Int_t i = 0;i<n_ADC_bins;i++){
      if(i==0)
	_h1_wf_ampl_ADC_int_arr_min[i] = 0.0;
      else
	_h1_wf_ampl_ADC_int_arr_min[i] = _h1_wf_ampl_ADC_int_arr_max[i-1];
    }    
    //
    //for(Int_t i = 0;i<n_ADC_bins;i++)
    //std::cout<<std::setw(30)<<i
    //	       <<std::setw(30)<<std::setprecision(20)<<_h1_wf_ampl_ADC_int_arr_min[i]
    //	       <<std::setw(30)<<std::setprecision(20)<<_h1_wf_ampl_ADC_int_arr_max[i]<<std::endl;
    //
    _ampl_ADC_arr = new unsigned char[_n_ADC_arr];
    unsigned int counter_per_bin = 0;
    unsigned int totcounter = 0;
    for( Int_t i = 1; i <= n_ADC_bins; i++){
      counter_per_bin = (unsigned int)(_h1_wf_ampl_ADC->GetBinContent(i)*_n_ADC_arr);
      //std::cout<<i<<" "<<_h1_wf_ampl_ADC->GetBinContent(i)<<" "<<counter_per_bin<<std::endl;
      for( unsigned int j = 0; j < counter_per_bin; j++)
	_ampl_ADC_arr[totcounter+j] = (unsigned char)i;
      totcounter += counter_per_bin;
    }
    //
    _n_ADC_max_for_generator = totcounter;
    std::cout<<"_n_ADC_arr               : "<<_n_ADC_arr<<std::endl
	     <<"totcounter               : "<<totcounter<<std::endl
      	     <<"_n_ADC_max_for_generator : "<<(Int_t)_n_ADC_max_for_generator<<std::endl;
    //
    if(totcounter<_n_ADC_arr)
      for( unsigned int j = 0; j < (_n_ADC_arr - totcounter); j++)
	_ampl_ADC_arr[totcounter + j] = (unsigned char)0;
  }
  else{
    std::cout<<"_h1_wf_ampl_ADC != NULL"<<std::endl
    	       <<"_ampl_ADC_arr   != NULL"<<std::endl;
    delete _h1_wf_ampl_ADC;
    delete _ampl_ADC_arr;
    _h1_wf_ampl_ADC = NULL;
    _ampl_ADC_arr = NULL;
    setupe_ADC_ampl( bits_per_pe, n_max_pe);
  }
}

void wfCamSim::getWF_ampl(TString name, Double_t &Ampl_Prompt_max, Double_t &Prompt_max){
  //
  _gr_wf_ampl->SetNameTitle("_gr_wf_ampl","_gr_wf_ampl");
  _h1_wf_ampl->SetNameTitle("_h1_wf_ampl","_h1_wf_ampl");
  //
  std::cout<<"     spe file : "<<name<<std::endl;
  std::ifstream fFile(name.Data());
  std::string mot;
  Double_t Ampl;
  Double_t Prompt;
  Ampl_Prompt_max = 0.0;
  Prompt_max = 0.0;
  Double_t Prompt_x_AP;
  Double_t Ampl_min;
  Double_t Ampl_max;
  Double_t d_Ampl;
  //
  if(fFile.is_open()){
    while(fFile>>mot){
      if(mot == "AP)")
        break;
    }
    fFile>>mot>>mot>>mot>>mot;
    while(fFile>>Ampl>>Prompt>>Prompt_x_AP){
      _gr_wf_ampl->SetPoint(_gr_wf_ampl->GetN(),Ampl,Prompt);
      if(Prompt_max<Prompt){
        Prompt_max = Prompt;
        Ampl_Prompt_max = Ampl;
      }
    }
    fFile.close();
  }
  _gr_wf_ampl->GetPoint(0,Ampl_min,Prompt);
  _gr_wf_ampl->GetPoint((_gr_wf_ampl->GetN()-1),Ampl_max,Prompt);
  d_Ampl = (Ampl_max - Ampl_min)/(_gr_wf_ampl->GetN()-1);
  //
  _h1_wf_ampl->SetBins(_gr_wf_ampl->GetN(),Ampl_min-d_Ampl/2.0,Ampl_max+d_Ampl/2.0);
  //
  for(Int_t i = 0;i<_gr_wf_ampl->GetN();i++){
    _gr_wf_ampl->GetPoint(i,Ampl,Prompt);
    _h1_wf_ampl->SetBinContent(i+1,Prompt);
  }  
}

void wfCamSim::getWF_tmpl(TString name){
  //
  Double_t stretchingFactor = 1.0;
  //
  std::ifstream fileIn(name.Data());
  std::cout<<"template file : "<<name<<std::endl;
  double x;
  double x_min;
  double x_max;
  double y;
  double ymax = -999.0;
  TGraph *grtmp = new TGraph();
  TGraph *grtmpstrached = new TGraph();
  _gr_wf_tmpl->SetName("_gr_wf_tmpl");
  _gr_wf_tmpl->SetTitle("_gr_wf_tmpl");
  if (fileIn.is_open()){
    while(fileIn>>x>>y){
      if(grtmp->GetN() == 0)
	x_min = x;
      x_max = x;
      grtmp->SetPoint(grtmp->GetN(),x,y);
      grtmpstrached->SetPoint(grtmpstrached->GetN(),x*stretchingFactor,y);
      if(ymax<y){
	ymax = y;
	_t_max_ampl_wf_tmpl = x;
      }
    }
    fileIn.close();
  }
  else {
    std::cout<<"Unable to open file"<<std::endl;
    assert(0);
  }
  //
  for(Int_t i = 0;i<grtmp->GetN();i++){
    grtmp->GetPoint(i,x,y);
    _gr_wf_tmpl->SetPoint(i,x,grtmpstrached->Eval(x));
  }
}

Double_t wfCamSim::generate_wf_ampl_from_file(){
  Double_t r_y;
  Int_t    r_x;
  if(_Ampl_Prompt_max <= 0){
    std::cout<<" ERROR--> : _Ampl_Prompt_max <= 0 "<<std::endl
	     <<"            _Ampl_Prompt_max  =   "<<_Ampl_Prompt_max<<std::endl;
    assert(0);
  }
  while(1){
    r_y=_rnd->Uniform(0.0,_Prompt_max);
    r_x=_rnd->Uniform(1,_h1_wf_ampl->GetNbinsX());
    if(r_y<=_h1_wf_ampl->GetBinContent(r_x))
      return _h1_wf_ampl->GetBinCenter(r_x)*_Ampl_Prompt_max;
  }
  return 0.0;
}

void wfCamSim::generate_gif_for_event(TString pathPref, Int_t event_id, const std::vector<std::vector<Int_t>> &wf){
  std::vector<std::vector<Int_t>> wf_ref;
  std::vector<std::vector<unsigned int>> trg_vector;
  generate_gif_for_event(pathPref, event_id, wf, wf_ref,trg_vector);
}

void wfCamSim::generate_gif_for_event(TString pathPref, Int_t event_id,
				      const std::vector<std::vector<Int_t>> &wf,
				      const std::vector<std::vector<Int_t>> &wf_ref,
				      const std::vector<std::vector<unsigned int>> &trg_vector,
				      const anabase *ab){
  //
  bool if_all_trig_superimposed = false;
  bool if_clean_with_trg_vector = true;
  //
  //gStyle->SetPalette(kCherry);
  //TColor::InvertPalette();
  //
  std::vector<unsigned int> trg_vector_all;
  for(Int_t i = 0;i<wf.at(0).size();i++){
    for(Int_t k = 0;k<trg_vector.at(i).size();k++){
      trg_vector_all.push_back(trg_vector.at(i).at(k));
    }
  }
  //
  Bool_t pdf_out = true;
  //Bool_t pdf_out = false;
  //
  TString ev_dir_name = pathPref;
  ev_dir_name += event_id; 
  mkdir(ev_dir_name.Data(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  //std::vector<sipmCameraHist*> sipm_cam_v; 
  std::ofstream merge_gif;
  TString merge_gif_name=ev_dir_name;
  merge_gif_name += "/mrg.sh";
  //cout<<"merge_gif_name "<<merge_gif_name<<endl;
  merge_gif.open(merge_gif_name.Data());
  //
  merge_gif<<"convert -delay 15 -loop 1000 ";
  //for(Int_t i = 20;i<(nn_fadc_point-20);i++){
  Int_t nn_fadc_point = wf.at(0).size();
  Int_t nChannels =wf.size();
  sipmCameraHist *sipm_cam = new sipmCameraHist("sipm_cam","sipm_cam","pixel_mapping.csv",0);
  sipm_cam->SetMinimum(299.0);
  //sipm_cam->SetMaximum(500.0);
  sipm_cam->SetMaximum(308.0);
  sipmCameraHist *sipm_cam_ref = new sipmCameraHist("sipm_cam_ref","sipm_cam_ref","pixel_mapping.csv",0);
  sipm_cam_ref->SetMinimum(299.0);
  //sipm_cam_ref->SetMaximum(500.0);
  sipm_cam_ref->SetMaximum(308.0);
  for(Int_t i = 0;i<nn_fadc_point;i++){
    TString sipm_cam_name = "sipm_cam_";
    TString gif_name = ev_dir_name;
    TString pdf_name = ev_dir_name;
    gif_name += "/sipm_cam_";
    gif_name += i;
    gif_name += ".gif";
    //
    pdf_name += "/sipm_cam_";
    pdf_name += i;
    pdf_name += ".pdf";
    TString gif_name_short = "sipm_cam_";
    gif_name_short += i;
    gif_name_short += ".gif";
    //gif_name_short += ".pdf";
    sipm_cam_name += i;
    //sipmCameraHist *sipm_cam = new sipmCameraHist(sipm_cam_name.Data(),sipm_cam_name.Data(),"pixel_mapping.csv",0);
    //sipm_cam->SetMinimum(0.0);
    //sipm_cam->SetMaximum(TMath::Power(2,14));
    for(Int_t j = 0;j<nChannels;j++){
      if(if_clean_with_trg_vector){
	sipm_cam->SetBinContent(j+1,300);
	for(Int_t kk = 0;kk<trg_vector.at(i).size();kk++){
	  if(trg_vector.at(i).at(kk)==j)
	    sipm_cam->SetBinContent(j+1,wf.at(j).at(i));
	}
	sipm_cam_ref->SetBinContent(j+1,wf_ref.at(j).at(i));
      }
      else{
	sipm_cam->SetBinContent(j+1,wf.at(j).at(i));
	sipm_cam_ref->SetBinContent(j+1,wf_ref.at(j).at(i));
      }
    }
    merge_gif<<gif_name_short<<" ";
    //
    //std::vector<Int_t> pixel_line_flower_vec;
    //for(unsigned int i = 0;i<trg_vector.size();i++)
    //std::cout<<<<" "<<trg_vector.at(i).at(1)<<std::endl;
    //
    //sipm_cam->Draw_cam("ZCOLOR",gif_name.Data(),"gamma",i,event_id,energy,xcore,ycore,ev_time,nphotons,n_pe,n_pixels);
    //sipm_cam->Draw_cam("ZCOLOR",gif_name.Data(),"proton",i,event_id,energy,xcore,ycore,ev_time,nphotons,n_pe,n_pixels);
    sipm_cam->set_wf_time_id(i);
    sipm_cam->set_anabase(ab);
    if(trg_vector.size()>0){
      //if(trg_vector.at(i).size()>0){
      if(trg_vector_all.size()>0){
	//std::cout<<"rr"<<std::endl;
	if(!if_all_trig_superimposed)
	  sipm_cam->Draw_cam("ZCOLOR",gif_name.Data(),sipm_cam_ref,trg_vector.at(i), ab);
	else
	  sipm_cam->Draw_cam("ZCOLOR",gif_name.Data(),sipm_cam_ref,trg_vector_all, ab);
	if(pdf_out)
	  sipm_cam->Draw_cam("ZCOLOR",pdf_name.Data(),sipm_cam_ref,trg_vector.at(i), ab);
      }
      else{
	sipm_cam->Draw_cam("ZCOLOR",gif_name.Data(), sipm_cam_ref, ab);
	if(pdf_out)
	  sipm_cam->Draw_cam("ZCOLOR",pdf_name.Data(), sipm_cam_ref, ab);
      }
    }
    else{
      //sipm_cam->Draw_cam("ZCOLOR",gif_name.Data(), sipm_cam_ref, ab);
      if(pdf_out)
	sipm_cam->Draw_cam("ZCOLOR",pdf_name.Data(), sipm_cam_ref, ab);
    }
    //
    //delete sipm_cam;
  }
  TString outtotGifName = "sipm_cam_ev";
  outtotGifName += event_id;
  outtotGifName += ".gif";
  merge_gif<<  outtotGifName.Data();
  merge_gif.close();
}
