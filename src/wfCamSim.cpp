//my
#include "wfCamSim.hh"
#include "sipmCameraHist.hh"

//root
#include "TRandom3.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TMath.h"
#include "TH1D.h"

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
		    const Float_t NGB_rate_in_MHz){
  //
  _h1_wf_ampl_ADC = NULL;
  _ampl_ADC_arr = NULL;
  //
  _rnd = rnd;
  _wf_tamplete = wf_tamplete;
  _spe_dat = spe_dat;
  //
  _fadc_electronic_noise_RMS = 3.94;
  _fadc_pedestal = 300;
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
  calculate_pedestal();
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

void wfCamSim::ger_gr_WF_tmpl_array(TGraph *gr){
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

void wfCamSim::test_single_pe_amplitude_generator( TH1D *h1, Int_t n_pe_to_sim){
  h1->SetBins(_h1_wf_ampl_ADC->GetNbinsX(),
	      _h1_wf_ampl_ADC->GetBinLowEdge(1),
	      _h1_wf_ampl_ADC->GetBinLowEdge(_h1_wf_ampl_ADC->GetNbinsX()) + _h1_wf_ampl_ADC->GetBinWidth(_h1_wf_ampl_ADC->GetNbinsX()));
  //for(unsigned int i = 0;i<((unsigned int)_n_ADC_max_for_generator);i++)
  //h1->Fill(((Double_t)_ampl_ADC_arr[i]-0.1),1.0/_n_ADC_max_for_generator);
  //_ampl_ADC_arr[((unsigned int)_rnd->Uniform(0.0,_n_ADC_max_for_generator))];
  for(Int_t i = 0;i<n_pe_to_sim;i++)
    h1->Fill(((Double_t)generate_single_pe_amplitude()-0.1),1.0/n_pe_to_sim);
}

void wfCamSim::test_single_pe_amplitude_from_hist_generator( TH1D *h1, Int_t n_pe_to_sim){
  h1->SetBins(_h1_wf_ampl_ADC->GetNbinsX(),
	      _h1_wf_ampl_ADC->GetBinLowEdge(1),
	      _h1_wf_ampl_ADC->GetBinLowEdge(_h1_wf_ampl_ADC->GetNbinsX()) + _h1_wf_ampl_ADC->GetBinWidth(_h1_wf_ampl_ADC->GetNbinsX()));
  for(Int_t i = 0;i<n_pe_to_sim;i++)
    h1->Fill(((Double_t)generate_single_pe_amplitude_from_hist()),1.0/n_pe_to_sim);
}

void wfCamSim::simulate_cam_event(const Int_t nn_fadc_point,
				  const Int_t nn_PMT_channels,
				  std::vector<std::vector<int>>& wf){
  for( unsigned int i = 0; i < wf.size(); i++ ){
    generate_zero_wf(wf.at(i),(_fadc_offset-_NGB_pedestal_mean));
    generateNGB(wf.at(i));
    generate_electronic_noise(wf.at(i));
  }
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
  unsigned int i0 = (pe_time - _wf_tmpl_t_left)/_dt_arr_wf_tmpl;
  for( unsigned int i = 0; i < wf.size(); i++)
    wf.at(i) += (_arr_wf_tmpl[i0+i*100]*generate_single_pe_amplitude());
}
  
void wfCamSim::generate_electronic_noise(std::vector<int> &wf){
  for( unsigned int i = 0; i < wf.size(); i++)
    wf.at(i) += (Int_t)_rnd->Gaus(0.0,_fadc_electronic_noise_RMS);
}

void wfCamSim::calculate_pedestal(){
  std::vector<std::vector<Int_t>> wf_pedestal(_nn_PMT_channels, std::vector<Int_t>(_nn_fadc_point));
  for( unsigned int i = 0; i < wf_pedestal.size(); i++){
    generate_zero_wf(wf_pedestal.at(i));
    generateNGB(wf_pedestal.at(i));
  }
  calculate_camera_ADC_mean_and_std( wf_pedestal, _NGB_pedestal_mean, _NGB_pedestal_std);
}

void wfCamSim::test_calculate_pedestal(TH1D *h1_adc, TH1D *h1_dadc){
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

void wfCamSim::simulate_cam_event(const Int_t nn_fadc_point,
				  const Int_t nn_PMT_channels,
				  std::vector<std::vector<Int_t>> &wf,
				  const Float_t ev_time,
				  const Float_t time_offset,
				  const Int_t n_pe,
				  const Int_t *pe_chID,
				  const Float_t *pe_time){
  simulate_cam_event( nn_fadc_point, nn_PMT_channels, wf); 
  //for(Int_t i = 0;i<n_pe;i++){
  //std::cout<<std::setw(20)<<i
  //	     <<std::setw(20)<<pe_chID[i]
  //	     <<std::setw(20)<<pe_time[i]<<std::endl;
  //}
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
    //
    for(Int_t i = 1;i<=n_ADC_bins;i++)
      _h1_wf_ampl_ADC->SetBinContent(i,integrate_spe(_gr_wf_ampl,
						     _h1_wf_ampl_ADC->GetBinLowEdge(i)/bits_per_pe,
						     _h1_wf_ampl_ADC->GetBinWidth(i)/bits_per_pe));
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
  std::ifstream fileIn(name.Data());
  std::cout<<"template file : "<<name<<std::endl;
  double x;
  double y;
  double ymax = -999.0;
  _gr_wf_tmpl->SetName("_gr_wf_tmpl");
  _gr_wf_tmpl->SetTitle("_gr_wf_tmpl");
  if (fileIn.is_open()){
    while(fileIn>>x>>y){
      _gr_wf_tmpl->SetPoint(_gr_wf_tmpl->GetN(),x,y);
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
  for(Int_t i = 0;i<nn_fadc_point;i++){
    TString sipm_cam_name = "sipm_cam_";
    TString gif_name = ev_dir_name;
    gif_name += "/sipm_cam_";
    gif_name += i;
    gif_name += ".gif";
    TString gif_name_short = "sipm_cam_";
    gif_name_short += i;
    gif_name_short += ".gif";
    sipm_cam_name += i;
    sipmCameraHist *sipm_cam = new sipmCameraHist(sipm_cam_name.Data(),sipm_cam_name.Data(),"pixel_mapping.csv",-10);
    sipm_cam->SetMinimum(0.0);
    sipm_cam->SetMaximum(TMath::Power(2,14));
    for(Int_t j = 0;j<nChannels;j++){
      sipm_cam->SetBinContent(j+1,wf.at(j).at(i));
    }
    merge_gif<<gif_name_short<<" ";
    //sipm_cam->Draw_cam("ZCOLOR",gif_name.Data(),"gamma",i,event_id,energy,xcore,ycore,ev_time,nphotons,n_pe,n_pixels);
    //sipm_cam->Draw_cam("ZCOLOR",gif_name.Data(),"proton",i,event_id,energy,xcore,ycore,ev_time,nphotons,n_pe,n_pixels);
    sipm_cam->Draw_cam("ZCOLOR",gif_name.Data());
    //
    delete sipm_cam;
  }
  TString outtotGifName = "sipm_cam_ev";
  outtotGifName += event_id;
  outtotGifName += ".gif";
  merge_gif<<  outtotGifName.Data();
  merge_gif.close();
}
