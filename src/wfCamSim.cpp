//my
#include "wfCamSim.hh"

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

wfCamSim::wfCamSim( TRandom3 *rnd, TString wf_tamplete, TString spe_dat){
  _rnd = rnd;
  _wf_tamplete = wf_tamplete;
  _spe_dat = spe_dat;
  _gr_wf_tmpl = new TGraph();
  _gr_wf_ampl = new TGraph();
  _h1_wf_ampl = new TH1D();
  getWF_ampl(spe_dat, _Ampl_Prompt_max, _Prompt_max);
  getWF_tmpl(wf_tamplete);
  setupe_ADC_ampl();
}

wfCamSim::~wfCamSim(){
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
				  std::vector<std::vector<int>>& wf,
				  const Float_t NGB_rate_in_MHz){
  for(unsigned int i = 0; i < wf.size(); i++)
    for(unsigned int j = 0; j < wf.at(i).size(); j++)
      wf.at(i).at(j) = 0;
}

void wfCamSim::simulate_cam_event(const Int_t nn_fadc_point,
				  const Int_t nn_PMT_channels,
				  std::vector<std::vector<Int_t>> &wf,
				  const Float_t NGB_rate_in_MHz,
				  const Float_t ev_time,
				  const Float_t time_offset,
				  const Int_t n_pe,
				  const Int_t *pe_chID,
				  const Float_t *pe_time){
  simulate_cam_event( nn_fadc_point, nn_PMT_channels, wf, NGB_rate_in_MHz); 
  for(Int_t i = 0;i<n_pe;i++){
    std::cout<<pe_chID[i]<<std::endl
	     <<pe_time[i]<<std::endl;
  }
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
    //std::cout<<"_h1_wf_ampl_ADC != NULL"<<std::endl
    //	       <<"_ampl_ADC_arr   != NULL"<<std::endl;
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
