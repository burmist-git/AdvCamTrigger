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
  getWF_ampl(spe_dat, _Ampl_Prompt_max, _Prompt_max);
}

wfCamSim::~wfCamSim(){
}

void wfCamSim::getWF_ampl(TString name, Double_t &Ampl_Prompt_max, Double_t &Prompt_max){
  //
  _gr_wf_ampl->SetNameTitle("_gr_wf_ampl","_gr_wf_ampl");
  _h1_wf_ampl->SetNameTitle("_h1_wf_ampl","_h1_wf_ampl");
  //
  std::cout<<name<<std::endl;
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

/*
void wfCamSim::gen_WF( TGraph *gr_wf, TGraph *gr_wf_sig, TGraph *gr_wf_sig_only, unsigned int n_signals, TH1D *h1_photon_time){
  //
  Double_t dT = (_wfConf->SimEndTime - _wfConf->SimStartTime);
  Double_t dT_s = dT*1.0e-9; 
  //unsigned int n_noise = (Int_t)(TMath::Floor(dT_s*_wfConf->DCRrate)+1);
  //Int_t n = (Int_t)(TMath::Floor(dT/_wfConf->SamplingTime)+1);
  unsigned int n_noise = (Int_t)(TMath::Floor(dT_s*_wfConf->DCRrate));
  Int_t n = (Int_t)(TMath::Floor(dT/_wfConf->SamplingTime));
  //
  Double_t first_pe_time;
  Double_t first_pe_ampl;
  Int_t parentID;
  Int_t typeID;
  Int_t parentGenerationID;
  Double_t probabilityCorrectionFactor = 1.0;
  //noise
  std::vector<std::vector<photoElectronInfo>> all_pe_vec;
  for(unsigned int i = 0; i<n_noise; i++){
    std::vector<photoElectronInfo> pe_vec;
    first_pe_time = _rnd->Uniform(_wfConf->SimEndTime,_wfConf->SimStartTime);
    if(_wfConf->amplDistFile != "NONE"){
      first_pe_ampl = generate_wf_ampl_from_file();
    }
    else{
      first_pe_ampl = _rnd->Gaus(_wfConf->single_p_e_ampl,_wfConf->SigmaGain);
    }
    _h1_first_pe_ampl->Fill(first_pe_ampl);
    parentGenerationID = -1;
    parentID = -1;
    typeID = 0;
    simPE(pe_vec, first_pe_time, first_pe_ampl, typeID, parentGenerationID, parentID, probabilityCorrectionFactor);
    all_pe_vec.push_back(pe_vec);
  }
  //print_pe_vec(all_pe_vec);
  //signal
  std::vector<std::vector<photoElectronInfo>> sig_pe_vec;
  Double_t dT_shower;
  //
  for(unsigned int i = 0; i<n_signals; i++){
    std::vector<photoElectronInfo> pe_vec;
    dT_shower = generateDistFromHist(h1_photon_time);
    //std::cout<<"dT_shower = "<<dT_shower<<std::endl;
    first_pe_time = _wfConf->signal_t0 + dT_shower;
    first_pe_ampl = _rnd->Gaus(_wfConf->single_p_e_ampl,_wfConf->SigmaGain);
    parentGenerationID = -1;
    parentID = -1;
    typeID = 0;
    simPE(pe_vec, first_pe_time, first_pe_ampl, typeID, parentGenerationID, parentID, probabilityCorrectionFactor);
    sig_pe_vec.push_back(pe_vec);
  }
  //
  //print_pe_vec(sig_pe_vec);
  //  
  Double_t t;
  Double_t ampl;
  Double_t ampl_p_signal;
  Double_t ampl_signal;
  Double_t electricNoise;
  for(Int_t i = 0;i<n;i++){
    //if(i%10000==0)
    //std::cout<<i<<std::endl;
    t = _wfConf->SimStartTime + _wfConf->SamplingTime*i;
    ampl = 0.0;
    ampl_p_signal = 0.0;
    ampl_signal = 0.0;
    //p.e. noise
    for(unsigned int j = 0;j<all_pe_vec.size();j++){
      for(unsigned int k = 0;k<all_pe_vec.at(j).size();k++){
	if((t - all_pe_vec.at(j).at(k).time )>=0 && (t - all_pe_vec.at(j).at(k).time)<1000){
    	  ampl += all_pe_vec.at(j).at(k).ampl*_gr_wf_tmpl->Eval(t - all_pe_vec.at(j).at(k).time);
	}
      }
    }
    ampl_p_signal = ampl;
    //signal
    for(unsigned int j = 0;j<sig_pe_vec.size();j++){
      for(unsigned int k = 0;k<sig_pe_vec.at(j).size();k++){
	if((t - sig_pe_vec.at(j).at(k).time )>=0 && (t - sig_pe_vec.at(j).at(k).time)<1000){   
	  ampl_p_signal += sig_pe_vec.at(j).at(k).ampl*_gr_wf_tmpl->Eval(t - sig_pe_vec.at(j).at(k).time);
	  ampl_signal += sig_pe_vec.at(j).at(k).ampl*_gr_wf_tmpl->Eval(t - sig_pe_vec.at(j).at(k).time);
	}
      }
    }
    //electric noise
    electricNoise = _rnd->Gaus(_wfConf->ElectronicBaseine,_wfConf->ElectronicNoiseSigm);
    ampl += electricNoise;
    ampl_p_signal += electricNoise;
    ampl_signal += electricNoise;
    gr_wf->SetPoint(gr_wf->GetN(),t,ampl);
    gr_wf_sig->SetPoint(gr_wf_sig->GetN(),t,ampl_p_signal);
    gr_wf_sig_only->SetPoint(gr_wf_sig_only->GetN(),t,ampl_signal);
  }
}
*/

/*
void wfCamSim::gen_WF( TGraph *gr_wf, TGraph *gr_wf_sig, TGraph *gr_wf_sig_only, unsigned int n_signals){
  //
  Double_t dT = (_wfConf->SimEndTime - _wfConf->SimStartTime);
  Double_t dT_s = dT*1.0e-9; 
  //unsigned int n_noise = (Int_t)(TMath::Floor(dT_s*_wfConf->DCRrate)+1);
  //Int_t n = (Int_t)(TMath::Floor(dT/_wfConf->SamplingTime)+1);
  unsigned int n_noise = (Int_t)(TMath::Floor(dT_s*_wfConf->DCRrate));
  Int_t n = (Int_t)(TMath::Floor(dT/_wfConf->SamplingTime));
  //
  Double_t first_pe_time;
  Double_t first_pe_ampl;
  Int_t parentID;
  Int_t typeID;
  Int_t parentGenerationID;
  Double_t probabilityCorrectionFactor = 1.0;
  //noise
  std::vector<std::vector<photoElectronInfo>> all_pe_vec;
  for(unsigned int i = 0; i<n_noise; i++){
    std::vector<photoElectronInfo> pe_vec;
    first_pe_time = _rnd->Uniform(_wfConf->SimEndTime,_wfConf->SimStartTime);
    first_pe_ampl = _rnd->Gaus(_wfConf->single_p_e_ampl,_wfConf->SigmaGain);
    parentGenerationID = -1;
    parentID = -1;
    typeID = 0;
    simPE(pe_vec, first_pe_time, first_pe_ampl, typeID, parentGenerationID, parentID, probabilityCorrectionFactor);
    all_pe_vec.push_back(pe_vec);
  }
  //print_pe_vec(all_pe_vec);
  //signal
  std::vector<std::vector<photoElectronInfo>> sig_pe_vec;
  for(unsigned int i = 0; i<n_signals; i++){
    std::vector<photoElectronInfo> pe_vec;
    first_pe_time = _wfConf->signal_t0;
    first_pe_ampl = _rnd->Gaus(_wfConf->single_p_e_ampl,_wfConf->SigmaGain);
    parentGenerationID = -1;
    parentID = -1;
    typeID = 0;
    simPE(pe_vec, first_pe_time, first_pe_ampl, typeID, parentGenerationID, parentID, probabilityCorrectionFactor);
    sig_pe_vec.push_back(pe_vec);
  }
  //
  //print_pe_vec(sig_pe_vec);
  //  
  Double_t t;
  Double_t ampl;
  Double_t ampl_p_signal;
  Double_t ampl_signal;
  Double_t electricNoise;
  for(Int_t i = 0;i<n;i++){
    //if(i%10000==0)
    //std::cout<<i<<std::endl;
    t = _wfConf->SimStartTime + _wfConf->SamplingTime*i;
    ampl = 0.0;
    ampl_p_signal = 0.0;
    ampl_signal = 0.0;
    //p.e. noise
    for(unsigned int j = 0;j<all_pe_vec.size();j++){
      for(unsigned int k = 0;k<all_pe_vec.at(j).size();k++){
	if((t - all_pe_vec.at(j).at(k).time )>=0 && (t - all_pe_vec.at(j).at(k).time)<1000){
    	  ampl += all_pe_vec.at(j).at(k).ampl*_gr_wf_tmpl->Eval(t - all_pe_vec.at(j).at(k).time);
	}
      }
    }
    ampl_p_signal = ampl;
    //signal
    for(unsigned int j = 0;j<sig_pe_vec.size();j++){
      for(unsigned int k = 0;k<sig_pe_vec.at(j).size();k++){
	if((t - sig_pe_vec.at(j).at(k).time )>=0 && (t - sig_pe_vec.at(j).at(k).time)<1000){   
	  ampl_p_signal += sig_pe_vec.at(j).at(k).ampl*_gr_wf_tmpl->Eval(t - sig_pe_vec.at(j).at(k).time);
	  ampl_signal += sig_pe_vec.at(j).at(k).ampl*_gr_wf_tmpl->Eval(t - sig_pe_vec.at(j).at(k).time);
	}
      }
    }
    //electric noise
    electricNoise = _rnd->Gaus(_wfConf->ElectronicBaseine,_wfConf->ElectronicNoiseSigm);
    ampl += electricNoise;
    ampl_p_signal += electricNoise;
    ampl_signal += electricNoise;
    gr_wf->SetPoint(gr_wf->GetN(),t,ampl);
    gr_wf_sig->SetPoint(gr_wf_sig->GetN(),t,ampl_p_signal);
    gr_wf_sig_only->SetPoint(gr_wf_sig_only->GetN(),t,ampl_signal);
  }  
}
*/
