//my
#include "anaTrgB.hh"
#include "ana.hh"
#include "anabase.hh"
#include "sipmCameraHist.hh"
#include "wfCamSim.hh"
#include "triggerSim.hh"
#include "evstHist.hh"
#include "dbscan.hh"
#include "anaConf.hh"

//root
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TString.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TRandom3.h>
#include <TVector3.h>

//C, C++
#include <iostream>
#include <stdlib.h>
#include <assert.h>
#include <fstream>
#include <iomanip>
#include <time.h>
#include <bits/stdc++.h>
#include <sys/stat.h>

using namespace std;

Double_t anaTrgB::get_n_tot_ev_after_cuts(){
  Double_t n_tot_ev_after_cuts = 0.0;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry); nbytes += nb;
    //
    if(cut())
      n_tot_ev_after_cuts++;
  }
  return n_tot_ev_after_cuts;
}

void anaTrgB::Loop(TString histOut,
		   Float_t NGB_rate_in_MHz, Float_t fadc_electronic_noise_RMS,
		   Int_t npe_min, Int_t npe_max, Int_t nEvSim_max, Int_t nEv_max,
		   Int_t rndseed, bool NGBsim){
  //
  cout<<"anaTrgB::Loop             "<<endl
      <<"NGB_rate_in_MHz           "<<NGB_rate_in_MHz<<endl
      <<"fadc_electronic_noise_RMS "<<fadc_electronic_noise_RMS<<endl
      <<"npe_min                   "<<npe_min<<endl
      <<"npe_max                   "<<npe_max<<endl
      <<"nEvSim_max                "<<nEvSim_max<<endl
      <<"nEv_max                   "<<nEv_max<<endl
      <<"rndseed                   "<<rndseed<<endl
      <<"NGBsim                    "<<NGBsim<<endl;
  //
  _anaConf.readFromFile(_name_ana_conf_file);
  _anaConf.printInfo();
  //
  Double_t flux_pers_persr_perm2 = anaTrgB::get_proton_flux_pers_persr_perm2(_anaConf.eMin, _anaConf.eMax);
  Double_t solidAngle = anaTrgB::get_solid_angle(_anaConf.thetaSimDeg);
  Double_t generationArea = TMath::Pi()*_anaConf.rsim*_anaConf.rsim;
  Double_t tot_flux = flux_pers_persr_perm2*solidAngle*generationArea;
  //
  Double_t simulation_sampling_factor = anaTrgB::calculate_simulation_sampling_factor(_anaConf.eMin, _anaConf.eMax);
  Double_t tot_number_of_sim_events = _anaConf.nEvPerFileSim*_anaConf.nShowerMultiplicationFactor*simulation_sampling_factor;
  cout<<"tot_flux                 "<<tot_flux<<endl
      <<"tot_number_of_sim_events "<<tot_number_of_sim_events<<endl;
  //
  _NGB_rate_in_MHz = NGB_rate_in_MHz;
  _fadc_electronic_noise_RMS = fadc_electronic_noise_RMS;
  //
  _npe_min = npe_min;
  _npe_max = npe_max;
  _nEv_max = nEv_max;
  _nEvSim_max = nEvSim_max;
  //
  print_cuts();
  //
  _n_tot_ev_after_cuts = get_n_tot_ev_after_cuts();
  //
  //
  Int_t nevsim = 0;
  Int_t nevcut = 0;
  Double_t x0_LST01 = -70.93;
  Double_t y0_LST01 = -52.07;  
  //
  clock_t start, finish;
  clock_t start_trg, finish_trg;
  clock_t start_sim, finish_sim;
  start = clock();
  //
  TH1D *h1_digital_sum = new TH1D("h1_digital_sum", "h1_digital_sum", 20001,-0.5,20000.5);
  TH1D *h1_digital_sum_rate = new TH1D("h1_digital_sum_rate", "h1_digital_sum_rate", 20001,-0.5,20000.5);
  TH1D *h1_fadc_val = new TH1D("h1_fadc_val","h1_fadc_val", (Int_t)(400.5-249.5), 249.5, 400.5);
  //
  TH1D *h1_energy = new TH1D("h1_energy","h1_energy", 100000, 0.0, 100000.0); 
  TH2D *h2_core_pos = new TH2D("h2_core_pos","h2_core_pos", 400, -2000.0, 2000.0, 400, -2000.0, 2000.0);
  TH1D *h1_npe = new TH1D("h1_npe","h1_npe", 10001, -0.5, 10000.5);
  TH1D *h1_npe_w = new TH1D("h1_npe_w","h1_npe_w", 10001, -0.5, 10000.5);
  TH1D *h1_rate_vs_npe = new TH1D("h1_rate_vs_npe","h1_rate_vs_npe", 10001, -0.5, 10000.5);
  TH1D *h1_ev_time = new TH1D("h1_ev_time","h1_ev_time", 10000, -5000.0, 5000.0);
  //  
  TH1D *h1_dbscan_N_points = new TH1D("h1_dbscan_N_points","h1_dbscan_N_points", 10001, -0.5, 10000.5);
  TH1D *h1_n_micro_clusters = new TH1D("h1_n_micro_clusters","h1_n_micro_clusters", 10001, -0.5, 10000.5);
  TH1D *h1_N_dbc = new TH1D("h1_N_dbc","h1_N_dbc", 21, -0.5, 20.5);
  TH1D *h1_dbc_number_of_points = new TH1D("h1_dbc_number_of_points","h1_dbc_number_of_points", 1001, -0.5, 1000.5);
  TH1D *h1_dbc_number_of_points_w = new TH1D("h1_dbc_number_of_points_w","h1_dbc_number_of_points_w", 1001, -0.5, 1000.5);
  TH1D *h1_rate_vs_dbc_number_of_points = new TH1D("h1_rate_vs_dbc_number_of_points","h1_rate_vs_dbc_number_of_points", 1001, -0.5, 1000.5);
  //  
  TH1D *h1_dbc_nsb_rate = new TH1D("h1_dbc_nsb_rate","h1_dbc_nsb_rate", 1001, -0.5, 1000.5);
  //
  TH1D *h1_dbc_mean_x = new TH1D("h1_dbc_mean_x","h1_dbc_mean_x", 300, -2.0, 2.0);
  TH1D *h1_dbc_mean_y = new TH1D("h1_dbc_mean_y","h1_dbc_mean_y", 300, -2.0, 2.0);
  TH1D *h1_dbc_mean_time_ii = new TH1D("h1_dbc_mean_time_ii","h1_dbc_mean_time_ii", 100, 0.0, 100.0);
  //
  const unsigned int nn_PMT_channels = nChannels;
  Int_t fadc_sum_offset = 15;
  Int_t fadc_MHz = 1024;
  Int_t fadc_offset = 300;
  Float_t fadc_sample_in_ns = 1000.0/fadc_MHz;
  Float_t time_offset = fadc_sum_offset*fadc_sample_in_ns;
  //
  TRandom3 *rnd = new TRandom3(rndseed);
  //
  vector<vector<Int_t>> wfcam(nn_PMT_channels, vector<Int_t>(nn_fadc_point));
  wfCamSim *wfc = new wfCamSim(rnd, "Template_CTA_SiPM.txt", "spe.dat",
			       nn_fadc_point, nn_PMT_channels, fadc_offset, fadc_sample_in_ns, _NGB_rate_in_MHz, _fadc_electronic_noise_RMS);
  wfc->print_wfCamSim_configure();  
  //
  sipmCameraHist *sipm_cam = new sipmCameraHist("sipm_cam","sipm_cam","pixel_mapping.csv",0);
  triggerSim *trg_sim = new triggerSim(sipm_cam);    
  cout<<"trg_sim->_trg_setup.load_trg_setup"<<endl
      <<"_trg_conf_file = "<<_trg_conf_file<<endl;
  trg_sim->_trg_setup.load_trg_setup(_trg_conf_file.Data());
  trg_sim->fill_trg_channel_mask_from_file(trg_sim->_trg_setup.trigger_channel_mask_file_list);
  //trg_sim->set_k_dist_graph_flag(true);
  //
  //cout<<"trg_sim->set_k_dist_graph_flag(false)"<<endl;
  //trg_sim->set_k_dist_graph_flag(false);
  //
  finish = clock();
  cout<<"initialization time : "<<((finish - start)/CLOCKS_PER_SEC)<<" (sec)"<<endl;	
  Long64_t nentries = fChain->GetEntriesFast();
  cout<<"nentries = "<<nentries<<endl;
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry); nbytes += nb;
    //
    if((nevsim >= _nEvSim_max) && (_nEvSim_max > -1))
      break;
    if(nevcut >= _nEv_max && (_nEv_max > -1))
      break;
    if(nevcut >= _anaConf.n_ev_cuts_max && (_anaConf.n_ev_cuts_max > -1))
      break;
    //
    if(cut()){
      if(NGBsim){
	h1_npe->Fill(0);
	h1_npe_w->Fill(0);
	h1_energy->Fill(0);
	h2_core_pos->Fill(0.0,0.0);
	h1_ev_time->Fill(0);
      }
      else{
	h1_npe->Fill(n_pe);
	h1_npe_w->Fill(n_pe,evstHist::get_Weight_ETeV(energy));
	h1_energy->Fill(energy*1000);
	h2_core_pos->Fill(xcore,ycore);
	h1_ev_time->Fill(ev_time);
      }
      //
      start_sim = clock();
      if(_anaConf.wf_trg_sim){
	if(NGBsim)
	  wfc->simulate_cam_event_NGB(wfcam);
	else
	  wfc->simulate_cam_event(nn_fadc_point,
				  nn_PMT_channels,
				  wfcam,
				  ev_time,
				  time_offset,
				  n_pe,
				  pe_chID,
				  pe_time);
	finish_sim = clock();
	start_trg = clock();
	//
	if(NGBsim)
	  cout<<"---------------------"<<endl<<"n_pe = "<<0<<endl;
	else
	  cout<<"---------------------"<<endl<<"n_pe = "<<n_pe<<endl;
	trg_sim->get_trigger(wfcam,h1_digital_sum,h1_fadc_val);
	finish_trg = clock();
	//
	cout<<nevcut
	    <<" "<<((finish_sim - start_sim)/(CLOCKS_PER_SEC/1000))<<" (msec)"
	    <<" "<<((finish_trg - start_trg)/(CLOCKS_PER_SEC/1000))<<" (msec)"<<endl;
      }
      //
      //
      h1_dbscan_N_points->Fill(trg_sim->get_dbscan_N_points());
      h1_n_micro_clusters->Fill(trg_sim->get_n_digital_sum_micro_clusters());      
      h1_N_dbc->Fill(trg_sim->get_dbclusters().size());
      //cout<<"trg_sim->get_dbscan_N_points()              "<<trg_sim->get_dbscan_N_points()<<endl
      //  <<"trg_sim->get_n_digital_sum_micro_clusters() "<<trg_sim->get_n_digital_sum_micro_clusters()<<endl      
      //  <<"trg_sim->get_dbclusters().size()            "<<trg_sim->get_dbclusters().size()<<endl;
      //      
      if(trg_sim->get_dbclusters().size() == 1){
	//
	h1_dbc_number_of_points->Fill(trg_sim->get_dbclusters().at(0).number_of_points);
	h1_dbc_number_of_points_w->Fill(trg_sim->get_dbclusters().at(0).number_of_points, evstHist::get_Weight_ETeV(energy));
	//
	h1_dbc_mean_x->Fill(trg_sim->get_dbclusters().at(0).mean_x);
	h1_dbc_mean_y->Fill(trg_sim->get_dbclusters().at(0).mean_y);
	h1_dbc_mean_time_ii->Fill(trg_sim->get_dbclusters().at(0).mean_time_ii);
	//
      }
      else if (trg_sim->get_dbclusters().size() > 1) {
	Int_t number_of_points_max = 0;
	unsigned int k_max;
	for(unsigned int k = 0; k<trg_sim->get_dbclusters().size(); k++){
	  if(number_of_points_max < trg_sim->get_dbclusters().at(k).number_of_points){
	    number_of_points_max = trg_sim->get_dbclusters().at(k).number_of_points;
	    k_max = k;
	  }
	}
	//
	h1_dbc_number_of_points->Fill(trg_sim->get_dbclusters().at(k_max).number_of_points);
	h1_dbc_number_of_points_w->Fill(trg_sim->get_dbclusters().at(k_max).number_of_points, evstHist::get_Weight_ETeV(energy));
	//
	h1_dbc_mean_x->Fill(trg_sim->get_dbclusters().at(k_max).mean_x);
	h1_dbc_mean_y->Fill(trg_sim->get_dbclusters().at(k_max).mean_y);
	h1_dbc_mean_time_ii->Fill(trg_sim->get_dbclusters().at(k_max).mean_time_ii);
	//
      }
      //
      _n_tot_ev_after_cuts_processed++;
      nevcut++;
    }
    nevsim++;
  }
  //
  //
  cout<<"nentries                       = "<<nentries<<endl
      <<"_n_tot_ev_after_cuts           = "<<_n_tot_ev_after_cuts<<endl
      <<"_n_tot_ev_after_cuts_processed = "<<_n_tot_ev_after_cuts_processed<<endl;
  //
  //
  //
  // rate culculations
  anaTrgB::digital_sum_rate_calculations( h1_digital_sum_rate, h1_digital_sum, h1_npe, nn_fadc_point, fadc_sample_in_ns);
  anaTrgB::dbc_nsb_rate_calculations( h1_dbc_nsb_rate, h1_dbc_number_of_points, h1_npe, nn_fadc_point, fadc_sample_in_ns);
  anaTrgB::rate_calculations( h1_rate_vs_npe, h1_npe_w, tot_flux, tot_number_of_sim_events, _n_tot_ev_after_cuts, _n_tot_ev_after_cuts_processed);
  anaTrgB::rate_calculations( h1_rate_vs_dbc_number_of_points, h1_dbc_number_of_points_w, tot_flux, tot_number_of_sim_events, _n_tot_ev_after_cuts, _n_tot_ev_after_cuts_processed);
  //
  //
  TFile* rootFile = new TFile(histOut.Data(), "RECREATE", " Histograms", 1);
  rootFile->cd();
  if (rootFile->IsZombie()){
    cout<<"  ERROR ---> file "<<histOut.Data()<<" is zombi"<<endl;
    assert(0);
  }
  else
    cout<<"  Output Histos file ---> "<<histOut.Data()<<endl;
  //
  h1_npe->Write();
  h1_npe_w->Write();
  h1_rate_vs_npe->Write();
  h1_energy->Write();
  h2_core_pos->Write();
  h1_ev_time->Write();
  //
  h1_dbscan_N_points->Write();
  h1_n_micro_clusters->Write();
  //
  h1_digital_sum->Write();
  h1_digital_sum_rate->Write();
  h1_fadc_val->Write();
  //
  h1_N_dbc->Write();
  //
  h1_dbc_number_of_points->Write();
  h1_dbc_number_of_points_w->Write();
  h1_rate_vs_dbc_number_of_points->Write();
  //
  h1_dbc_nsb_rate->Write();
  //
  h1_dbc_mean_x->Write();
  h1_dbc_mean_y->Write();
  h1_dbc_mean_time_ii->Write();
  //
  rootFile->Close();
}

Bool_t anaTrgB::cut(){
  if(n_pe>=_npe_min && n_pe<=_npe_max){
    return true;
  }
  return false;
}

const void anaTrgB::print_cuts(){
  cout<<" _npe_min    "<<_npe_min<<endl
      <<" _npe_max    "<<_npe_max<<endl
      <<" _nEv_max    "<<_nEv_max<<endl
      <<" _nEvSim_max "<<_nEvSim_max<<endl;
}

void anaTrgB::digital_sum_rate_calculations( TH1D *h1_digital_sum_rate, TH1D *h1_digital_sum, TH1D *h1_npe, Int_t nn_fadc_point, Float_t fadc_sample_in_ns){
  Double_t integral_val;
  Double_t rate;
  Double_t time_in_s = nn_fadc_point*fadc_sample_in_ns*1.0e-9*h1_npe->GetEntries();
  //cout<<"time_in_s                                               "<<time_in_s<<endl
  //    <<"h1_digital_sum->Integral(1,h1_digital_sum->GetNbinsX()) "<<h1_digital_sum->Integral(1,h1_digital_sum->GetNbinsX())<<endl;
  for( Int_t i = 1; i <= h1_digital_sum->GetNbinsX(); i++){
    integral_val = h1_digital_sum->Integral(i,h1_digital_sum->GetNbinsX());
    if(time_in_s>0.0){
      rate = integral_val/time_in_s;
      h1_digital_sum_rate->SetBinContent(i,rate);    
    }
  }  
}   

void anaTrgB::dbc_nsb_rate_calculations( TH1D *h1_dbc_nsb_rate, TH1D *h1_dbc_number_of_points, TH1D *h1_npe, Int_t nn_fadc_point, Float_t fadc_sample_in_ns){
  Double_t integral_val;
  Double_t rate;
  Double_t time_in_s = nn_fadc_point*fadc_sample_in_ns*1.0e-9*h1_npe->GetEntries();
  for( Int_t i = 1; i <= h1_dbc_number_of_points->GetNbinsX(); i++){
    integral_val = h1_dbc_number_of_points->Integral(i,h1_dbc_number_of_points->GetNbinsX());
    if(time_in_s>0.0){
      rate = integral_val/time_in_s;
      h1_dbc_nsb_rate->SetBinContent(i,rate);
    }
  }
}

void anaTrgB::rate_calculations( TH1D *h1_rate_vs_th, TH1D *h1_w, Double_t tot_flux, Double_t tot_number_of_sim_events,
				 Double_t n_tot_ev_after_cuts, Double_t n_tot_ev_after_cuts_processed){
  Double_t integral_val;
  Double_t rate;
  Double_t tot_number_of_sim_events_correction_factor = 1.0;
  if(n_tot_ev_after_cuts>0)
    tot_number_of_sim_events_correction_factor = n_tot_ev_after_cuts_processed/n_tot_ev_after_cuts;
  for( Int_t i = 1; i <= h1_w->GetNbinsX(); i++){
    integral_val = h1_w->Integral(i,h1_w->GetNbinsX());
    if(tot_number_of_sim_events>0 && tot_number_of_sim_events_correction_factor>0){
      rate = integral_val/(tot_number_of_sim_events*tot_number_of_sim_events_correction_factor)*tot_flux;
    }
    else{
      rate = 0.0;
    }
    h1_rate_vs_th->SetBinContent(i,rate);
  }
}

Double_t anaTrgB::get_proton_flux_pers_persr_perm2(Double_t e_min_GeV, Double_t e_max_GeV){
  return 5781.68*(1.0/TMath::Power(e_min_GeV,1.674)-1.0/TMath::Power(e_max_GeV,1.674));
}

Double_t anaTrgB::get_solid_angle(Double_t theta_deg){
  return 2.0*TMath::Pi()*(1.0 - TMath::Cos(theta_deg/180.0*TMath::Pi()));
}

Double_t anaTrgB::calculate_simulation_sampling_factor(Double_t e_min_GeV, Double_t e_max_GeV){
  //
  TH1D *h1 = new TH1D("h1","h1",10000, e_min_GeV, e_max_GeV);
  TH1D *h1w = new TH1D("h1w","h1w",10000, e_min_GeV, e_max_GeV);
  Double_t esim;
  //
  TRandom3 *rnd = new TRandom3(123123);
  for(Int_t i = 0;i<100000000;i++){
    esim = 1.0/(rnd->Uniform(1.0/e_max_GeV,1.0/e_min_GeV));
    h1->Fill(esim);
    h1w->Fill(esim,evstHist::get_Weight_ETeV(esim/1000.0));
  }
  Double_t integral_val = h1->Integral(1,h1->GetNbinsX());
  Double_t integralw_val = h1w->Integral(1,h1w->GetNbinsX());
  cout<<"integral_val  "<<integral_val<<endl
      <<"integralw_val "<<integralw_val<<endl;
  if(integral_val>0)
    cout<<"factor      "<<integralw_val/integral_val<<endl;
  //
  if(integral_val>0)
    return integralw_val/integral_val;
  return 0.0;
}
