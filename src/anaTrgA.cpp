//my
#include "anaTrgA.hh"
#include "ana.hh"
#include "anabase.hh"
#include "sipmCameraHist.hh"
#include "wfCamSim.hh"
#include "triggerSim.hh"
#include "evstHist.hh"
#include "dbscan.hh"

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

Double_t anaTrgA::get_theta_p_t( const TVector3 &v_det, Double_t altitude_v, Double_t azimuth_v){
  TVector3 v_prot;
  v_prot.SetMagThetaPhi(1.0,TMath::Pi()/2.0-altitude_v,TMath::Pi() - azimuth_v);
  TVector3 v_prot_inv(v_prot.x(),v_prot.y(),v_prot.z());
  return TMath::ACos(v_prot_inv.Dot(v_det)/v_prot_inv.Mag()/v_det.Mag());
}

void anaTrgA::transform_SiPM_distadd( TGraph *gr_tr, const TGraph *gr){
  Double_t x,y;
  Double_t x_new,y_new; 
  Double_t ymax = 0.0;
  for(Int_t i = 0;i<gr->GetN();i++){
    gr->GetPoint(i,x,y);
    if(y>ymax)
      ymax = y;
  }
  Double_t x_min = -5.0*1.0e+6;
  Double_t x_max = 30.0*1.0e+6;
  Int_t nn = 1000;
  for(Int_t i = 0;i<nn;i++){
    //gr->GetPoint(i,x,y);
    x_new = x_min + (x_max - x_min)/(nn-1)*i;
    x = x_new/(1.0/0.920148*1.0e+6*7.5);
    y = gr->Eval(x);
    //x_new=x*1.0/0.920148*1.0e+6*7.5;
    y_new=y/ymax+20.0*TMath::Gaus(x_new,0.0,1.0e+6*7.5*0.12);
    y_new *= 0.0016;
    gr_tr->SetPoint(i,x_new,y_new);
  }
}

void anaTrgA::SiPM_dist(TString histOut){
  //
  const unsigned int nn_PMT_channels = nChannels;
  Int_t fadc_sum_offset = 15;
  Int_t fadc_MHz = 1024;
  Int_t fadc_offset = 300;
  Float_t fadc_sample_in_ns = 1000.0/fadc_MHz;
  Float_t time_offset = fadc_sum_offset*fadc_sample_in_ns;
  Float_t NGB_rate_in_MHz = 0.0;
  Float_t fadc_electronic_noise_RMS = 0.01;
  //
  TRandom3 *rnd = new TRandom3(123123);
  //
  vector<vector<Int_t>> wfcam(nn_PMT_channels, vector<Int_t>(nn_fadc_point));
  wfCamSim *wfc = new wfCamSim(rnd, "Template_CTA_SiPM.txt", "spe.dat",
			       nn_fadc_point, nn_PMT_channels, fadc_offset, fadc_sample_in_ns, NGB_rate_in_MHz, fadc_electronic_noise_RMS);
  wfc->print_wfCamSim_configure();
  //
  TGraph *gr_aml_gain_pedestal = new TGraph();
  transform_SiPM_distadd( gr_aml_gain_pedestal, wfc->get_gr_wf_ampl());
  //
  //----------------------
  TFile* rootFile = new TFile(histOut.Data(), "RECREATE", " Histograms", 1);
  rootFile->cd();
  if (rootFile->IsZombie()){
    cout<<"  ERROR ---> file "<<histOut.Data()<<" is zombi"<<endl;
    assert(0);
  }
  else
    cout<<"  Output Histos file ---> "<<histOut.Data()<<endl;
  //
  wfc->getTemplate()->Write();
  wfc->get_gr_wf_ampl()->Write();
  wfc->get_h1_wf_ampl_ADC()->Write();
  wfc->get_h1_wf_ampl()->Write();
  wfc->get_h1_adc_NGB_pedestal()->Write();
  wfc->get_h1_dadc_NGB_pedestal()->Write();
  //
  gr_aml_gain_pedestal->SetNameTitle("gr_aml_gain_pedestal","gr_aml_gain_pedestal");
  gr_aml_gain_pedestal->Write();
  //
  rootFile->Close();
}

void anaTrgA::test_single_pe_amplitude_generator(TString histOut){
  //
  //----------------------
  //
  const unsigned int nn_PMT_channels = nChannels;
  Int_t fadc_sum_offset = 15;
  Int_t fadc_MHz = 1024;
  Int_t fadc_offset = 300;
  Float_t fadc_sample_in_ns = 1000.0/fadc_MHz;
  Float_t time_offset = fadc_sum_offset*fadc_sample_in_ns;
  Float_t NGB_rate_in_MHz = 0.0;
  Float_t fadc_electronic_noise_RMS = 0.01;
  //
  TRandom3 *rnd = new TRandom3(123123);
  //
  wfCamSim *wfc = new wfCamSim(rnd, "Template_CTA_SiPM.txt", "spe.dat",
			       nn_fadc_point, nn_PMT_channels, fadc_offset, fadc_sample_in_ns, NGB_rate_in_MHz, fadc_electronic_noise_RMS);
  wfc->print_wfCamSim_configure();
  //
  //----------------------
  //
  Int_t n_pe_to_sim = 10000000;
  TH1D *h1_single_pe_amplitude_generator = new TH1D("h1_single_pe_amplitude_generator","h1_single_pe_amplitude_generator",100, 0.0, 100.0);
  TH1D *h1_single_pe_amplitude_from_hist_generator = new TH1D("h1_single_pe_amplitude_from_hist_generator","h1_single_pe_amplitude_from_hist_generator",100, 0.0, 100.0);
  TH1D *h1_single_pe_amplitude_invf_generate = new TH1D("h1_single_pe_amplitude_invf_generate","h1_single_pe_amplitude_invf_generate",100, 0.0, 100.0);
  wfc->test_single_pe_amplitude_generator( h1_single_pe_amplitude_generator,  n_pe_to_sim);
  wfc->test_single_pe_amplitude_from_hist_generator( h1_single_pe_amplitude_from_hist_generator, n_pe_to_sim);
  wfc->test_single_pe_amplitude_invf_generate(h1_single_pe_amplitude_invf_generate, n_pe_to_sim);
  //
  //----------------------
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
  wfc->getTemplate()->Write();
  wfc->get_gr_wf_ampl()->Write();
  wfc->get_h1_wf_ampl_ADC()->Write();
  wfc->get_h1_wf_ampl_ADC_integral()->Write();  
  wfc->get_h1_wf_ampl()->Write();
  wfc->get_h1_adc_NGB_pedestal()->Write();
  wfc->get_h1_dadc_NGB_pedestal()->Write();
  //
  h1_single_pe_amplitude_generator->Write();
  h1_single_pe_amplitude_from_hist_generator->Write();
  h1_single_pe_amplitude_invf_generate->Write();
  //
  rootFile->Close();
}

void anaTrgA::Loop(TString histOut, Int_t binE, Int_t binTheta, Int_t binDist, Int_t npe_min, Int_t npe_max, Int_t nEv_max, Int_t rndseed, Int_t data_chunk_ID, bool NGBsim){
  //
  cout<<"anaTrgA::Loop"<<endl
      <<"data_chunk_ID                            "<<data_chunk_ID<<endl
      <<"NGBsim                                   "<<NGBsim<<endl
      <<"_n_data_chunks                           "<<_n_data_chunks<<endl
      <<"_disable_energy_theta_rcore_binwise_cuts "<<_disable_energy_theta_rcore_binwise_cuts<<endl
      <<"_rsimulation                             "<<_rsimulation<<endl;
  //
  TVector3 v_det;
  v_det.SetXYZ(1.0*TMath::Sin(20.0/180.0*TMath::Pi()),0,1.0*TMath::Cos(20.0/180.0*TMath::Pi()));
  //
  evstHist *evH = new evstHist("evH","evH");
  evH->Print_hist_BinsInfo();
  evH->Print_hist_BinsInfo(binE, binTheta, binDist);
  //
  evH->get_Bin_Edge(evH->get_E_hist(),binE,_E_min, _E_max);
  evH->get_Bin_Edge(evH->get_theta_hist(),binTheta, _theta_deg_min, _theta_deg_max);
  evH->get_Bin_Edge(evH->get_v_r().at(0),binDist, _dist_min, _dist_max);
  //    
  _npe_min = npe_min;
  _npe_max = npe_max;
  _nEv_max = nEv_max;
  //
  print_cuts();
  //
  clock_t start, finish;
  clock_t start_trg, finish_trg;
  clock_t start_sim, finish_sim;
  start = clock();
  //
  TH1D *h1_digital_sum     = new TH1D("h1_digital_sum",     "h1_digital_sum",    1001,-0.5,1000.5);
  TH1D *h1_digital_sum_3ns = new TH1D("h1_digital_sum_3ns", "h1_digital_sum_3ns",1001,-0.5,1000.5);
  TH1D *h1_digital_sum_5ns = new TH1D("h1_digital_sum_5ns", "h1_digital_sum_5ns",1001,-0.5,1000.5);
  TH1D *h1_fadc_val        = new TH1D("h1_fadc_val",        "h1_fadc_val",       1001,-0.5,1000.5);
  //
  evstHist *evH_h1_energy_all = new evstHist("evH_h1_energy_all","evH_h1_energy_all");
  evstHist *evH_h1_energy_trg = new evstHist("evH_h1_energy_trg","evH_h1_energy_trg");
  evstHist *evH_h1_energy_eff_r = new evstHist("evH_h1_energy_eff_r","evH_h1_energy_eff_r");
  //
  evH_h1_energy_all->get_E_hist()->SetNameTitle("_h1_energy_all","_h1_energy_all");
  evH_h1_energy_trg->get_E_hist()->SetNameTitle("_h1_energy_trg","_h1_energy_trg");
  evH_h1_energy_eff_r->get_E_hist()->SetNameTitle("_h1_energy_eff_r","_h1_energy_eff_r");
  //
  TH1D *h1_energy = new TH1D("h1_energy","h1_energy", 100000, 0.0, 100000.0); 
  TH1D *h1_rcore = new TH1D("h1_rcore","h1_rcore", 1000, 0.0, 2000.0);
  TH1D *h1_theta_p_t_deg = new TH1D("h1_theta_p_t_deg","h1_theta_p_t_deg", 1000, 0.0, 10.0);
  TH1D *h1_npe = new TH1D("h1_npe","h1_npe", 1000, 0.0, 10000.0);
  TH1D *h1_n_micro_clusters = new TH1D("h1_n_micro_clusters","h1_n_micro_clusters", 1000, 0.0, 10000.0);
  //
  TH1D *h1_N_dbc = new TH1D("h1_N_dbc","h1_N_dbc", 21, -0.5, 20.5);
  TH1D *h1_dbc_number_of_points = new TH1D("h1_dbc_number_of_points","h1_dbc_number_of_points", 1001, -0.5, 1000.5);
  TH1D *h1_dbc_number_of_points_w = new TH1D("h1_dbc_number_of_points_w","h1_dbc_number_of_points_w", 1001, -0.5, 1000.5);
  TH1D *h1_dbc_number_of_points_tot = new TH1D("h1_dbc_number_of_points_tot","h1_dbc_number_of_points_tot", 1001, -0.5, 1000.5);
  TH1D *h1_dbc_number_of_points_tot_w = new TH1D("h1_dbc_number_of_points_tot_w","h1_dbc_number_of_points_tot_w", 1001, -0.5, 1000.5);
  TH1D *h1_dbc_number_of_CORE_POINT = new TH1D("h1_dbc_number_of_CORE_POINT","h1_dbc_number_of_CORE_POINT", 1001, -0.5, 1000.5);
  TH1D *h1_dbc_number_of_BORDER_POINT = new TH1D("h1_dbc_number_of_BORDER_POINT","h1_dbc_number_of_BORDER_POINT", 1001, -0.5, 1000.5);
  //
  TH2D *h2_dbc_number_of_points_vs_npe = new TH2D("h2_dbc_number_of_points_vs_npe","h2_dbc_number_of_points_vs_npe", 101, -0.5, 100.5, 201, -0.5, 200.5);
  TH2D *h2_dbc_number_of_points_vs_npe_w = new TH2D("h2_dbc_number_of_points_vs_npe_w","h2_dbc_number_of_points_vs_npe_w", 101, -0.5, 100.5, 201, -0.5, 200.5);
  TH1D *h1_dbc_number_of_points_npe_norm = new TH1D("h1_dbc_number_of_points_npe_norm","h1_dbc_number_of_points_npe_norm", 201, -0.5, 200.5);
  //
  TH1D *h1_dbc_mean_x = new TH1D("h1_dbc_mean_x","h1_dbc_mean_x", 300, -2.0, 2.0);
  TH1D *h1_dbc_mean_y = new TH1D("h1_dbc_mean_y","h1_dbc_mean_y", 300, -2.0, 2.0);
  TH1D *h1_dbc_mean_time_ii = new TH1D("h1_dbc_mean_time_ii","h1_dbc_mean_time_ii", 100, 0.0, 100.0);
  //
  TH2D *h2_dbc_number_of_points_vs_mean_time_ii = new TH2D("h2_dbc_number_of_points_vs_mean_time_ii","h2_dbc_number_of_points_vs_mean_time_ii", 101, -0.5, 100.5, 101, -0.5, 100.5);
  //
  TH2D *h2_fadc_val_vs_time_ii = new TH2D("h2_fadc_val_vs_time_ii","h2_fadc_val_vs_time_ii", 101, -0.5, 100.5, 201, 199.5, 400.5);  
  //
  TGraph *gr_event_id_vs_jentry = new TGraph();
  gr_event_id_vs_jentry->SetNameTitle("gr_event_id_vs_jentry","gr_event_id_vs_jentry");
  //
  TGraph *gr_dbscan_time_vs_npe = new TGraph();
  gr_dbscan_time_vs_npe->SetNameTitle("gr_dbscan_time_vs_npe","gr_dbscan_time_vs_npe");  
  //
  TGraph *gr_dbscan_time_vs_n_micro_clusters = new TGraph();
  gr_dbscan_time_vs_n_micro_clusters->SetNameTitle("gr_dbscan_time_vs_n_micro_clusters","gr_dbscan_time_vs_n_micro_clusters");  
  //
  Int_t number_of_points_tot = 0;  
  //
  /////////////
  //
  Double_t x0_LST01 = -70.93;
  Double_t y0_LST01 = -52.07;
  //
  Int_t nevsim = 0;
  //
  const unsigned int nn_PMT_channels = nChannels;
  Int_t fadc_sum_offset = 15;
  Int_t fadc_MHz = 1024;
  Int_t fadc_offset = 300;
  Float_t fadc_sample_in_ns = 1000.0/fadc_MHz;
  Float_t time_offset = fadc_sum_offset*fadc_sample_in_ns;
  //Float_t NGB_rate_in_MHz = 386.0;
  Float_t NGB_rate_in_MHz = 268.0;
  //Float_t NGB_rate_in_MHz = 0.0;
  //Float_t fadc_electronic_noise_RMS = 3.94;
  //Float_t fadc_electronic_noise_RMS = 0.01;
  Float_t fadc_electronic_noise_RMS = 3.8436441; //takes into account 3.0/sqrt(12)
  //
  bool if_dbc_trg = false;
  bool if_digital_sum_trg = true;
  //
  TRandom3 *rnd = new TRandom3(rndseed);
  //
  bool iftrg = false;
  //
  ////////////////////////////////////  
  //static const Int_t nChannels = 7987;
  //static const Int_t nn_fadc_point = 75;
  //
  vector<vector<Int_t>> wfcam(nn_PMT_channels, vector<Int_t>(nn_fadc_point));
  //
  wfCamSim *wfc = new wfCamSim(rnd, "Template_CTA_SiPM.txt", "spe.dat",
			       nn_fadc_point, nn_PMT_channels, fadc_offset, fadc_sample_in_ns, NGB_rate_in_MHz, fadc_electronic_noise_RMS);
  wfc->print_wfCamSim_configure();  
  //
  sipmCameraHist *sipm_cam = new sipmCameraHist("sipm_cam","sipm_cam","pixel_mapping.csv",0);
  triggerSim *trg_sim = new triggerSim(sipm_cam);    
  trg_sim->set_k_dist_graph_flag(_k_dist_graph_flag);
  cout<<"_k_dist_graph_flag = "<<_k_dist_graph_flag<<endl;
  //
  finish = clock();
  cout<<"initialization time : "<<((finish - start)/CLOCKS_PER_SEC)<<" (sec)"<<endl;	
  Long64_t nentries = fChain->GetEntriesFast();
  cout<<"nentries = "<<nentries<<endl;
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    //if(jentry%1000000 == 0)
    //cout<<jentry<<endl;
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry); nbytes += nb;
    //
    if( get_current_data_chunk_ID(nentries, jentry) == data_chunk_ID || data_chunk_ID == -999){
      if(nevsim >= _nEv_max)
	break;
      //
      Double_t rcore = TMath::Sqrt((x0_LST01 - xcore)*(x0_LST01 - xcore) + (y0_LST01 - ycore)*(y0_LST01 - ycore));
      //
      Double_t theta_p_t = get_theta_p_t( v_det, altitude, azimuth);
      Double_t theta_p_t_deg = theta_p_t*180/TMath::Pi();
      //
      //
      if(cut( nevsim, theta_p_t_deg, rcore)){
	//
	//
	gr_event_id_vs_jentry->SetPoint(gr_event_id_vs_jentry->GetN(),(Double_t)jentry,(Double_t)event_id);
	//
	//
	h1_theta_p_t_deg->Fill(theta_p_t_deg);
	h1_npe->Fill(n_pe);
	h1_energy->Fill(energy*1000);
	h1_rcore->Fill(rcore);
	//
	start_sim = clock();
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
	trg_sim->get_trigger(wfcam,h1_digital_sum,h1_digital_sum_3ns,h1_digital_sum_5ns,h1_fadc_val);
	finish_trg = clock();
	cout<<nevsim
	    <<" "<<((finish_sim - start_sim)/(CLOCKS_PER_SEC/1000))<<" (msec)"
	    <<" "<<((finish_trg - start_trg)/(CLOCKS_PER_SEC/1000))<<" (msec)"<<endl;
	//
	gr_dbscan_time_vs_npe->SetPoint(gr_dbscan_time_vs_npe->GetN(), n_pe, trg_sim->get_dbscan_run_time_musec());
	gr_dbscan_time_vs_n_micro_clusters->SetPoint(gr_dbscan_time_vs_n_micro_clusters->GetN(), trg_sim->get_dbscan_N_points(), trg_sim->get_dbscan_run_time_musec());
	//
	h1_n_micro_clusters->Fill(trg_sim->get_dbscan_N_points());
	//
	h1_N_dbc->Fill(trg_sim->get_dbclusters().size());
	number_of_points_tot = 0;
	if(if_dbc_trg){
	  iftrg = false;
	  for(unsigned int k = 0; k<trg_sim->get_dbclusters().size(); k++){
	    h1_dbc_number_of_points->Fill(trg_sim->get_dbclusters().at(k).number_of_points);
	    number_of_points_tot += trg_sim->get_dbclusters().at(k).number_of_points;
	    h1_dbc_number_of_points_w->Fill(trg_sim->get_dbclusters().at(k).number_of_points, evstHist::get_Weight_ETeV(energy));
	    h2_dbc_number_of_points_vs_npe->Fill(trg_sim->get_dbclusters().at(k).number_of_points,n_pe);
	    h2_dbc_number_of_points_vs_npe_w->Fill(trg_sim->get_dbclusters().at(k).number_of_points,n_pe, evstHist::get_Weight_ETeV(energy));
	    h1_dbc_number_of_points_npe_norm->Fill(n_pe);
	    h1_dbc_number_of_CORE_POINT->Fill(trg_sim->get_dbclusters().at(k).number_of_CORE_POINT);
	    h1_dbc_number_of_BORDER_POINT->Fill(trg_sim->get_dbclusters().at(k).number_of_BORDER_POINT);
	    //
	    h1_dbc_mean_x->Fill(trg_sim->get_dbclusters().at(k).mean_x);
	    h1_dbc_mean_y->Fill(trg_sim->get_dbclusters().at(k).mean_y);
	    h1_dbc_mean_time_ii->Fill(trg_sim->get_dbclusters().at(k).mean_time_ii);
	    //
	    h2_dbc_number_of_points_vs_mean_time_ii->Fill(trg_sim->get_dbclusters().at(k).mean_time_ii,
							  trg_sim->get_dbclusters().at(k).number_of_points);
	    if(trg_sim->get_dbclusters().at(k).number_of_points >= 39)
	      iftrg = true;
	  }
	  h1_dbc_number_of_points_tot->Fill(number_of_points_tot);
	  h1_dbc_number_of_points_tot_w->Fill(number_of_points_tot, evstHist::get_Weight_ETeV(energy));
	}
	else if(if_digital_sum_trg){
	  iftrg = false;
	  if(trg_sim->get_digital_sum_max()>312)
	    iftrg = true;
	  h1_dbc_number_of_points->Fill(trg_sim->get_digital_sum_max());
	  h1_dbc_number_of_points_w->Fill(trg_sim->get_digital_sum_max(), evstHist::get_Weight_ETeV(energy));
	  //
	  h1_dbc_number_of_points_tot->Fill(trg_sim->get_digital_sum_max());
	  h1_dbc_number_of_points_tot_w->Fill(trg_sim->get_digital_sum_max(), evstHist::get_Weight_ETeV(energy));
	}
	//
	//
	//
	evH_h1_energy_all->get_E_hist()->Fill(energy*1000);
	if(iftrg)
	  evH_h1_energy_trg->get_E_hist()->Fill(energy*1000);
	//
	//
	//trg_sim->fill_fadc_val_vs_time(wfcam,h2_fadc_val_vs_time_ii);
	//
	nevsim++;
      }
    }
  }
  //
  //----------------------
  evH_h1_energy_eff_r->Divideh1(evH_h1_energy_trg, evH_h1_energy_all, TMath::Pi()*_rsimulation*_rsimulation);
  //----------------------
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
  wfc->getTemplate()->Write();
  wfc->get_gr_wf_ampl()->Write();
  wfc->get_h1_wf_ampl_ADC()->Write();
  wfc->get_h1_wf_ampl()->Write();
  wfc->get_h1_adc_NGB_pedestal()->Write();
  wfc->get_h1_dadc_NGB_pedestal()->Write();
  //
  h1_digital_sum->Write();
  h1_digital_sum_3ns->Write();
  h1_digital_sum_5ns->Write();
  h1_fadc_val->Write();
  //
  h1_energy->Write();
  h1_rcore->Write();
  h1_theta_p_t_deg->Write();
  h1_npe->Write();
  //
  h1_n_micro_clusters->Write();
  h1_N_dbc->Write();
  h1_dbc_number_of_points->Write();
  h1_dbc_number_of_points_w->Write();
  h1_dbc_number_of_points_tot->Write();
  h1_dbc_number_of_points_tot_w->Write();
  h2_dbc_number_of_points_vs_npe->Write();
  h1_dbc_number_of_points_npe_norm->Write();
  h1_dbc_number_of_CORE_POINT->Write();
  h1_dbc_number_of_BORDER_POINT->Write();
  //
  h1_dbc_mean_x->Write();
  h1_dbc_mean_y->Write();
  h1_dbc_mean_time_ii->Write();
  //
  h2_dbc_number_of_points_vs_mean_time_ii->Write();
  gr_event_id_vs_jentry->Write();
  gr_dbscan_time_vs_npe->Write();
  gr_dbscan_time_vs_n_micro_clusters->Write();
  //
  evH_h1_energy_all->get_E_hist()->Write();
  evH_h1_energy_trg->get_E_hist()->Write();
  evH_h1_energy_eff_r->get_E_hist()->Write();
  //
  //h2_fadc_val_vs_time_ii->Write();
  //
  //cout<<"nevsim = "<<nevsim<<endl;
  //
  rootFile->Close();
}

Bool_t anaTrgA::cut( Int_t nevsim, Double_t theta_p_t_deg, Double_t rcore){
  if(nevsim<_nEv_max){
    if((energy>=_E_min/1000.0 && energy<_E_max/1000.0) || _disable_energy_theta_rcore_binwise_cuts){
      if((rcore>=_dist_min && rcore<_dist_max) || _disable_energy_theta_rcore_binwise_cuts){
	if((theta_p_t_deg>=_theta_deg_min && theta_p_t_deg<_theta_deg_max) || _disable_energy_theta_rcore_binwise_cuts){
	  if(n_pe>=_npe_min && n_pe<_npe_max){
	    return true;
	  }
	}
      }
    }
  }
  return false;
}

const void anaTrgA::print_cuts(){
  cout<<"_E_min         "<<_E_min<<endl
      <<"_E_max         "<<_E_max<<endl
      <<"_theta_deg_min "<<_theta_deg_min<<endl
      <<"_theta_deg_max "<<_theta_deg_max<<endl
      <<"_dist_min      "<<_dist_min<<endl
      <<"_dist_max      "<<_dist_max<<endl
      <<"_npe_min       "<<_npe_min<<endl
      <<"_npe_max       "<<_npe_max<<endl
      <<"_nEv_max       "<<_nEv_max<<endl;
}
