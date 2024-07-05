//my
#include "anaTrg.hh"
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

void anaTrg::Loop(TString histOut){
  //
  _v_det.SetXYZ(1.0*TMath::Sin(20.0/180.0*TMath::Pi()),0,1.0*TMath::Cos(20.0/180.0*TMath::Pi()));
  //
  clock_t start, finish;
  clock_t start_trg, finish_trg;
  clock_t start_sim, finish_sim;
  start = clock();
  //
  TH1D *h1_n_pe = new TH1D("h1_n_pe","h1_n_pe",1000,0.0,1000000);
  TH1D *h1_n_pe_zoom = new TH1D("h1_n_pe_zoom","h1_n_pe_zoom",1001,-0.5,1000.5);
  TH1D *h1_n_pe_9bin = new TH1D("h1_n_pe_9bin","h1_n_pe_9bin",1000, 1.0,9001);
  TH1D *h1_n_pe_bins = new TH1D();
  h1_n_pe_bins->SetNameTitle("h1_n_pe_bins","h1_n_pe_bins");
  set_n_pe_bins(h1_n_pe_bins,_npe_max,_npe_min);
  dump_n_pe_bins(h1_n_pe_bins);
  //
  TH1D *h1_digital_sum     = new TH1D("h1_digital_sum","h1_digital_sum",1001,-0.5,1000.5);
  TH1D *h1_digital_sum_3ns = new TH1D("h1_digital_sum_3ns","h1_digital_sum_3ns",1001,-0.5,1000.5);
  TH1D *h1_digital_sum_5ns = new TH1D("h1_digital_sum_5ns","h1_digital_sum_5ns",1001,-0.5,1000.5);
  TH1D *h1_fadc_val      = new TH1D("h1_fadc_val","h1_fadc_val",1001,-0.5,1000.5);
  //TH1D *h1_fadc_val = new TH1D("h1_fadc_val","h1_fadc_val", 249.5, 400.5, (Int_t)(400.5-249.5));
  //
  //
  TH1D *h1_n_pe_vs_chID    = new TH1D("h1_n_pe_vs_chID","h1_n_pe_vs_chID",nChannels+1,-0.5,nChannels+0.5);
  //
  cout<<"_npe_max = "<<_npe_max<<endl
      <<"_npe_min = "<<_npe_min<<endl;
  //
  /////////////
  //
  //
  Double_t x0_LST01 = -70.93;
  Double_t y0_LST01 = -52.07;
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
  TRandom3 *rnd = new TRandom3(123123);
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
  evstHist *evH = new evstHist("evH","evH");
  evH->test();
  evstHist *evH02 = new evstHist("evH02","evH02");
  evH02->test_get_bin(1.1, 0.1, 1);
  evH02->test_get_bin(10.1, 0.1, 2);
  evH02->test_get_bin(10.1, 9.1, 3);
  evH02->test_get_bin(100.1, 9.1, 4);
  evH02->test_get_bin(1000.1, 9.1, 5);
  evH02->test_get_bin(10000.1, 8.1, 6);
  evH02->test_get_bin(99999.1, 7.1, 7);
  //
  TGraph *gr_WF_tmpl_array = new TGraph();
  gr_WF_tmpl_array->SetNameTitle("gr_WF_tmpl_array","gr_WF_tmpl_array");
  wfc->get_gr_WF_tmpl_array(gr_WF_tmpl_array);
  //
  sipmCameraHist *sipm_cam = new sipmCameraHist("sipm_cam","sipm_cam","pixel_mapping.csv",0);
  triggerSim *trg_sim = new triggerSim(sipm_cam);    
  //
  vector<TGraph*> v_gr;
  for(Int_t j = 0;j<nChannels;j++){
    TString gr_name_title = "gr_ch_";
    gr_name_title += j;
    TGraph *gr = new TGraph();
    gr->SetNameTitle(gr_name_title.Data(),gr_name_title.Data());
    v_gr.push_back(gr);
  }
  //
  Int_t nevsim = 0;
  //
  finish = clock();
  cout<<"initialization time : "<<((finish - start)/CLOCKS_PER_SEC)<<" (sec)"<<endl;	
  Long64_t nentries = fChain->GetEntriesFast();
  cout<<"nentries = "<<nentries<<endl;
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    if(jentry%1000000 == 0)
      cout<<jentry<<endl;
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    //if (ientry > 10) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    //
    Double_t rcore = TMath::Sqrt((x0_LST01 - xcore)*(x0_LST01 - xcore) + (y0_LST01 - ycore)*(y0_LST01 - ycore));
    //if(cut(h1_n_pe_bins)){
    Double_t theta_p_t = get_theta_p_t();
    Double_t theta_p_t_deg = theta_p_t*180/TMath::Pi();
    if(cut()){
      if(nevsim<100 && n_pe < 49 ){
	if( rcore<150.0){
	  if(theta_p_t_deg>1.0 && theta_p_t_deg<=2.0){
	  //if(nevsim<100){
	  cout<<"jentry = "<<jentry<<endl;
	  h1_n_pe->Fill(n_pe);
	  h1_n_pe_zoom->Fill(n_pe);
	  h1_n_pe_9bin->Fill(n_pe);
	  h1_n_pe_bins->Fill(n_pe);
	  start_sim = clock();
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
	  TString trg_vector_out_file = "ev_synthetic_trg_v_";
	  trg_vector_out_file += (Int_t)jentry;
	  trg_vector_out_file += "ev.csv";
	  //triggerSim::print_trigger_vec_to_csv(trg_sim->get_trigger(wfcam,h1_digital_sum,h1_digital_sum_3ns,h1_digital_sum_5ns,h1_fadc_val),
	  //				     sipm_cam,
	  //				     trg_vector_out_file);
	  //
	  //std::vector<std::vector<unsigned int>> trg_vec = trg_sim->get_trigger(wfcam);
	  //triggerSim::print_trigger_vec(trg_vec);
	  //trg_sim->get_trigger(wfcam,h1_digital_sum,h1_digital_sum_3ns,h1_digital_sum_5ns,h1_fadc_val);
	  trg_sim->get_trigger(wfcam,h1_digital_sum,h1_fadc_val);
	  finish_trg = clock();
	  cout<<nevsim
	      <<" "<<((finish_sim - start_sim)/(CLOCKS_PER_SEC/1000))<<" (msec)"
	      <<" "<<((finish_trg - start_trg)/(CLOCKS_PER_SEC/1000))<<" (msec)"<<endl;
	  //
	  //
	  //for(Int_t i = 0;i<n_pe;i++)
	  //h1_n_pe_vs_chID->Fill(pe_chID[i]);
	  //
	  //
	  //for(Int_t i = 0;i<nn_PMT_channels;i++){
	  //for(Int_t j = 0;j<nn_fadc_point;j++){
	  // v_gr.at(i)->SetPoint(j,j*fadc_sample_in_ns,wfcam[i][j]);
	  //}
	  //}
	  nevsim++;
	  }
	}
      }
    }
    //
  }
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
  for(Int_t i = 0;i<nn_PMT_channels;i++)
    v_gr.at(i)->Write();
  //
  h1_n_pe->Write();
  h1_n_pe_zoom->Write();
  h1_n_pe_9bin->Write();
  h1_n_pe_bins->Write();
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
  h1_n_pe_vs_chID->Write();
  evH->Draw_hist("")->Write();
  evH02->Draw_hist("")->Write();
  //
  rootFile->Close();
}

void anaTrg::set_n_pe_bins(TH1D *h1, Int_t &npe_max, Int_t &npe_min){
  npe_min = 1;
  npe_max = 500;
  const Int_t nn = 25;
  Double_t* xBins = new Double_t[nn];
  xBins[0] = npe_min;
  xBins[1] = 2;
  for(Int_t i = 0;i<20;i++)
    xBins[1+i+1] = (i+1)*10;
  //
  xBins[22] = 300;
  xBins[23] = 400;
  xBins[24] = npe_max;
  //
  h1->SetBins(nn-1,xBins);
}

void anaTrg::dump_n_pe_bins(TH1D *h1){
  for(Int_t i = 0;i<h1->GetNbinsX();i++)
    cout<<setw(10)<<h1->GetBinLowEdge(i+1)<<setw(10)<<(h1->GetBinLowEdge(i+1)+h1->GetBinWidth(i+1))<<endl;
}

Bool_t anaTrg::cut(TH1D *h1){
  //if(n_pe<=_npe_max && n_pe>=_npe_min){
  //if(h1 != NULL){
  //if(h1->GetBinContent(h1->FindBin(n_pe))<1000)
  //return true;
  //}
  //return false;
  //}
  //
  //if(TMath::Abs(xcore)<150)
  //if(TMath::Abs(ycore)<150)
  if(n_pe==1)
    return true;
  return false;
}

Double_t anaTrg::get_theta_p_t(){
  TVector3 v_prot;
  v_prot.SetMagThetaPhi(1.0,TMath::Pi()/2.0-altitude,TMath::Pi() - azimuth);
  TVector3 v_prot_inv(v_prot.x(),v_prot.y(),v_prot.z());
  return TMath::ACos(v_prot_inv.Dot(_v_det)/v_prot_inv.Mag()/_v_det.Mag());
}
