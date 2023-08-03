//my
#include "anaPCA.hh"
#include "sipmCameraHist.hh"
#include "wfCamSim.hh"
#include "triggerSim.hh"

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

void anaPCA::Loop(TString histOut){
  //
  const unsigned int nn_fadc_point = 75;
  const unsigned int nn_PMT_channels = 7987;
  //
  Int_t fadc_sum_offset = 15;
  //Float_t fadc_amplitude = 8.25;
  Int_t fadc_MHz = 1024;
  Int_t fadc_offset = 300;
  //Int_t fadc_offset = 0;
  Float_t fadc_sample_in_ns = 1000.0/fadc_MHz;
  Float_t time_offset = fadc_sum_offset*fadc_sample_in_ns;
  //Float_t NGB_rate_in_MHz = 0.0001;
  //Float_t fadc_electronic_noise_RMS = 0.01;
  Float_t NGB_rate_in_MHz = 386.0;
  Float_t fadc_electronic_noise_RMS = 3.94;
  //Float_t NGB_rate_in_MHz = 0.0;
  //Float_t fadc_electronic_noise_RMS = 0.0;
  //
  Int_t n_ev_cuts = 0;
  //
  vector<vector<Int_t>> wfcam(nn_PMT_channels, vector<Int_t>(nn_fadc_point));
  //
  TRandom3 *rnd = new TRandom3(123123);
  wfCamSim *wf = new wfCamSim( rnd, "Template_CTA_SiPM.txt", "spe.dat",
  			       nn_fadc_point, nn_PMT_channels, fadc_offset, fadc_sample_in_ns, NGB_rate_in_MHz, fadc_electronic_noise_RMS);
  wf->print_wfCamSim_configure();
  wf->simulate_cam_event(nn_fadc_point,
			 nn_PMT_channels,
			 wfcam,
			 0,
			 time_offset,
			 0,
			 NULL,
			 NULL);
  //
  sipmCameraHist *sipm_cam_trg = new sipmCameraHist("sipm_cam_trg","sipm_cam_trg","pixel_mapping.csv",0);
  triggerSim *trg_sim = new triggerSim(sipm_cam_trg);
  //
  sipmCameraHist *sipm_cam = new sipmCameraHist("sipm_cam","sipm_cam","pixel_mapping.csv",0);
  //
  TH1D *h1_digital_sum = new TH1D("h1_digital_sum","h1_digital_sum",1000,0.0,1000);
  TH1D *h1_digital_sum_3ns = new TH1D("h1_digital_sum_3ns","h1_digital_sum_3ns",1000,0.0,1000);
  TH1D *h1_digital_sum_5ns = new TH1D("h1_digital_sum_5ns","h1_digital_sum_5ns",1000,0.0,1000);
  TH1D *h1_fadc_val = new TH1D("h1_fadc_val","h1_fadc_val",1000,0.0,1000);
  //
  TH1D *h1_nphotons = new TH1D("h1_nphotons","h1_nphotons",1000,0.0,10000);
  TH1D *h1_n_pe = new TH1D("h1_n_pe","h1_n_pe",1000,0.0,1000);
  TH1D *h1_n_pixels = new TH1D("h1_n_pixels","h1_n_pixels",1000,0.0,10000);
  //
  TH1D *h1_azimuth = new TH1D("h1_azimuth","azimuth",400,-4,4);
  TH1D *h1_altitude = new TH1D("h1_altitude","h1_altitude",400,4,4);
  TH1D *h1_h_first_int = new TH1D("h1_h_first_int","h1_h_first_int",4000,0.0,100000);
  TH1D *h1_xmax = new TH1D("h1_xmax","h1_xmax",400,0.0,1000);
  TH1D *h1_hmax = new TH1D("h1_hmax","h1_hmax",400,0.0,100000);
  TH1D *h1_emax = new TH1D("h1_emax","h1_emax",400,0.0,1000);
  TH1D *h1_cmax = new TH1D("h1_cmax","h1_cmax",400,0.0,1000);
  //
  TH1D *h1_xcore = new TH1D("h1_xcore","h1_xcore",400,-1000,1000);
  TH1D *h1_ycore = new TH1D("h1_ycore","h1_ycore",400,-1000,1000);
  TH1D *h1_r_core = new TH1D("h1_r_core","h1_r_core",400, 0,1100);
  TH1D *h1_theta_core = new TH1D("h1_theta_core","h1_theta_core",400,-0.1,2*TMath::Pi()+0.1);
  //
  TH2D *h2_theta_core_vs_pix_theta = new TH2D("h2_theta_core_vs_pix_theta","h2_theta_core_vs_pix_theta",
					      100,-0.1,2*TMath::Pi()+0.1,
					      100,-0.1,2*TMath::Pi()+0.1); 
  //
  TH1D *h1_pix_r = new TH1D("h1_pix_r","h1_pix_r",400,-2,2);
  //
  TH2D *h2_ycore_vs_xcore = new TH2D("h2_ycore_vs_xcore","h2_ycore_vs_xcore",220,-1100,1100,220,-1100,1100);
  TH2D *h2_ycore_vs_xcore_w = new TH2D("h2_ycore_vs_xcore_w","h2_ycore_vs_xcore_w",220,-1100,1100,220,-1100,1100);
  TH2D *h2_ycore_vs_xcore_norm = new TH2D("h2_ycore_vs_xcore_norm","h2_ycore_vs_xcore_norm",220,-1100,1100,220,-1100,1100);
  TH2D *h2_ycore_vs_xcore_cut = new TH2D("h2_ycore_vs_xcore_cut","h2_ycore_vs_xcore_cut",220,-1100,1100,220,-1100,1100);
  //
  Double_t x0_LST01 = -70.93;
  Double_t y0_LST01 = -52.07;
  //
  Double_t r_core = 0.0;
  Double_t theta_core = 0.0;
  //
  Long64_t nentries = fChain->GetEntriesFast();
  cout<<"nentries = "<<nentries<<endl;
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    //for (Long64_t jentry=0; jentry<100000;jentry++) {
    if(jentry%10 == 0)
      cout<<jentry<<endl;
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    //if (ientry > 10) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    //
    getCore_rel_R_theta( x0_LST01, y0_LST01, xcore, ycore, r_core, theta_core);
    _r_core = r_core;
    _theta_core = theta_core;
    if(cuts()){
      n_ev_cuts++;
      if(n_ev_cuts>100)
	break;
      if(cuts()){
	h2_ycore_vs_xcore_cut->Fill(xcore,ycore);
      }
      h1_nphotons->Fill(nphotons);
      h1_n_pe->Fill(n_pe);
      h1_n_pixels->Fill(n_pixels);
      //
      h1_azimuth->Fill(azimuth);
      h1_altitude->Fill(altitude);
      h1_h_first_int->Fill(h_first_int);
      //
      h1_xmax->Fill(xmax);
      h1_hmax->Fill(hmax);
      h1_emax->Fill(emax);
      h1_cmax->Fill(cmax);
      //
      h1_xcore->Fill(xcore);
      h1_ycore->Fill(ycore);
      //
      h1_r_core->Fill(r_core);
      h1_theta_core->Fill(theta_core);
      //
      for(Int_t jj = 0;jj<n_pe;jj++){
	sipm_cam->get_pixel_vec().at(pe_chID[jj]).pix_phi;
	h1_pix_r->Fill(sipm_cam->get_pixel_vec().at(pe_chID[jj]).pix_r);
	if(sipm_cam->get_pixel_vec().at(pe_chID[jj]).pix_r>0.2)
	  h2_theta_core_vs_pix_theta->Fill(theta_core, sipm_cam->get_pixel_vec().at(pe_chID[jj]).pix_phi);
      }
      //
      //if(cuts()){
      h2_ycore_vs_xcore->Fill(xcore,ycore);
      h2_ycore_vs_xcore_w->Fill(xcore,ycore,n_pe);
      //
      ///*
      wf->simulate_cam_event(nn_fadc_point,
      			   nn_PMT_channels,
      			   wfcam,
      			   ev_time,
      			   time_offset,
      			   n_pe,
      			   pe_chID,
      			   pe_time);
      //
      trg_sim->get_trigger( wfcam, h1_digital_sum, h1_digital_sum_3ns, h1_digital_sum_5ns, h1_fadc_val);
      std::vector<std::vector<unsigned int>> trg_vec;
      //*/
      //sipm_cam->Fill_pe(n_pe, pe_chID);
      //sipm_cam->Fill_pe(n_pe, pe_chID, -theta_core);
      sipm_cam->Fill_pe_center(n_pe, pe_chID);
      //sipm_cam->Fill(0,0);
      //printEv();    
    }
  }
  //
  //----------------------
  TH2D_divide( h2_ycore_vs_xcore_w, h2_ycore_vs_xcore, h2_ycore_vs_xcore_norm);
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
  h1_nphotons->Write();
  h1_n_pe->Write();
  h1_n_pixels->Write();
  //
  h1_azimuth->Write();
  h1_altitude->Write();
  h1_h_first_int->Write();
  h1_xmax->Write();
  h1_hmax->Write();
  h1_emax->Write();
  h1_cmax->Write();
  //
  h1_xcore->Write();
  h1_ycore->Write();
  h1_r_core->Write();
  h1_theta_core->Write();
  h2_theta_core_vs_pix_theta->Write();
  h1_pix_r->Write();
  //
  h2_ycore_vs_xcore->Write();
  h2_ycore_vs_xcore_w->Write();
  h2_ycore_vs_xcore_norm->Write();
  h2_ycore_vs_xcore_cut->Write();
  //
  h1_digital_sum->Write();
  h1_digital_sum_3ns->Write();
  h1_digital_sum_5ns->Write();
  h1_fadc_val->Write();
  //
  sipm_cam->Write();
  //
  rootFile->Close();
  //
  cout<<"n_ev_cuts = "<<n_ev_cuts<<endl;
}

bool anaPCA::cuts(){
  //if((_theta_core*180/TMath::Pi()<10.0) ||
  //  ((_theta_core*180/TMath::Pi()>125.0) && (_theta_core*180/TMath::Pi()<135.0)) ||
  // ((_theta_core*180/TMath::Pi()>245.0) && (_theta_core*180/TMath::Pi()<255.0)))
  //return true;
  //if(energy>0.1)
  if(n_pe == 100)
    return true;
  //if(n_pe>10 && n_pe<40)
  //return true;
  //if(n_pe>10)
  //return true;
  //if(n_pe>40)
  //return true;
  //if(energy>1.0)
  //return true;
  //if(n_pe<20 && n_pe>10)
  //if(n_pe>100)
  //if(xcore>-10 && xcore<10)
  //if(ycore>-10 && ycore<10)
  //return true;
  return false; 
}
