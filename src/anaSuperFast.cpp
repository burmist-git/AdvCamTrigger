//my
#include "anaSuperFast.hh"
#include "sipmCameraHist.hh"
#include "sipmCameraHistCropped.hh"
#include "wfCamSim.hh"
#include "triggerSim.hh"
#include "anaConf.hh"
#include "evstHist.hh"

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
#include "TPrincipal.h"
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

anaSuperFast::anaSuperFast(TString fileList, TString anaConfFile) : anashort(fileList)
{
  _anaConf.readFromFile(anaConfFile);
  _v_det.SetXYZ(1.0*TMath::Sin(20.0/180.0*TMath::Pi()),0,1.0*TMath::Cos(20.0/180.0*TMath::Pi()));
}

void anaSuperFast::Loop(TString histOut){
  //
  cout<<"anaSuperFast::Loop"<<endl;
  _anaConf.printInfo();
  //assert(0);
  //
  TH1D *h1_energy = new TH1D("h1_energy","h1_energy",110000,0.0,110.0);
  //
  TH1D *h1_nphotons = new TH1D("h1_nphotons","h1_nphotons",10000,0.0,1000000);
  TH1D *h1_n_pe = new TH1D("h1_n_pe","h1_n_pe",1001,-0.5,10000.5);
  TH1D *h1_n_pixels = new TH1D("h1_n_pixels","h1_n_pixels",1000,0.0,10000);
  //
  TH1D *h1_theta_p_t_deg = new TH1D("theta_p_t_deg","theta_p_t_deg",4000,-370,370);
  TH1D *h1_azimuth_deg = new TH1D("h1_azimuth_deg","azimuth_deg",400,0.0,360);
  TH1D *h1_altitude_deg = new TH1D("h1_altitude_deg","h1_altitude_deg",400,0.0,180);
  //
  TH1D *h1_xcore_km = new TH1D("h1_xcore_km","h1_xcore_km",200,-1.8,1.8);
  TH1D *h1_ycore_km = new TH1D("h1_ycore_km","h1_ycore_km",200,-1.8,1.8);
  //
  TH2D *h2_ycore_vs_xcore_km_norm = new TH2D("h2_ycore_vs_xcore_km_norm","h2_ycore_vs_xcore_km_norm",
					     200,-1.8,1.8,200,-1.8,1.8);
  TH2D *h2_ycore_vs_xcore_km_ev_time_w = new TH2D("h2_ycore_vs_xcore_km_ev_time_w","h2_ycore_vs_xcore_km_ev_time_w",
						  200,-1.8,1.8,200,-1.8,1.8);
  TH2D *h2_ycore_vs_xcore_km_ev_time = new TH2D("h2_ycore_vs_xcore_km_ev_time","h2_ycore_vs_xcore_km_ev_time",
						200,-1.8,1.8,200,-1.8,1.8);
  //
  TH1D *h1_xcore_km_ev_time_cut_ycore = new TH1D("h1_xcore_km_ev_time_cut_ycore","h1_xcore_km_ev_time_cut_ycore", 200,-1.8,1.8);
  TH1D *h1_xcore_km_ev_time_cut_ycore_w = new TH1D("h1_xcore_km_ev_time_cut_ycore_w","h1_xcore_km_ev_time_cut_ycore_w", 200,-1.8,1.8);
  TH1D *h1_xcore_km_ev_time_cut_ycore_norm = new TH1D("h1_xcore_km_ev_time_cut_ycore_norm","h1_xcore_km_ev_time_cut_ycore_norm", 200,-1.8,1.8);
  //
  TH1D *h1_ev_time = new TH1D("h1_ev_time","h1_ev_time",4000,-4000,4000);
  TH1D *h1_ev_time_cut_ycore_xcore = new TH1D("h1_ev_time_cut_ycore_xcore","h1_ev_time_cut_ycore_xcore",4000,-4000,4000);
  TH1D *h1_pe_time = new TH1D("h1_pe_time","h1_pe_time",4000,-4000,4000);
  TH1D *h1_pe_time_shift = new TH1D("h1_pe_time_shift","h1_pe_time_shift",4000,-4000,4000);
  //
  Double_t x0_LST01 = -70.93;
  Double_t y0_LST01 = -52.07;
  //
  Double_t r_core = 0.0;
  Double_t theta_core = 0.0;
  //
  Double_t azimuth_deg;
  Double_t altitude_deg;
  //
  Double_t theta_p_t_deg;
  //
  Int_t fadc_sum_offset = 15;
  Int_t fadc_MHz = 1024;
  Float_t fadc_sample_in_ns = 1000.0/fadc_MHz;
  Float_t time_offset = fadc_sum_offset*fadc_sample_in_ns;  
  //
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  if(_anaConf.nentries_max>0)
    nentries = _anaConf.nentries_max;
  cout<<"nentries = "<<nentries<<endl;
  //
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    if(jentry%_anaConf.jentry_modulo == 0)
      cout<<jentry<<endl;
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    //if (ientry > 10) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    //
    //
    theta_p_t_deg = get_theta_p_t()*180/TMath::Pi();
    //
    if(n_pe>50){
      //
      getCore_rel_R_theta( x0_LST01, y0_LST01, xcore, ycore, r_core, theta_core);
      //
      azimuth_deg = azimuth*180.0/TMath::Pi();
      altitude_deg = altitude*180.0/TMath::Pi();
      //
      h2_ycore_vs_xcore_km_norm->Fill( xcore/1000.0, ycore/1000.0);
      h2_ycore_vs_xcore_km_ev_time_w->Fill( xcore/1000.0, ycore/1000.0, ev_time);
      //
      if(ycore<(y0_LST01+100.0)){
	if(ycore>(y0_LST01-100.0)){
	  h1_xcore_km_ev_time_cut_ycore_w->Fill(xcore/1000.0,ev_time);
	  h1_xcore_km_ev_time_cut_ycore_norm->Fill(xcore/1000.0);
	}
      }
      //
      if(ycore<(y0_LST01+5)){
	if(ycore>(y0_LST01-5)){
	  if(xcore<(x0_LST01+5)){
	    if(xcore>(x0_LST01-5)){
	      h1_ev_time_cut_ycore_xcore->Fill(ev_time);
	    }
	  }
	}
      }
      //
      h1_theta_p_t_deg->Fill(theta_p_t_deg);
      //
      h1_energy->Fill(energy);
      h1_nphotons->Fill(nphotons);
      h1_n_pe->Fill(n_pe);
      h1_n_pixels->Fill(n_pixels);
      //
      h1_azimuth_deg->Fill(azimuth_deg);
      h1_altitude_deg->Fill(altitude_deg);
      //
      h1_xcore_km->Fill(xcore/1000.0);
      h1_ycore_km->Fill(ycore/1000.0);
      //
      h1_ev_time->Fill(ev_time);
      //
      for(Int_t jj = 0;jj<n_pe;jj++){
	h1_pe_time->Fill(pe_time[jj]);
	h1_pe_time_shift->Fill(pe_time[jj]-ev_time+time_offset);
      }
    }
    //
  }
  //  
  //
  TH2D_divide( h2_ycore_vs_xcore_km_ev_time_w, h2_ycore_vs_xcore_km_norm, h2_ycore_vs_xcore_km_ev_time);  
  TH1D_divide( h1_xcore_km_ev_time_cut_ycore_w, h1_xcore_km_ev_time_cut_ycore_norm, h1_xcore_km_ev_time_cut_ycore);
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
  h1_energy->Write();
  h1_nphotons->Write();
  h1_n_pe->Write();
  h1_n_pixels->Write();
  //
  h1_azimuth_deg->Write();
  h1_altitude_deg->Write();
  //
  h1_xcore_km->Write();
  h1_ycore_km->Write();
  //
  h1_theta_p_t_deg->Write();
  //
  h1_ev_time->Write();
  h1_ev_time_cut_ycore_xcore->Write();
  h1_pe_time->Write();
  h1_pe_time_shift->Write();
  //
  h2_ycore_vs_xcore_km_norm->Write();
  h2_ycore_vs_xcore_km_ev_time_w->Write();
  h2_ycore_vs_xcore_km_ev_time->Write();
  //
  h1_xcore_km_ev_time_cut_ycore_w->Write();
  h1_xcore_km_ev_time_cut_ycore_norm->Write();
  h1_xcore_km_ev_time_cut_ycore->Write();  
  //
  rootFile->Close();
}

Double_t anaSuperFast::get_theta_p_t(){
  TVector3 v_prot;
  v_prot.SetMagThetaPhi(1.0,TMath::Pi()/2.0-altitude,TMath::Pi() - azimuth);
  TVector3 v_prot_inv(v_prot.x(),v_prot.y(),v_prot.z());
  return TMath::ACos(v_prot_inv.Dot(_v_det)/v_prot_inv.Mag()/_v_det.Mag());
}
