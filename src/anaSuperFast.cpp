//my
#include "anaSuperFast.hh"
#include "anabase.hh"
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
  Double_t val_Emin = 1.0;      // GeV
  Double_t val_Emax = 100000;   // GeV
  Int_t val_N_bins_E = 25;
  //
  Double_t val_Thetamin = 0.0;  //deg
  Double_t val_Thetamax = 10.0; //deg
  Int_t val_N_bins_t = 10;
  //
  evstHist *evH_sim_more_then_100pe_all = new evstHist("evH_sim_more_then_100pe_all","evH_sim_more_then_100pe_all",
						       val_Emin, val_Emax, val_N_bins_E,
						       val_Thetamin, val_Thetamax, val_N_bins_t);
  evstHist *evH_sim_more_then_100pe_100peminperch = new evstHist("evH_sim_more_then_100pe_100peminperch","evH_sim_more_then_100pe_100peminperch",
						   val_Emin, val_Emax, val_N_bins_E,
						   val_Thetamin, val_Thetamax, val_N_bins_t);
  //
  evstHist *evH_sim_more_then_100pe_all_soft = new evstHist("evH_sim_more_then_100pe_all_soft","evH_sim_more_then_100pe_all_soft",
							    val_Emin, val_Emax, val_N_bins_E,
							    val_Thetamin, val_Thetamax, val_N_bins_t);
  evstHist *evH_sim_more_then_100pe_100peminperch_soft = new evstHist("evH_sim_more_then_100pe_100peminperch_soft","evH_sim_more_then_100pe_100peminperch_soft",
								      val_Emin, val_Emax, val_N_bins_E,
								      val_Thetamin, val_Thetamax, val_N_bins_t); 
  //
  evstHist *evH_sim_more_then_100pe_100peminperch_eff = new evstHist("evH_sim_more_then_100pe_100peminperch_eff","evH_sim_more_then_100pe_100peminperch_eff",
								     val_Emin, val_Emax, val_N_bins_E,
								     val_Thetamin, val_Thetamax, val_N_bins_t);
  evstHist *evH_sim_more_then_100pe_100peminperch_soft_eff = new evstHist("evH_sim_more_then_100pe_100peminperch_soft_eff","evH_sim_more_then_100pe_100peminperch_soft_eff",
									  val_Emin, val_Emax, val_N_bins_E,
									  val_Thetamin, val_Thetamax, val_N_bins_t);
  
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
  //TH1D *h1_if_more_than_100pe = new TH1D("h1_if_more_than_100pe","h1_if_more_than_100pe",2,-0.5,1.5);
  TH1D *h1_npe_per_ch_max = new TH1D("h1_npe_per_ch_max","h1_npe_per_ch_max",100000,0.0,100000);
  TH1D *h1_npe_per_ch_max_integrated = new TH1D("h1_npe_per_ch_max_integrated","h1_npe_per_ch_max_integrated",100000,0.0,100000);
  TH1D *h1_npe_per_ch_max_nsim_integrated = new TH1D("h1_npe_per_ch_max_nsim_integrated","h1_npe_per_ch_max_nsim_integrated",100000,0.0,100000);
  //
  TH1D *h1_npe_per_ch_max_norm = new TH1D("h1_npe_per_ch_max_norm","h1_npe_per_ch_max_norm",100000,0.0,100000);
  TH1D *h1_npe_per_ch_max_norm_integrated = new TH1D("h1_npe_per_ch_max_norm_integrated","h1_npe_per_ch_max_norm_integrated",100000,0.0,100000);
  TH1D *h1_npe_per_ch_max_normsim_integrated = new TH1D("h1_npe_per_ch_max_normsim_integrated","h1_npe_per_ch_max_normsim_integrated",100000,0.0,100000);
  //
  TH1D *h1_npe_per_ch_max_norm_soft = new TH1D("h1_npe_per_ch_max_norm_soft","h1_npe_per_ch_max_norm_soft",100000,0.0,100000);
  TH1D *h1_npe_per_ch_max_norm_soft_integrated = new TH1D("h1_npe_per_ch_max_norm_soft_integrated","h1_npe_per_ch_max_norm_soft_integrated",100000,0.0,100000);
  TH1D *h1_npe_per_ch_max_normsim_soft_integrated = new TH1D("h1_npe_per_ch_max_normsim_soft_integrated","h1_npe_per_ch_max_normsim_soft_integrated",100000,0.0,100000);
  //
  Double_t npe_per_ch_max = 0.0;
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
  //
  Int_t npe_inchannel[nChannels];
  //
  Int_t npe_per_ch_saturation = 100;
  Int_t n_ch_saturation = 0;
  //
  Double_t normsim_ev = 0.0;
  //
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
    //theta_p_t_deg = get_theta_p_t()*180/TMath::Pi();
    //
    h1_energy->Fill(energy);
    h1_nphotons->Fill(nphotons);
    h1_n_pe->Fill(n_pe);
    h1_n_pixels->Fill(n_pixels);
    //
    normsim_ev++;
    //
    if(n_pe>100){
      for(Int_t jj = 0;jj<nChannels;jj++)
	npe_inchannel[jj] = 0;
      //
      for( Int_t jj = 0; jj<n_pe; jj++)
	if(pe_chID[jj]>=0 && pe_chID[jj]<nChannels)
	  npe_inchannel[pe_chID[jj]] = npe_inchannel[pe_chID[jj]] + 1;
      //
      npe_per_ch_max = 0.0;
      n_ch_saturation = 0;
      for(Int_t jj = 0;jj<nChannels;jj++){
	if(npe_inchannel[jj]>npe_per_ch_max)
	  npe_per_ch_max = npe_inchannel[jj];
	if(npe_inchannel[jj]>=npe_per_ch_saturation)
	  n_ch_saturation++;	
      }
      //
      h1_npe_per_ch_max->Fill(npe_per_ch_max);
      h1_npe_per_ch_max_norm->Fill( npe_per_ch_max, evstHist::get_Weight_ETeV(energy));
      h1_npe_per_ch_max_norm_soft->Fill( npe_per_ch_max, evstHist::get_Weight_ETeV_soft(energy));
      //
      //
      evH_sim_more_then_100pe_all->get_E_hist()->Fill(energy*1000.0, evstHist::get_Weight_ETeV(energy));
      if(n_ch_saturation>=10)
	if(npe_per_ch_max>=npe_per_ch_saturation)
	  evH_sim_more_then_100pe_100peminperch->get_E_hist()->Fill(energy*1000.0, evstHist::get_Weight_ETeV(energy));
      //
      //
      evH_sim_more_then_100pe_all_soft->get_E_hist()->Fill(energy*1000.0, evstHist::get_Weight_ETeV_soft(energy));
      if(n_ch_saturation>=10)
	if(npe_per_ch_max>=npe_per_ch_saturation)
	  evH_sim_more_then_100pe_100peminperch_soft->get_E_hist()->Fill(energy*1000.0, evstHist::get_Weight_ETeV_soft(energy));      
    }
    //}
      //h1_if_more_than_100pe->Fill(0.0, evstHist::get_Weight_ETeV(energy));
    //}
    //
    //if(n_pe>50){
      //
      /*
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
      */
    //}
    //
  }
  //  
  //
  //TH2D_divide( h2_ycore_vs_xcore_km_ev_time_w, h2_ycore_vs_xcore_km_norm, h2_ycore_vs_xcore_km_ev_time);  
  //TH1D_divide( h1_xcore_km_ev_time_cut_ycore_w, h1_xcore_km_ev_time_cut_ycore_norm, h1_xcore_km_ev_time_cut_ycore);
  //
  anabase::TH1D_divide( evH_sim_more_then_100pe_100peminperch->get_E_hist(),
			evH_sim_more_then_100pe_all->get_E_hist(),
			evH_sim_more_then_100pe_100peminperch_eff->get_E_hist());
  anabase::TH1D_divide( evH_sim_more_then_100pe_100peminperch_soft->get_E_hist(),
			evH_sim_more_then_100pe_all_soft->get_E_hist(),
			evH_sim_more_then_100pe_100peminperch_soft_eff->get_E_hist());
  //
  get_cumulative(h1_npe_per_ch_max,h1_npe_per_ch_max_integrated);
  get_cumulative(h1_npe_per_ch_max_norm,h1_npe_per_ch_max_norm_integrated);
  get_cumulative(h1_npe_per_ch_max_norm_soft,h1_npe_per_ch_max_norm_soft_integrated);
  //
  get_cumulative(h1_npe_per_ch_max,          h1_npe_per_ch_max_nsim_integrated, normsim_ev);
  get_cumulative(h1_npe_per_ch_max_norm,     h1_npe_per_ch_max_normsim_integrated, normsim_ev);
  get_cumulative(h1_npe_per_ch_max_norm_soft,h1_npe_per_ch_max_normsim_soft_integrated, normsim_ev);
  //
  save_hist_to_csv("h1_npe_per_ch_max_norm_integrated_E2p7spectrum.csv", h1_npe_per_ch_max_norm_integrated);
  save_hist_to_csv("h1_npe_per_ch_max_norm_integrated_E2p2spectrum.csv", h1_npe_per_ch_max_norm_soft_integrated);
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
  h1_npe_per_ch_max->Write();
  h1_npe_per_ch_max_norm->Write();
  h1_npe_per_ch_max_norm_soft->Write();  
  //
  h1_npe_per_ch_max_integrated->Write();
  h1_npe_per_ch_max_norm_integrated->Write();
  h1_npe_per_ch_max_norm_soft_integrated->Write();
  //
  //
  h1_npe_per_ch_max_nsim_integrated->Write();
  h1_npe_per_ch_max_normsim_integrated->Write();
  h1_npe_per_ch_max_normsim_soft_integrated->Write();
  //
  //
  evH_sim_more_then_100pe_all->get_E_hist()->Write();
  evH_sim_more_then_100pe_100peminperch->get_E_hist()->Write();
  evH_sim_more_then_100pe_all_soft->get_E_hist()->Write();
  evH_sim_more_then_100pe_100peminperch_soft->get_E_hist()->Write();
  //
  evH_sim_more_then_100pe_100peminperch_eff->get_E_hist()->Write();
  evH_sim_more_then_100pe_100peminperch_soft_eff->get_E_hist()->Write();  
  //
  //
  //
  //
  //
  //
  //h1_azimuth_deg->Write();
  //h1_altitude_deg->Write();
  //
  //h1_xcore_km->Write();
  //h1_ycore_km->Write();
  //
  //h1_theta_p_t_deg->Write();
  //
  //h1_ev_time->Write();
  //h1_ev_time_cut_ycore_xcore->Write();
  //h1_pe_time->Write();
  //h1_pe_time_shift->Write();
  //
  //h2_ycore_vs_xcore_km_norm->Write();
  //h2_ycore_vs_xcore_km_ev_time_w->Write();
  //h2_ycore_vs_xcore_km_ev_time->Write();
  //
  //h1_xcore_km_ev_time_cut_ycore_w->Write();
  //h1_xcore_km_ev_time_cut_ycore_norm->Write();
  //h1_xcore_km_ev_time_cut_ycore->Write();  
  //
  //h1_if_more_than_100pe->Write();
  //
  rootFile->Close();
}

Double_t anaSuperFast::get_theta_p_t(){
  TVector3 v_prot;
  v_prot.SetMagThetaPhi(1.0,TMath::Pi()/2.0-altitude,TMath::Pi() - azimuth);
  TVector3 v_prot_inv(v_prot.x(),v_prot.y(),v_prot.z());
  return TMath::ACos(v_prot_inv.Dot(_v_det)/v_prot_inv.Mag()/_v_det.Mag());
}

void anaSuperFast::get_cumulative( TH1D *h1, TH1D *h1_int){
  Double_t norm = h1->Integral(1,h1->GetNbinsX());
  get_cumulative( h1, h1_int, norm);
}

void anaSuperFast::get_cumulative( TH1D *h1, TH1D *h1_int, Double_t norm){
  for(Int_t i = 1;i<=h1->GetNbinsX();i++){
    Double_t val = h1->Integral(i,h1->GetNbinsX());
    h1_int->SetBinContent(i,val/norm);
  }
}

void anaSuperFast::save_hist_to_csv(TString outfilename, TH1D *h1){
  ofstream outfile;
  outfile.open(outfilename.Data());
  //
  outfile<<"x,y"<<endl;
  for(Int_t i = 1;i<=h1->GetNbinsX();i++){
    Double_t val = h1->GetBinContent(i);
    Double_t valx = h1->GetBinCenter(i);
    outfile<<valx<<","<<val<<endl;
  }
  outfile.close();
}
