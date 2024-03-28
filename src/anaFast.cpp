//my
#include "anaFast.hh"
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

anaFast::anaFast(TString fileList, TString anaConfFile) : anashort(fileList)
{
  _anaConf.readFromFile(anaConfFile);
  _v_det.SetXYZ(1.0*TMath::Sin(20.0/180.0*TMath::Pi()),0,1.0*TMath::Cos(20.0/180.0*TMath::Pi()));
}

void anaFast::Loop_scan( TString histOut, TString simtel_all_dat, TString flux_dat){
  //
  cout<<"Loop_scan"<<endl;  
  _anaConf.printInfo();
  //
  Double_t val_Emin = 1.0;    // GeV
  Double_t val_Emax = 100000; // GeV
  Int_t val_N_bins_E = 25;
  //
  Double_t val_Thetamin = 0.0;  //deg
  Double_t val_Thetamax = 10.0; //deg
  Int_t val_N_bins_t = 10;
  //
  Int_t n_ev_cuts = 0;
  //

  evstHist *evH_nhe = new evstHist("evH_nhe","evH_nhe",
				   val_Emin, val_Emax, val_N_bins_E,
				   val_Thetamin, val_Thetamax, val_N_bins_t);
  evH_nhe->init_h1_arb_v( "h1_npe", "h1_npe", 50001,-0.5,50000.5);
  evstHist *evH_nhe_w = new evstHist("evH_nhe_w","evH_nhe_w",
				     val_Emin, val_Emax, val_N_bins_E,
				     val_Thetamin, val_Thetamax, val_N_bins_t);
  evH_nhe_w->init_h1_arb_v( "h1_npe_w", "h1_npe_w", 10001,-0.5,10000.5);
  //
  evstHist *evH = new evstHist("evH","evH",
			       val_Emin, val_Emax, val_N_bins_E,
			       val_Thetamin, val_Thetamax, val_N_bins_t);
  //
  evstHist *evH_10npe = new evstHist("evH_10npe","evH_10npe",
				     val_Emin, val_Emax, val_N_bins_E,
				     val_Thetamin, val_Thetamax, val_N_bins_t);
  evstHist *evH_20npe = new evstHist("evH_20npe","evH_20npe",
				     val_Emin, val_Emax, val_N_bins_E,
				     val_Thetamin, val_Thetamax, val_N_bins_t);
  evstHist *evH_30npe = new evstHist("evH_30npe","evH_30npe",
				     val_Emin, val_Emax, val_N_bins_E,
				     val_Thetamin, val_Thetamax, val_N_bins_t);
  evstHist *evH_40npe = new evstHist("evH_40npe","evH_40npe",
				     val_Emin, val_Emax, val_N_bins_E,
				     val_Thetamin, val_Thetamax, val_N_bins_t);
  evstHist *evH_50npe = new evstHist("evH_50npe","evH_50npe",
				     val_Emin, val_Emax, val_N_bins_E,
				     val_Thetamin, val_Thetamax, val_N_bins_t);
  evstHist *evH_60npe = new evstHist("evH_60npe","evH_60npe",
				     val_Emin, val_Emax, val_N_bins_E,
				     val_Thetamin, val_Thetamax, val_N_bins_t);
  evstHist *evH_70npe = new evstHist("evH_70npe","evH_70npe",
				     val_Emin, val_Emax, val_N_bins_E,
				     val_Thetamin, val_Thetamax, val_N_bins_t);
  evstHist *evH_80npe = new evstHist("evH_80npe","evH_80npe",
				     val_Emin, val_Emax, val_N_bins_E,
				     val_Thetamin, val_Thetamax, val_N_bins_t);
  evstHist *evH_90npe = new evstHist("evH_90npe","evH_90npe",
				     val_Emin, val_Emax, val_N_bins_E,
				     val_Thetamin, val_Thetamax, val_N_bins_t);
  evstHist *evH_100npe = new evstHist("evH_100npe","evH_100npe",
				     val_Emin, val_Emax, val_N_bins_E,
				     val_Thetamin, val_Thetamax, val_N_bins_t);
  //
  //
  evstHist *evH_10npe_w = new evstHist("evH_10npe_w","evH_10npe_w",
				     val_Emin, val_Emax, val_N_bins_E,
				     val_Thetamin, val_Thetamax, val_N_bins_t);
  evstHist *evH_20npe_w = new evstHist("evH_20npe_w","evH_20npe_w",
				     val_Emin, val_Emax, val_N_bins_E,
				     val_Thetamin, val_Thetamax, val_N_bins_t);
  evstHist *evH_30npe_w = new evstHist("evH_30npe_w","evH_30npe_w",
				     val_Emin, val_Emax, val_N_bins_E,
				     val_Thetamin, val_Thetamax, val_N_bins_t);
  evstHist *evH_40npe_w = new evstHist("evH_40npe_w","evH_40npe_w",
				     val_Emin, val_Emax, val_N_bins_E,
				     val_Thetamin, val_Thetamax, val_N_bins_t);
  evstHist *evH_50npe_w = new evstHist("evH_50npe_w","evH_50npe_w",
				     val_Emin, val_Emax, val_N_bins_E,
				     val_Thetamin, val_Thetamax, val_N_bins_t);
  evstHist *evH_60npe_w = new evstHist("evH_60npe_w","evH_60npe_w",
				     val_Emin, val_Emax, val_N_bins_E,
				     val_Thetamin, val_Thetamax, val_N_bins_t);
  evstHist *evH_70npe_w = new evstHist("evH_70npe_w","evH_70npe_w",
				     val_Emin, val_Emax, val_N_bins_E,
				     val_Thetamin, val_Thetamax, val_N_bins_t);
  evstHist *evH_80npe_w = new evstHist("evH_80npe_w","evH_80npe_w",
				     val_Emin, val_Emax, val_N_bins_E,
				     val_Thetamin, val_Thetamax, val_N_bins_t);
  evstHist *evH_90npe_w = new evstHist("evH_90npe_w","evH_90npe_w",
				     val_Emin, val_Emax, val_N_bins_E,
				     val_Thetamin, val_Thetamax, val_N_bins_t);
  evstHist *evH_100npe_w = new evstHist("evH_100npe_w","evH_100npe_w",
				     val_Emin, val_Emax, val_N_bins_E,
				     val_Thetamin, val_Thetamax, val_N_bins_t);
  //
  evstHist *evH_simtel = new evstHist("evH_simtel","evH_simtel",
				      val_Emin, val_Emax, val_N_bins_E,
				      val_Thetamin, val_Thetamax, val_N_bins_t);
  evstHist *evH_simtel_weight= new evstHist("evH_simtel_weight","evH_simtel_weight",
					    val_Emin, val_Emax, val_N_bins_E,
					    val_Thetamin, val_Thetamax, val_N_bins_t);
  evstHist *evH_tot_flux_nocore = new evstHist("evH_tot_flux_nocore","evH_tot_flux_nocore",
					       val_Emin, val_Emax, val_N_bins_E,
					       val_Thetamin, val_Thetamax, val_N_bins_t);
  evstHist *evH_tot_flux = new evstHist("evH_tot_flux","evH_tot_flux",
					val_Emin, val_Emax, val_N_bins_E,
					val_Thetamin, val_Thetamax, val_N_bins_t);
  evstHist *evH_eff = new evstHist("evH_eff","evH_eff",
				   val_Emin, val_Emax, val_N_bins_E,
				   val_Thetamin, val_Thetamax, val_N_bins_t);
  //
  evstHist *evH_eff_10npe = new evstHist("evH_eff_10npe","evH_eff_10npe",
					 val_Emin, val_Emax, val_N_bins_E,
					 val_Thetamin, val_Thetamax, val_N_bins_t);
  evstHist *evH_eff_20npe = new evstHist("evH_eff_20npe","evH_eff_20npe",
					 val_Emin, val_Emax, val_N_bins_E,
					 val_Thetamin, val_Thetamax, val_N_bins_t);
  evstHist *evH_eff_30npe = new evstHist("evH_eff_30npe","evH_eff_30npe",
					 val_Emin, val_Emax, val_N_bins_E,
					 val_Thetamin, val_Thetamax, val_N_bins_t);
  evstHist *evH_eff_40npe = new evstHist("evH_eff_40npe","evH_eff_40npe",
					 val_Emin, val_Emax, val_N_bins_E,
					 val_Thetamin, val_Thetamax, val_N_bins_t);
  evstHist *evH_eff_50npe = new evstHist("evH_eff_50npe","evH_eff_50npe",
					 val_Emin, val_Emax, val_N_bins_E,
					 val_Thetamin, val_Thetamax, val_N_bins_t);
  evstHist *evH_eff_60npe = new evstHist("evH_eff_60npe","evH_eff_60npe",
					 val_Emin, val_Emax, val_N_bins_E,
					 val_Thetamin, val_Thetamax, val_N_bins_t);
  evstHist *evH_eff_70npe = new evstHist("evH_eff_70npe","evH_eff_70npe",
					 val_Emin, val_Emax, val_N_bins_E,
					 val_Thetamin, val_Thetamax, val_N_bins_t);
  evstHist *evH_eff_80npe = new evstHist("evH_eff_80npe","evH_eff_80npe",
					 val_Emin, val_Emax, val_N_bins_E,
					 val_Thetamin, val_Thetamax, val_N_bins_t);
  evstHist *evH_eff_90npe = new evstHist("evH_eff_90npe","evH_eff_90npe",
					 val_Emin, val_Emax, val_N_bins_E,
					 val_Thetamin, val_Thetamax, val_N_bins_t);
  evstHist *evH_eff_100npe = new evstHist("evH_eff_100npe","evH_eff_100npe",
					 val_Emin, val_Emax, val_N_bins_E,
					 val_Thetamin, val_Thetamax, val_N_bins_t);
  //
  //
  evstHist *evH_eff_10npe_w = new evstHist("evH_eff_10npe_w","evH_eff_10npe_w",
					 val_Emin, val_Emax, val_N_bins_E,
					 val_Thetamin, val_Thetamax, val_N_bins_t);
  evstHist *evH_eff_20npe_w = new evstHist("evH_eff_20npe_w","evH_eff_20npe_w",
					 val_Emin, val_Emax, val_N_bins_E,
					 val_Thetamin, val_Thetamax, val_N_bins_t);
  evstHist *evH_eff_30npe_w = new evstHist("evH_eff_30npe_w","evH_eff_30npe_w",
					 val_Emin, val_Emax, val_N_bins_E,
					 val_Thetamin, val_Thetamax, val_N_bins_t);
  evstHist *evH_eff_40npe_w = new evstHist("evH_eff_40npe_w","evH_eff_40npe_w",
					 val_Emin, val_Emax, val_N_bins_E,
					 val_Thetamin, val_Thetamax, val_N_bins_t);
  evstHist *evH_eff_50npe_w = new evstHist("evH_eff_50npe_w","evH_eff_50npe_w",
					 val_Emin, val_Emax, val_N_bins_E,
					 val_Thetamin, val_Thetamax, val_N_bins_t);
  evstHist *evH_eff_60npe_w = new evstHist("evH_eff_60npe_w","evH_eff_60npe_w",
					 val_Emin, val_Emax, val_N_bins_E,
					 val_Thetamin, val_Thetamax, val_N_bins_t);
  evstHist *evH_eff_70npe_w = new evstHist("evH_eff_70npe_w","evH_eff_70npe_w",
					 val_Emin, val_Emax, val_N_bins_E,
					 val_Thetamin, val_Thetamax, val_N_bins_t);
  evstHist *evH_eff_80npe_w = new evstHist("evH_eff_80npe_w","evH_eff_80npe_w",
					 val_Emin, val_Emax, val_N_bins_E,
					 val_Thetamin, val_Thetamax, val_N_bins_t);
  evstHist *evH_eff_90npe_w = new evstHist("evH_eff_90npe_w","evH_eff_90npe_w",
					 val_Emin, val_Emax, val_N_bins_E,
					 val_Thetamin, val_Thetamax, val_N_bins_t);
  evstHist *evH_eff_100npe_w = new evstHist("evH_eff_100npe_w","evH_eff_100npe_w",
					 val_Emin, val_Emax, val_N_bins_E,
					 val_Thetamin, val_Thetamax, val_N_bins_t);
  //
  //
  evstHist *evH_flux = new evstHist("evH_flux","evH_flux",
				    val_Emin, val_Emax, val_N_bins_E,
				    val_Thetamin, val_Thetamax, val_N_bins_t);
  evstHist *evH_flux_10npe = new evstHist("evH_flux_10npe","evH_flux_10npe",
					  val_Emin, val_Emax, val_N_bins_E,
					  val_Thetamin, val_Thetamax, val_N_bins_t);
  evstHist *evH_flux_20npe = new evstHist("evH_flux_20npe","evH_flux_20npe",
					  val_Emin, val_Emax, val_N_bins_E,
					  val_Thetamin, val_Thetamax, val_N_bins_t);
  evstHist *evH_flux_30npe = new evstHist("evH_flux_30npe","evH_flux_30npe",
					  val_Emin, val_Emax, val_N_bins_E,
					  val_Thetamin, val_Thetamax, val_N_bins_t);
  evstHist *evH_flux_40npe = new evstHist("evH_flux_40npe","evH_flux_40npe",
					  val_Emin, val_Emax, val_N_bins_E,
					  val_Thetamin, val_Thetamax, val_N_bins_t);
  evstHist *evH_flux_50npe = new evstHist("evH_flux_50npe","evH_flux_50npe",
					  val_Emin, val_Emax, val_N_bins_E,
					  val_Thetamin, val_Thetamax, val_N_bins_t);
  evstHist *evH_flux_60npe = new evstHist("evH_flux_60npe","evH_flux_60npe",
					  val_Emin, val_Emax, val_N_bins_E,
					  val_Thetamin, val_Thetamax, val_N_bins_t);
  evstHist *evH_flux_70npe = new evstHist("evH_flux_70npe","evH_flux_70npe",
					  val_Emin, val_Emax, val_N_bins_E,
					  val_Thetamin, val_Thetamax, val_N_bins_t);
  evstHist *evH_flux_80npe = new evstHist("evH_flux_80npe","evH_flux_80npe",
					  val_Emin, val_Emax, val_N_bins_E,
					  val_Thetamin, val_Thetamax, val_N_bins_t);
  evstHist *evH_flux_90npe = new evstHist("evH_flux_90npe","evH_flux_90npe",
					  val_Emin, val_Emax, val_N_bins_E,
					  val_Thetamin, val_Thetamax, val_N_bins_t);
  evstHist *evH_flux_100npe = new evstHist("evH_flux_100npe","evH_flux_100npe",
					  val_Emin, val_Emax, val_N_bins_E,
					  val_Thetamin, val_Thetamax, val_N_bins_t);
  //
  //
  evstHist *evH_flux_10npe_w = new evstHist("evH_flux_10npe_w","evH_flux_10npe_w",
					  val_Emin, val_Emax, val_N_bins_E,
					  val_Thetamin, val_Thetamax, val_N_bins_t);
  evstHist *evH_flux_20npe_w = new evstHist("evH_flux_20npe_w","evH_flux_20npe_w",
					  val_Emin, val_Emax, val_N_bins_E,
					  val_Thetamin, val_Thetamax, val_N_bins_t);
  evstHist *evH_flux_30npe_w = new evstHist("evH_flux_30npe_w","evH_flux_30npe_w",
					  val_Emin, val_Emax, val_N_bins_E,
					  val_Thetamin, val_Thetamax, val_N_bins_t);
  evstHist *evH_flux_40npe_w = new evstHist("evH_flux_40npe_w","evH_flux_40npe_w",
					  val_Emin, val_Emax, val_N_bins_E,
					  val_Thetamin, val_Thetamax, val_N_bins_t);
  evstHist *evH_flux_50npe_w = new evstHist("evH_flux_50npe_w","evH_flux_50npe_w",
					  val_Emin, val_Emax, val_N_bins_E,
					  val_Thetamin, val_Thetamax, val_N_bins_t);
  evstHist *evH_flux_60npe_w = new evstHist("evH_flux_60npe_w","evH_flux_60npe_w",
					  val_Emin, val_Emax, val_N_bins_E,
					  val_Thetamin, val_Thetamax, val_N_bins_t);
  evstHist *evH_flux_70npe_w = new evstHist("evH_flux_70npe_w","evH_flux_70npe_w",
					  val_Emin, val_Emax, val_N_bins_E,
					  val_Thetamin, val_Thetamax, val_N_bins_t);
  evstHist *evH_flux_80npe_w = new evstHist("evH_flux_80npe_w","evH_flux_80npe_w",
					  val_Emin, val_Emax, val_N_bins_E,
					  val_Thetamin, val_Thetamax, val_N_bins_t);
  evstHist *evH_flux_90npe_w = new evstHist("evH_flux_90npe_w","evH_flux_90npe_w",
					  val_Emin, val_Emax, val_N_bins_E,
					  val_Thetamin, val_Thetamax, val_N_bins_t);
  evstHist *evH_flux_100npe_w = new evstHist("evH_flux_100npe_w","evH_flux_100npe_w",
					  val_Emin, val_Emax, val_N_bins_E,
					  val_Thetamin, val_Thetamax, val_N_bins_t);
  //
  evH_simtel->LoadBinContent(simtel_all_dat.Data(), true);
  evH_simtel_weight->LoadBinContent(simtel_all_dat.Data(), true);
  evH_simtel_weight->Weight_TeV();
  evH_tot_flux->LoadBinContent(simtel_all_dat.Data(), true);
  evH_tot_flux->selfNorm();
  evH_tot_flux_nocore->LoadBinContent(flux_dat.Data());
  evH_tot_flux->Multiply_core(evH_tot_flux_nocore);
  //
  evstHist::PrintBinsInfo(evH->get_E_hist());
  //
  Double_t theta_p_t;
  Double_t theta_p_t_deg;
  TH1D *h1_n_pe = new TH1D("h1_n_pe","h1_n_pe",10001,-0.5,10000.5);
  //
  Double_t x0_LST01 = -70.93;
  Double_t y0_LST01 = -52.07;
  Double_t r_core;
  Double_t theta_core;
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
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    //
    theta_p_t = get_theta_p_t();
    theta_p_t_deg = theta_p_t*180/TMath::Pi();
    getCore_rel_R_theta( x0_LST01, y0_LST01, xcore, ycore, r_core, theta_core);
    //
    evH->Fill(theta_p_t_deg,energy*1000.0);
    evH->Fill_rcore(theta_p_t_deg,energy*1000.0,r_core);
    //
    evH_nhe->fill_h1_arb_v(energy*1000.0, theta_p_t_deg, r_core, n_pe);
    evH_nhe_w->fill_h1_arb_v(energy*1000.0, theta_p_t_deg, r_core, n_pe,evstHist::get_Weight_ETeV(energy));
    //
    h1_n_pe->Fill(n_pe);
    //
    if(n_pe>9){
      evH_10npe->Fill(theta_p_t_deg,energy*1000.0);
      evH_10npe->Fill_rcore(theta_p_t_deg,energy*1000.0,r_core);
      //
      evH_10npe_w->Fill(theta_p_t_deg,energy*1000.0, evstHist::get_Weight_ETeV(energy));
      evH_10npe_w->Fill_rcore(theta_p_t_deg,energy*1000.0,r_core, evstHist::get_Weight_ETeV(energy));
    }
    if(n_pe>19){
      evH_20npe->Fill(theta_p_t_deg,energy*1000.0);
      evH_20npe->Fill_rcore(theta_p_t_deg,energy*1000.0,r_core);
      //
      evH_20npe_w->Fill(theta_p_t_deg,energy*1000.0, evstHist::get_Weight_ETeV(energy));
      evH_20npe_w->Fill_rcore(theta_p_t_deg,energy*1000.0,r_core, evstHist::get_Weight_ETeV(energy));
    }
    if(n_pe>29){
      evH_30npe->Fill(theta_p_t_deg,energy*1000.0);
      evH_30npe->Fill_rcore(theta_p_t_deg,energy*1000.0,r_core);
      //
      evH_30npe_w->Fill(theta_p_t_deg,energy*1000.0, evstHist::get_Weight_ETeV(energy));
      evH_30npe_w->Fill_rcore(theta_p_t_deg,energy*1000.0,r_core, evstHist::get_Weight_ETeV(energy));
    }
    if(n_pe>39){
      evH_40npe->Fill(theta_p_t_deg,energy*1000.0);
      evH_40npe->Fill_rcore(theta_p_t_deg,energy*1000.0,r_core);
      //
      evH_40npe_w->Fill(theta_p_t_deg,energy*1000.0, evstHist::get_Weight_ETeV(energy));
      evH_40npe_w->Fill_rcore(theta_p_t_deg,energy*1000.0,r_core, evstHist::get_Weight_ETeV(energy));
    }
    if(n_pe>49){
      evH_50npe->Fill(theta_p_t_deg,energy*1000.0);
      evH_50npe->Fill_rcore(theta_p_t_deg,energy*1000.0,r_core);
      //
      evH_50npe_w->Fill(theta_p_t_deg,energy*1000.0, evstHist::get_Weight_ETeV(energy));
      evH_50npe_w->Fill_rcore(theta_p_t_deg,energy*1000.0,r_core, evstHist::get_Weight_ETeV(energy));
    }
    if(n_pe>59){
      evH_60npe->Fill(theta_p_t_deg,energy*1000.0);
      evH_60npe->Fill_rcore(theta_p_t_deg,energy*1000.0,r_core);
      //
      evH_60npe_w->Fill(theta_p_t_deg,energy*1000.0, evstHist::get_Weight_ETeV(energy));
      evH_60npe_w->Fill_rcore(theta_p_t_deg,energy*1000.0,r_core, evstHist::get_Weight_ETeV(energy));
    }
    if(n_pe>69){
      evH_70npe->Fill(theta_p_t_deg,energy*1000.0);
      evH_70npe->Fill_rcore(theta_p_t_deg,energy*1000.0,r_core);
      //
      evH_70npe_w->Fill(theta_p_t_deg,energy*1000.0, evstHist::get_Weight_ETeV(energy));
      evH_70npe_w->Fill_rcore(theta_p_t_deg,energy*1000.0,r_core, evstHist::get_Weight_ETeV(energy));
    }
    if(n_pe>79){
      evH_80npe->Fill(theta_p_t_deg,energy*1000.0);
      evH_80npe->Fill_rcore(theta_p_t_deg,energy*1000.0,r_core);
      //
      evH_80npe_w->Fill(theta_p_t_deg,energy*1000.0, evstHist::get_Weight_ETeV(energy));
      evH_80npe_w->Fill_rcore(theta_p_t_deg,energy*1000.0,r_core, evstHist::get_Weight_ETeV(energy));
    }
    if(n_pe>89){
      evH_90npe->Fill(theta_p_t_deg,energy*1000.0);
      evH_90npe->Fill_rcore(theta_p_t_deg,energy*1000.0,r_core);
      //
      evH_90npe_w->Fill(theta_p_t_deg,energy*1000.0, evstHist::get_Weight_ETeV(energy));
      evH_90npe_w->Fill_rcore(theta_p_t_deg,energy*1000.0,r_core, evstHist::get_Weight_ETeV(energy));
    }
    if(n_pe>99){
      evH_100npe->Fill(theta_p_t_deg,energy*1000.0);
      evH_100npe->Fill_rcore(theta_p_t_deg,energy*1000.0,r_core);
      //
      evH_100npe_w->Fill(theta_p_t_deg,energy*1000.0, evstHist::get_Weight_ETeV(energy));
      evH_100npe_w->Fill_rcore(theta_p_t_deg,energy*1000.0,r_core, evstHist::get_Weight_ETeV(energy));
    }
    //
    n_ev_cuts++;
    if(n_ev_cuts>_anaConf.n_ev_cuts_max && _anaConf.n_ev_cuts_max>0)
      break;
  }
  //
  evH_eff->Divide( evH, evH_simtel, true);
  evH_eff_10npe->Divide( evH_10npe, evH_simtel, true);
  evH_eff_20npe->Divide( evH_20npe, evH_simtel, true);
  evH_eff_30npe->Divide( evH_30npe, evH_simtel, true);
  evH_eff_40npe->Divide( evH_40npe, evH_simtel, true);
  evH_eff_50npe->Divide( evH_50npe, evH_simtel, true);
  evH_eff_60npe->Divide( evH_60npe, evH_simtel, true);
  evH_eff_70npe->Divide( evH_70npe, evH_simtel, true);
  evH_eff_80npe->Divide( evH_80npe, evH_simtel, true);
  evH_eff_90npe->Divide( evH_90npe, evH_simtel, true);
  evH_eff_100npe->Divide( evH_100npe, evH_simtel, true);
  //
  evH_eff_10npe_w->Divide( evH_10npe_w, evH_simtel_weight, true);
  evH_eff_20npe_w->Divide( evH_20npe_w, evH_simtel_weight, true);
  evH_eff_30npe_w->Divide( evH_30npe_w, evH_simtel_weight, true);
  evH_eff_40npe_w->Divide( evH_40npe_w, evH_simtel_weight, true);
  evH_eff_50npe_w->Divide( evH_50npe_w, evH_simtel_weight, true);
  evH_eff_60npe_w->Divide( evH_60npe_w, evH_simtel_weight, true);
  evH_eff_70npe_w->Divide( evH_70npe_w, evH_simtel_weight, true);
  evH_eff_80npe_w->Divide( evH_80npe_w, evH_simtel_weight, true);
  evH_eff_90npe_w->Divide( evH_90npe_w, evH_simtel_weight, true);
  evH_eff_100npe_w->Divide( evH_100npe_w, evH_simtel_weight, true);
  //
  evH_flux->Multiply(evH_tot_flux, evH_eff, true);
  evH_flux_10npe->Multiply(evH_tot_flux, evH_eff_10npe, true);
  evH_flux_20npe->Multiply(evH_tot_flux, evH_eff_20npe, true);
  evH_flux_30npe->Multiply(evH_tot_flux, evH_eff_30npe, true);
  evH_flux_40npe->Multiply(evH_tot_flux, evH_eff_40npe, true);
  evH_flux_50npe->Multiply(evH_tot_flux, evH_eff_50npe, true);
  evH_flux_60npe->Multiply(evH_tot_flux, evH_eff_60npe, true);
  evH_flux_70npe->Multiply(evH_tot_flux, evH_eff_70npe, true);
  evH_flux_80npe->Multiply(evH_tot_flux, evH_eff_80npe, true);
  evH_flux_90npe->Multiply(evH_tot_flux, evH_eff_90npe, true);
  evH_flux_100npe->Multiply(evH_tot_flux, evH_eff_100npe, true);
  //
  evH_flux_10npe_w->Multiply(evH_tot_flux, evH_eff_10npe_w, true);
  evH_flux_20npe_w->Multiply(evH_tot_flux, evH_eff_20npe_w, true);
  evH_flux_30npe_w->Multiply(evH_tot_flux, evH_eff_30npe_w, true);
  evH_flux_40npe_w->Multiply(evH_tot_flux, evH_eff_40npe_w, true);
  evH_flux_50npe_w->Multiply(evH_tot_flux, evH_eff_50npe_w, true);
  evH_flux_60npe_w->Multiply(evH_tot_flux, evH_eff_60npe_w, true);
  evH_flux_70npe_w->Multiply(evH_tot_flux, evH_eff_70npe_w, true);
  evH_flux_80npe_w->Multiply(evH_tot_flux, evH_eff_80npe_w, true);
  evH_flux_90npe_w->Multiply(evH_tot_flux, evH_eff_90npe_w, true);
  evH_flux_100npe_w->Multiply(evH_tot_flux, evH_eff_100npe_w, true);  
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
  //
  for(unsigned int ii = 0;ii<evH_eff->get_v_r().size();ii++)
    evH_eff->get_v_r().at(ii)->Write();
  for(unsigned int ii = 0;ii<evH_eff_50npe->get_v_r().size();ii++)
    evH_eff_50npe->get_v_r().at(ii)->Write();
  for(unsigned int ii = 0;ii<evH_eff_50npe->get_v_r().size();ii++)
    evH_eff_100npe->get_v_r().at(ii)->Write();
  for(unsigned int ii = 0;ii<evH_flux->get_v_r().size();ii++)
    evH_flux->get_v_r().at(ii)->Write();
  for(unsigned int ii = 0;ii<evH_20npe->get_v_r().size();ii++)
    evH_20npe->get_v_r().at(ii)->Write();
  //
  //
  TDirectory *cd_evH_nhe = rootFile->mkdir("evH_nhe");
  cd_evH_nhe->cd();
  for(unsigned int ii = 0;ii<evH_nhe->get_h1_arb_v().size();ii++)
    evH_nhe->get_h1_arb_v().at(ii)->Write();
  rootFile->cd();
  //
  TDirectory *cd_evH_nhe_w = rootFile->mkdir("evH_nhe_w");
  cd_evH_nhe_w->cd();
  for(unsigned int ii = 0;ii<evH_nhe_w->get_h1_arb_v().size();ii++)
    evH_nhe_w->get_h1_arb_v().at(ii)->Write();
  rootFile->cd();  
  //
  //
  //
  evH->Write();
  evH_simtel->Write();
  evH_tot_flux->Write();
  evH_eff->Write();
  evH_eff_10npe->Write();
  evH_eff_20npe->Write();
  evH_eff_30npe->Write();
  evH_eff_40npe->Write();
  evH_eff_50npe->Write();
  evH_eff_60npe->Write();
  evH_eff_70npe->Write();
  evH_eff_80npe->Write();
  evH_eff_90npe->Write();
  evH_eff_100npe->Write();
  evH_flux->Write();
  evH_flux_10npe->Write();
  evH_flux_20npe->Write();
  evH_flux_30npe->Write();
  evH_flux_40npe->Write();
  evH_flux_50npe->Write();
  evH_flux_60npe->Write();
  evH_flux_70npe->Write();
  evH_flux_80npe->Write();
  evH_flux_90npe->Write();
  evH_flux_100npe->Write();
  //
  //
  evH_eff_10npe_w->Write();
  evH_eff_20npe_w->Write();
  evH_eff_30npe_w->Write();
  evH_eff_40npe_w->Write();
  evH_eff_50npe_w->Write();
  evH_eff_60npe_w->Write();
  evH_eff_70npe_w->Write();
  evH_eff_80npe_w->Write();
  evH_eff_90npe_w->Write();
  evH_eff_100npe_w->Write();
  //
  evH_flux_10npe_w->Write();
  evH_flux_20npe_w->Write();
  evH_flux_30npe_w->Write();
  evH_flux_40npe_w->Write();
  evH_flux_50npe_w->Write();
  evH_flux_60npe_w->Write();
  evH_flux_70npe_w->Write();
  evH_flux_80npe_w->Write();
  evH_flux_90npe_w->Write();
  evH_flux_100npe_w->Write();
  //
  //
  evH->Draw_hist("")->Write();
  evH_simtel->Draw_hist("")->Write();
  evH_tot_flux->Draw_hist("")->Write();
  evH_eff->Draw_hist("")->Write();
  evH_eff_10npe->Draw_hist("")->Write();
  evH_eff_20npe->Draw_hist("")->Write();
  evH_eff_30npe->Draw_hist("")->Write();
  evH_eff_40npe->Draw_hist("")->Write();
  evH_eff_50npe->Draw_hist("")->Write();
  evH_eff_60npe->Draw_hist("")->Write();
  evH_eff_70npe->Draw_hist("")->Write();
  evH_eff_80npe->Draw_hist("")->Write();
  evH_eff_90npe->Draw_hist("")->Write();
  evH_eff_100npe->Draw_hist("")->Write();
  evH_flux->Draw_hist("")->Write();
  evH_flux_10npe->Draw_hist("")->Write();
  evH_flux_20npe->Draw_hist("")->Write();
  evH_flux_30npe->Draw_hist("")->Write();
  evH_flux_40npe->Draw_hist("")->Write();
  evH_flux_50npe->Draw_hist("")->Write();
  evH_flux_60npe->Draw_hist("")->Write();
  evH_flux_70npe->Draw_hist("")->Write();
  evH_flux_80npe->Draw_hist("")->Write();
  evH_flux_90npe->Draw_hist("")->Write();
  evH_flux_100npe->Draw_hist("")->Write();
  //
  //
  evH_eff_10npe_w->Draw_hist("")->Write();
  evH_eff_20npe_w->Draw_hist("")->Write();
  evH_eff_30npe_w->Draw_hist("")->Write();
  evH_eff_40npe_w->Draw_hist("")->Write();
  evH_eff_50npe_w->Draw_hist("")->Write();
  evH_eff_60npe_w->Draw_hist("")->Write();
  evH_eff_70npe_w->Draw_hist("")->Write();
  evH_eff_80npe_w->Draw_hist("")->Write();
  evH_eff_90npe_w->Draw_hist("")->Write();
  evH_eff_100npe_w->Draw_hist("")->Write();
  //
  evH_flux_10npe_w->Draw_hist("")->Write();
  evH_flux_20npe_w->Draw_hist("")->Write();
  evH_flux_30npe_w->Draw_hist("")->Write();
  evH_flux_40npe_w->Draw_hist("")->Write();
  evH_flux_50npe_w->Draw_hist("")->Write();
  evH_flux_60npe_w->Draw_hist("")->Write();
  evH_flux_70npe_w->Draw_hist("")->Write();
  evH_flux_80npe_w->Draw_hist("")->Write();
  evH_flux_90npe_w->Draw_hist("")->Write();
  evH_flux_100npe_w->Draw_hist("")->Write();
  //
  //
  evH_10npe->Write();
  evH_20npe->Write();
  evH_30npe->Write();
  evH_40npe->Write();
  evH_50npe->Write();
  evH_60npe->Write();
  evH_70npe->Write();
  evH_80npe->Write();
  evH_90npe->Write();
  evH_100npe->Write();
  //
  //
  evH_10npe_w->Write();
  evH_20npe_w->Write();
  evH_30npe_w->Write();
  evH_40npe_w->Write();
  evH_50npe_w->Write();
  evH_60npe_w->Write();
  evH_70npe_w->Write();
  evH_80npe_w->Write();
  evH_90npe_w->Write();
  evH_100npe_w->Write();
  //
  //
  evH_10npe->Draw_hist("")->Write();
  evH_20npe->Draw_hist("")->Write();
  evH_30npe->Draw_hist("")->Write();
  evH_40npe->Draw_hist("")->Write();
  evH_50npe->Draw_hist("")->Write();
  evH_60npe->Draw_hist("")->Write();
  evH_70npe->Draw_hist("")->Write();
  evH_80npe->Draw_hist("")->Write();
  evH_90npe->Draw_hist("")->Write();
  evH_100npe->Draw_hist("")->Write();
  //
  //
  evH_10npe_w->Draw_hist("")->Write();
  evH_20npe_w->Draw_hist("")->Write();
  evH_30npe_w->Draw_hist("")->Write();
  evH_40npe_w->Draw_hist("")->Write();
  evH_50npe_w->Draw_hist("")->Write();
  evH_60npe_w->Draw_hist("")->Write();
  evH_70npe_w->Draw_hist("")->Write();
  evH_80npe_w->Draw_hist("")->Write();
  evH_90npe_w->Draw_hist("")->Write();
  evH_100npe_w->Draw_hist("")->Write();
  //
  h1_n_pe->Write();
  //
  //
  rootFile->Close();
  //
  cout<<"n_ev_cuts = "<<n_ev_cuts<<endl;  
}

void anaFast::Loop(TString histOut){
  //
  cout<<"anaFast::Loop"<<endl;
  _anaConf.printInfo();
  //assert(0);
  //
  const unsigned int nn_fadc_point = 75;
  const unsigned int nn_PMT_channels = 7987;
  //
  Int_t fadc_sum_offset = 15;
  Int_t fadc_MHz = 1024;
  Int_t fadc_offset = 300;
  Float_t fadc_sample_in_ns = 1000.0/fadc_MHz;
  Float_t time_offset = fadc_sum_offset*fadc_sample_in_ns;
  Float_t NGB_rate_in_MHz = 0.0;
  Float_t fadc_electronic_noise_RMS = 0.0;
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
  Double_t val_Emin = 1.0;    // GeV
  Double_t val_Emax = 100000; // GeV
  Int_t val_N_bins_E = 25;
  //
  Double_t val_Thetamin = 0.0;  //deg
  Double_t val_Thetamax = 10.0; //deg
  Int_t val_N_bins_t = 10;
  //
  evstHist *evH_E_all = new evstHist("evH_E_all","evH_E_all",
				     val_Emin, val_Emax, val_N_bins_E,
				     val_Thetamin, val_Thetamax, val_N_bins_t);
  //
  evstHist *evH_E_cut = new evstHist("evH_E_cut","evH_E_cut",
				     val_Emin, val_Emax, val_N_bins_E,
				     val_Thetamin, val_Thetamax, val_N_bins_t);
  evstHist::PrintBinsInfo(evH_E_all->get_E_hist());
  //
  
  //
  evstHist *evH_all = new evstHist("evH_all","evH_all",
				   val_Emin, val_Emax, val_N_bins_E,
				   val_Thetamin, val_Thetamax, val_N_bins_t);
  //
  evstHist *evH_cut = new evstHist("evH_cut","evH_cut",
				   val_Emin, val_Emax, val_N_bins_E,
				   val_Thetamin, val_Thetamax, val_N_bins_t);
  //
  evstHist *evH_eff = new evstHist("evH_eff","evH_eff",
				   val_Emin, val_Emax, val_N_bins_E,
				   val_Thetamin, val_Thetamax, val_N_bins_t);
  //
  evstHist *evH_eff_one = new evstHist("evH_eff_one","evH_eff_one",
				       val_Emin, val_Emax, val_N_bins_E,
				       val_Thetamin, val_Thetamax, val_N_bins_t);
  //
  evstHist *evH_pe_cut = new evstHist("evH_pe_cut","evH_pe_cut",
				      val_Emin, val_Emax, val_N_bins_E,
				      val_Thetamin, val_Thetamax, val_N_bins_t);
  //
  evstHist *evH_flux_tot = new evstHist("evH_flux_tot","evH_flux_tot",
					val_Emin, val_Emax, val_N_bins_E,
					val_Thetamin, val_Thetamax, val_N_bins_t);
  evH_flux_tot->LoadBinContent("../cosmique_gamma_hadron_generator/flux.dat");
  //
  //
  evstHist *evH_all_simtel = new evstHist("evH_all_simtel","evH_all_simtel",
					  val_Emin, val_Emax, val_N_bins_E,
					  val_Thetamin, val_Thetamax, val_N_bins_t);
  evH_all_simtel->LoadBinContent("../cosmique_gamma_hadron_generator/flux_simtel_all.dat");
  //evH_all_simtel->LoadBinContent("../cosmique_gamma_hadron_generator/flux_simtel_all_gamma_diff.dat");
  //
  evstHist *evH_flux_one_final = new evstHist("evH_flux_one_final","evH_flux_one_final",
					val_Emin, val_Emax, val_N_bins_E,
					val_Thetamin, val_Thetamax, val_N_bins_t);
  //
  evstHist *evH_flux_final = new evstHist("evH_flux_final","evH_flux_final",
					val_Emin, val_Emax, val_N_bins_E,
					val_Thetamin, val_Thetamax, val_N_bins_t);
  Double_t theta_p_t;
  Double_t theta_p_t_deg;
  //
  //
  sipmCameraHist *sipm_cam_trg = new sipmCameraHist("sipm_cam_trg","sipm_cam_trg","pixel_mapping.csv",0);
  triggerSim *trg_sim = new triggerSim(sipm_cam_trg);
  //
  sipmCameraHist *sipm_cam = new sipmCameraHist("sipm_cam","sipm_cam","pixel_mapping.csv",0);
  //
  TH1D *h1_energy = new TH1D("h1_energy","h1_energy",110000,0.0,110.0);
  //
  TH1D *h1_nphotons = new TH1D("h1_nphotons","h1_nphotons",10000,0.0,1000000);
  TH1D *h1_n_pe = new TH1D("h1_n_pe","h1_n_pe",1001,-0.5,10000.5);
  TH1D *h1_n_pe_notcuts = new TH1D("h1_n_pe_notcuts","h1_n_pe_notcuts",1001,-0.5,1000.5);
  TH1D *h1_n_pixels = new TH1D("h1_n_pixels","h1_n_pixels",1000,0.0,10000);
  //
  TH1D *h1_azimuth = new TH1D("h1_azimuth","azimuth",400,-4,4);
  TH1D *h1_altitude = new TH1D("h1_altitude","h1_altitude",400,-4,4);
  //
  TH1D *h1_theta_p_t_deg = new TH1D("theta_p_t_deg","theta_p_t_deg",4000,-370,370);
  //
  TH1D *h1_azimuth_deg = new TH1D("h1_azimuth_deg","azimuth_deg",400,0.0,360);
  TH1D *h1_altitude_deg = new TH1D("h1_altitude_deg","h1_altitude_deg",400,0.0,180);
  //
  TH1D *h1_azimuth_notcuts = new TH1D("h1_azimuth_notcuts","azimuth notcuts",400,-4,4);
  TH1D *h1_altitude_notcuts = new TH1D("h1_altitude_notcuts","altitude notcuts",400,-4,4);
  //
  TH1D *h1_h_first_int = new TH1D("h1_h_first_int","h1_h_first_int",4000,0.0,100000);
  TH1D *h1_xmax = new TH1D("h1_xmax","h1_xmax",400,0.0,1000);
  TH1D *h1_hmax = new TH1D("h1_hmax","h1_hmax",400,0.0,100000);
  TH1D *h1_hmax_notcuts = new TH1D("h1_hmax_notcuts","h1_hmax_notcuts",400,0.0,100000);
  TH1D *h1_emax = new TH1D("h1_emax","h1_emax",400,0.0,1000);
  TH1D *h1_cmax = new TH1D("h1_cmax","h1_cmax",400,0.0,1000);
  //
  TH2D *h2_hmax_vs_h_first_int = new TH2D("h2_hmax_vs_h_first_int","h2_hmax_vs_h_first_int",400,0.0,20000,400,0.0,50000);
  TH2D *h2_n_pe_vs_hmax = new TH2D("h2_n_pe_vs_hmax","h2_n_pe_vs_hmax",400,0.0,10000,400,0.0,20000);
  //
  TH1D *h1_xcore = new TH1D("h1_xcore","h1_xcore",800,-2000,2000);
  TH1D *h1_ycore = new TH1D("h1_ycore","h1_ycore",800,-2000,2000);
  TH1D *h1_xcore_km = new TH1D("h1_xcore_km","h1_xcore_km",200,-1.8,1.8);
  TH1D *h1_ycore_km = new TH1D("h1_ycore_km","h1_ycore_km",200,-1.8,1.8);
  TH1D *h1_r_core = new TH1D("h1_r_core","h1_r_core",400, 0,_anaConf.ellipse_A*1.2);
  TH1D *h1_theta_core = new TH1D("h1_theta_core","h1_theta_core",400,-0.1,2*TMath::Pi()+0.1);
  //
  TH2D *h2_ycore_vs_xcore = new TH2D("h2_ycore_vs_xcore","h2_ycore_vs_xcore",400,-2000,2000, 400,-2000,2000);  
  //
  TH2D *h2_ycore_vs_xcore_ring = new TH2D("h2_ycore_vs_xcore_ring","h2_ycore_vs_xcore_ring",400,-2000,2000, 400,-2000,2000);  
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
  sipmCameraHistCropped *tmp_cam_hist = new sipmCameraHistCropped("tmp","tmp",sipm_cam,"sipmCameraHistCropped_pix.map");  
  //
  Double_t angleSim;
  //
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  if(_anaConf.nentries_max>0)
    nentries = _anaConf.nentries_max;
  cout<<"nentries = "<<nentries<<endl;
  //
  //
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    //for (Long64_t jentry=0; jentry<100000;jentry++) {
    //if(jentry%100 == 0)
    if(jentry%_anaConf.jentry_modulo == 0)
      cout<<jentry<<endl;
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    //if (ientry > 10) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    //
    //
    getCore_rel_R_theta( x0_LST01, y0_LST01, xcore, ycore, r_core, theta_core);
    _r_core = r_core;
    _theta_core = theta_core;
    //
    //
    h1_azimuth_notcuts->Fill(azimuth);
    h1_altitude_notcuts->Fill(altitude);
    //
    azimuth_deg = azimuth*180.0/TMath::Pi();
    altitude_deg = altitude*180.0/TMath::Pi();
    //
    h1_n_pe_notcuts->Fill(n_pe);
    h1_hmax_notcuts->Fill(hmax);
    //
    theta_p_t = get_theta_p_t();
    theta_p_t_deg = theta_p_t*180/TMath::Pi();
    //
    h1_theta_p_t_deg->Fill(theta_p_t_deg);
    if(theta_p_t_deg<10.0){
      evH_all->Fill(theta_p_t_deg,energy*1000.0);
      evH_all->Fill_rcore(theta_p_t_deg,energy*1000.0,_r_core);
      //
      evH_E_all->get_E_hist()->Fill(energy*1000.0);
      //
      if(cuts(theta_p_t_deg)){
	n_ev_cuts++;
	//
	if(n_ev_cuts>_anaConf.n_ev_cuts_max && _anaConf.n_ev_cuts_max>0)
	  break;
	//
	//
	if(_anaConf.if_cout_jentry){
	  //cout<<"jentry = "<<jentry<<"   "<<n_pe<<energy<<endl;
	  //./runana 4 $inRootFile $outHistF  588 gamma  #   68
	  //cout<<"./runana 4 $inRootFile $outHistF "<<jentry<<" gamma # n_pe "<<n_pe<<"  energy "<<energy*1000.0<<" GeV"<<" azi. "<<azimuth_deg<<" alt. "<<altitude_deg<<endl;
	  cout<<"./runana 4 $inRootFile $outHistF "<<jentry<<" proton # n_pe "<<n_pe<<"  energy "<<energy*1000.0<<" GeV"<<" azi. "<<azimuth_deg<<" alt. "<<altitude_deg<<endl;
	}
	//
	evH_E_cut->get_E_hist()->Fill(energy*1000.0);
	//
	h1_energy->Fill(energy);
	h1_nphotons->Fill(nphotons);
	h1_n_pe->Fill(n_pe);
	h1_n_pixels->Fill(n_pixels);
	//
	h1_azimuth->Fill(azimuth);
	h1_altitude->Fill(altitude);
	//
	h1_azimuth_deg->Fill(azimuth_deg);
	h1_altitude_deg->Fill(altitude_deg);
	//
	h1_h_first_int->Fill(h_first_int);
	//
	h2_hmax_vs_h_first_int->Fill(hmax,h_first_int);
	h2_n_pe_vs_hmax->Fill(n_pe,hmax);
	//
	h1_xmax->Fill(xmax);
	h1_hmax->Fill(hmax);
	h1_emax->Fill(emax);
	h1_cmax->Fill(cmax);
	//
	h1_xcore->Fill(xcore);
	h1_ycore->Fill(ycore);
	//
	h1_xcore_km->Fill(xcore/1000.0);
	h1_ycore_km->Fill(ycore/1000.0);
	//
	h2_ycore_vs_xcore->Fill(xcore,ycore);
	angleSim = rnd->Uniform(0.0,2.0*TMath::Pi());
	h2_ycore_vs_xcore_ring->Fill(_anaConf.ellipse_A*TMath::Cos(angleSim),_anaConf.ellipse_B*TMath::Sin(angleSim));
	//
	//h2_ycore_vs_xcore->Fill(xcore,ycore);
	//
	h1_r_core->Fill(r_core);
	h1_theta_core->Fill(theta_core);
	//
	evH_cut->Fill(theta_p_t_deg,energy*1000.0);
	evH_cut->Fill_rcore(theta_p_t_deg,energy*1000.0,_r_core);
	evH_pe_cut->Fill(theta_p_t_deg,energy*1000.0,n_pe);
      }
    }
  }
  //
  // Do the actual analysis
  //
  evH_eff->Divide(evH_cut,evH_all_simtel);
  evH_eff_one->Divide(evH_all,evH_all_simtel);
  evH_flux_one_final->Multiply(evH_eff_one,evH_flux_tot);
  evH_flux_final->Multiply(evH_eff,evH_flux_tot);
  //evH_flux_final->SetMaximum(1.0e+4);
  //evH_flux_final->SetMinimum(1.0e-4);
  //
  //  
  //proton
  evH_flux_tot->SetMinimum(1.0e-1);
  evH_flux_tot->SetMaximum(1.0e+9);
  evH_all->SetMinimum(1.0e+0);
  evH_all->SetMaximum(1.0e+7);
  evH_cut->SetMinimum(1.0e+0);
  evH_cut->SetMaximum(1.0e+7);
  evH_eff->SetMinimum(1.0e-8);
  evH_eff->SetMaximum(1.0e+0);
  evH_eff_one->SetMinimum(1.0e-3);
  evH_eff_one->SetMaximum(1.0e+0);
  evH_flux_final->SetMinimum(1.0e-2);
  evH_flux_final->SetMaximum(5.0*1.0e+2);
  evH_flux_one_final->SetMinimum(1.0e-1);
  evH_flux_one_final->SetMaximum(2.5*1.0e+4);
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
  //
  for(unsigned int ii = 0;ii<evH_all->get_v_r().size();ii++)
    evH_all->get_v_r().at(ii)->Write();
  for(unsigned int ii = 0;ii<evH_cut->get_v_r().size();ii++)
    evH_cut->get_v_r().at(ii)->Write();
  //
  //
  h1_energy->Write();
  h1_nphotons->Write();
  h1_n_pe->Write();
  h1_n_pe_notcuts->Write();
  h1_n_pixels->Write();
  //
  h1_azimuth->Write();
  h1_altitude->Write();
  //
  h1_azimuth_deg->Write();
  h1_altitude_deg->Write();
  //
  h1_azimuth_notcuts->Write();
  h1_altitude_notcuts->Write();
  //
  h1_h_first_int->Write();
  h1_xmax->Write();
  h1_hmax->Write();
  h1_hmax_notcuts->Write();
  h1_emax->Write();
  h1_cmax->Write();
  //
  h2_hmax_vs_h_first_int->Write();
  h2_n_pe_vs_hmax->Write();
  //
  h1_xcore->Write();
  h1_ycore->Write();
  h1_r_core->Write();
  h1_theta_core->Write();
  //
  h2_ycore_vs_xcore->Write();
  h2_ycore_vs_xcore_ring->Write();
  //
  h1_xcore_km->Write();
  h1_ycore_km->Write();
  //
  h1_theta_p_t_deg->Write();
  //
  evH_all->Write();
  evH_cut->Write();
  evH_pe_cut->Write();
  evH_eff->Write();
  evH_flux_tot->Write();
  evH_flux_final->Write();
  evH_flux_one_final->Write();
  evH_all_simtel->Write();
  evH_eff_one->Write();
  //
  evH_all->_title="All events";
  evH_all->Draw_hist("c1_evH_all.pdf")->Write();
  evH_cut->_title="Events with more than 100 p.e.";
  evH_cut->Draw_hist("c1_evH_cut.pdf")->Write();
  evH_pe_cut->Draw_hist("")->Write();
  evH_eff->_title="Efficiency";
  evH_eff->Draw_hist("c1_evH_eff.pdf")->Write();
  evH_flux_tot->_title="Total flux, Hz";
  evH_flux_tot->Draw_hist("c1_evH_flux_tot.pdf")->Write();
  evH_flux_final->_title="Flux (with more then 100 p.e.), Hz";
  evH_flux_final->Draw_hist("c1_evH_flux_final.pdf")->Write();
  evH_flux_one_final->_title="Flux (with more then 1 p.e.), Hz";
  evH_flux_one_final->Draw_hist("c1_evH_flux_one_final.pdf")->Write();
  evH_all_simtel->Draw_hist("c1_evH_all_simtel.pdf")->Write();
  evH_eff_one->Draw_hist("c1_evH_eff_one.pdf")->Write();
  //
  evH_E_all->get_E_hist()->Write();
  evH_E_cut->get_E_hist()->Write();
  //
  rootFile->Close();
  //
  cout<<"n_ev_cuts = "<<n_ev_cuts<<endl;
}

bool anaFast::cuts(Double_t theta_p_t_deg){
  //
  if(_anaConf.disable_all_cuts)
    return true;
  //
  if(_anaConf.cuts_set_to_false)
    return false;
  Double_t azimuth_min = (180.0 - 0.2)/180.0*TMath::Pi();
  Double_t azimuth_max = (180.0 + 0.2)/180.0*TMath::Pi();
  Double_t altitude_min = (90.0 - 20.0 - 1.0 - 0.2)/180.0*TMath::Pi();
  Double_t altitude_max = (90.0 - 20.0 - 1.0 + 0.2)/180.0*TMath::Pi();
  //
  //Double_t x0_LST01 = -70.93;
  //Double_t y0_LST01 = -52.07;
  //Double_t r = TMath::Sqrt((x0_LST01 - xcore)*(x0_LST01 - xcore) + (y0_LST01 - ycore)*(y0_LST01 - ycore));
  //
  if(energy>=75.0)
    if(n_pe>=12000)
      if(n_pe<13000)
	return true;
  //if(n_pe>=0)
  //return true;
  //if(n_pe>=500){
  //if(azimuth>=azimuth_min && azimuth<=azimuth_max){
  //if(altitude>=altitude_min && altitude<=altitude_max){
  //if(r<=100)
  //return true;
  //}
  //}
  //}
  //
  Double_t x0_LST01 = -70.93;
  Double_t y0_LST01 = -52.07;
  Double_t r = TMath::Sqrt((x0_LST01 - xcore)*(x0_LST01 - xcore) + (y0_LST01 - ycore)*(y0_LST01 - ycore));
  //if(energy>=0.010 && energy<0.02)
  //if(n_pe>=40)
  //if(r<=300)
  //return true;
  //proton
  //if(theta_p_t_deg<3.0)
  //if(energy>=3.98107)
  //if(r<=150)
  //return true;     
  //gamma diff, gamma, electron
  //if(theta_p_t_deg<2.0)
  //if(energy>=0.20)
  //if(r<=70)
  //return true;
  //
  //TVector2 v_core(xcore,ycore);
  //Double_t r_core = TMath::Sqrt(xcore*xcore + ycore*ycore);
  //Double_t r_ellipse = TMath::Sqrt(_anaConf.ellipse_A*TMath::Cos(v_core.Phi())*_anaConf.ellipse_A*TMath::Cos(v_core.Phi()) +
  //			  _anaConf.ellipse_B*TMath::Sin(v_core.Phi())*_anaConf.ellipse_B*TMath::Sin(v_core.Phi()));
  //if(r_core<=r_ellipse)
  //return true;
  return false;
}

Double_t anaFast::get_theta_p_t(){
  TVector3 v_prot;
  v_prot.SetMagThetaPhi(1.0,TMath::Pi()/2.0-altitude,TMath::Pi() - azimuth);
  TVector3 v_prot_inv(v_prot.x(),v_prot.y(),v_prot.z());
  return TMath::ACos(v_prot_inv.Dot(_v_det)/v_prot_inv.Mag()/_v_det.Mag());
}
