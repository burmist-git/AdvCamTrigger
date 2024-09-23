//my
#include "anastereo.hh"
#include "anabase.hh"
#include "sipmCameraHist.hh"
#include "evstHist.hh"
#include "anaTrgB.hh"
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

anastereo::anastereo(TString fileList) : anabasestereo(fileList),
					 _name_ana_conf_file("./anaStereo_p.conf")
{
}

anastereo::anastereo(TString file, Int_t key) : anabasestereo(file, key),
						_name_ana_conf_file("./anaStereo_p.conf")
{
}

void anastereo::Loop(TString histOut){
  //
  _anaConf.readFromFile(_name_ana_conf_file);
  _anaConf.printInfo();
  //
  TH1D *h1_n_tot_event   = new TH1D("h1_n_tot_event",   "h1_n_tot_event",   3, -1.5, 1.5);
  TH1D *h1_n_tot_event_w = new TH1D("h1_n_tot_event_w", "h1_n_tot_event_w", 3, -1.5, 1.5);
  //
  TH1D *h1_n_pe_LST_all       = new TH1D( "h1_n_pe_LST_all",    "h1_n_pe_LST_all",    3, -1.5, 1.5);
  TH1D *h1_n_pe_LST_all_LST01 = new TH1D( "h1_n_pe_LST_all_LST01",    "h1_n_pe_LST_all_LST01",    3, -1.5, 1.5);
  TH1D *h1_n_pe_LST_all_LST02 = new TH1D( "h1_n_pe_LST_all_LST02",    "h1_n_pe_LST_all_LST02",    3, -1.5, 1.5);
  TH1D *h1_n_pe_LST_all_LST03 = new TH1D( "h1_n_pe_LST_all_LST03",    "h1_n_pe_LST_all_LST03",    3, -1.5, 1.5);
  TH1D *h1_n_pe_LST_all_LST04 = new TH1D( "h1_n_pe_LST_all_LST04",    "h1_n_pe_LST_all_LST04",    3, -1.5, 1.5);
  TH1D *h1_n_pe_LST_single    = new TH1D( "h1_n_pe_LST_single", "h1_n_pe_LST_single", 6, -1.5, 4.5);
  TH1D *h1_n_pe_LST_double    = new TH1D( "h1_n_pe_LST_double", "h1_n_pe_LST_double", 8, -1.5, 6.5);
  TH1D *h1_n_pe_LST_triple    = new TH1D( "h1_n_pe_LST_triple", "h1_n_pe_LST_triple", 6, -1.5, 4.5);
  TH1D *h1_n_pe_LST_four      = new TH1D( "h1_n_pe_LST_four",   "h1_n_pe_LST_four",   3, -1.5, 1.5);
  //
  TH2D *h2_n_pe_LST2_vs_n_pe_LST1 = new TH2D("h2_n_pe_LST2_vs_n_pe_LST1","h2_n_pe_LST2_vs_n_pe_LST1",400,0.0,1000,400,0.0,1000);
  //
  TH1D *h1_xcore = new TH1D("h1_xcore","h1_xcore",1000,-2000.0,2000);
  TH1D *h1_ycore = new TH1D("h1_ycore","h1_ycore",1000,-2000.0,2000);
  TH1D *h1_xcore_LST1 = new TH1D("h1_xcore_LST1","h1_xcore_LST1",1000,-2000.0,2000);
  TH1D *h1_ycore_LST1 = new TH1D("h1_ycore_LST1","h1_ycore_LST1",1000,-2000.0,2000);
  TH1D *h1_xcore_LST2 = new TH1D("h1_xcore_LST2","h1_xcore_LST2",1000,-2000.0,2000);
  TH1D *h1_ycore_LST2 = new TH1D("h1_ycore_LST2","h1_ycore_LST2",1000,-2000.0,2000);
  TH1D *h1_xcore_LST3 = new TH1D("h1_xcore_LST3","h1_xcore_LST3",1000,-2000.0,2000);
  TH1D *h1_ycore_LST3 = new TH1D("h1_ycore_LST3","h1_ycore_LST3",1000,-2000.0,2000);
  TH1D *h1_xcore_LST4 = new TH1D("h1_xcore_LST4","h1_xcore_LST4",1000,-2000.0,2000);
  TH1D *h1_ycore_LST4 = new TH1D("h1_ycore_LST4","h1_ycore_LST4",1000,-2000.0,2000);
  //
  TH2D *h2_ycore_vs_xcore_all = new TH2D("h2_ycore_vs_xcore_all","h2_ycore_vs_xcore_all", 50,-2000.0,2000, 50,-2000.0,2000);
  TH2D *h2_ycore_vs_xcore_single_LST1 = new TH2D("h2_ycore_vs_xcore_single_LST1","h2_ycore_vs_xcore_single_LST1",50,-2000.0,2000, 50,-2000.0,2000);
  TH2D *h2_ycore_vs_xcore_double_LST1LST4 = new TH2D("h2_ycore_vs_xcore_double_LST1LST4","h2_ycore_vs_xcore_double_LST1LST4", 50,-2000.0,2000, 50,-2000.0,2000);
  TH2D *h2_ycore_vs_xcore_triple_LST1LST2LST3 = new TH2D("h2_ycore_vs_xcore_triple_LST1LST2LST3","h2_ycore_vs_xcore_triple_LST1LST2LST3", 50,-2000.0,2000, 50,-2000.0,2000);
  TH2D *h2_ycore_vs_xcore_four_LST = new TH2D("h2_ycore_vs_xcore_four_LST","h2_ycore_vs_xcore_four_LST", 50,-2000.0,2000, 50,-2000.0,2000);
  //
  TH2D *h2_ycore_vs_xcore_single_LST1_norm = new TH2D("h2_ycore_vs_xcore_single_LST1_norm","h2_ycore_vs_xcore_single_LST1_norm",50,-2000.0,2000, 50,-2000.0,2000);
  TH2D *h2_ycore_vs_xcore_double_LST1LST4_norm = new TH2D("h2_ycore_vs_xcore_double_LST1LST4_norm","h2_ycore_vs_xcore_double_LST1LST4_norm", 50,-2000.0,2000, 50,-2000.0,2000);
  TH2D *h2_ycore_vs_xcore_triple_LST1LST2LST3_norm = new TH2D("h2_ycore_vs_xcore_triple_LST1LST2LST3_norm","h2_ycore_vs_xcore_triple_LST1LST2LST3_norm", 50,-2000.0,2000, 50,-2000.0,2000);
  TH2D *h2_ycore_vs_xcore_four_LST_norm = new TH2D("h2_ycore_vs_xcore_four_LST_norm","h2_ycore_vs_xcore_four_LST_norm", 50,-2000.0,2000, 50,-2000.0,2000);
  //
  TH1D *h1_ev_time_LST1 = new TH1D("h1_ev_time_LST1","h1_ev_time_LST1",1000,-2000.0,2000);
  TH1D *h1_ev_time_LST2 = new TH1D("h1_ev_time_LST2","h1_ev_time_LST2",1000,-2000.0,2000);
  TH1D *h1_ev_time_LST3 = new TH1D("h1_ev_time_LST3","h1_ev_time_LST3",1000,-2000.0,2000);
  TH1D *h1_ev_time_LST4 = new TH1D("h1_ev_time_LST4","h1_ev_time_LST4",1000,-2000.0,2000);
  //
  TH1D *h1_dtime_LST1_m_LST2 = new TH1D("h1_dtime_LST1_m_LST2","h1_dtime_LST1_m_LST2",200,-400.0,400);
  TH1D *h1_dtime_LST1_m_LST3 = new TH1D("h1_dtime_LST1_m_LST3","h1_dtime_LST1_m_LST3",200,-400.0,400);
  TH1D *h1_dtime_LST1_m_LST4 = new TH1D("h1_dtime_LST1_m_LST4","h1_dtime_LST1_m_LST4",200,-400.0,400);
  TH1D *h1_dtime_LST2_m_LST3 = new TH1D("h1_dtime_LST2_m_LST3","h1_dtime_LST2_m_LST3",200,-400.0,400);
  TH1D *h1_dtime_LST2_m_LST4 = new TH1D("h1_dtime_LST2_m_LST4","h1_dtime_LST2_m_LST4",200,-400.0,400);
  TH1D *h1_dtime_LST3_m_LST4 = new TH1D("h1_dtime_LST3_m_LST4","h1_dtime_LST3_m_LST4",200,-400.0,400);
  //
  TH1D *h1_dtime_LST1_m_LST2_norm = new TH1D("h1_dtime_LST1_m_LST2_norm","h1_dtime_LST1_m_LST2_norm",200,-400.0,400);
  TH1D *h1_dtime_LST1_m_LST3_norm = new TH1D("h1_dtime_LST1_m_LST3_norm","h1_dtime_LST1_m_LST3_norm",200,-400.0,400);
  TH1D *h1_dtime_LST1_m_LST4_norm = new TH1D("h1_dtime_LST1_m_LST4_norm","h1_dtime_LST1_m_LST4_norm",200,-400.0,400);
  TH1D *h1_dtime_LST2_m_LST3_norm = new TH1D("h1_dtime_LST2_m_LST3_norm","h1_dtime_LST2_m_LST3_norm",200,-400.0,400);
  TH1D *h1_dtime_LST2_m_LST4_norm = new TH1D("h1_dtime_LST2_m_LST4_norm","h1_dtime_LST2_m_LST4_norm",200,-400.0,400);
  TH1D *h1_dtime_LST3_m_LST4_norm = new TH1D("h1_dtime_LST3_m_LST4_norm","h1_dtime_LST3_m_LST4_norm",200,-400.0,400);
  //
  TH1D *h1_dtime_LST1_LST3_m_exp = new TH1D("h1_dtime_LST1_LST3_m_exp","h1_dtime_LST1_m_LST3_norm",300,0.0,300);
  TH1D *h1_dtime_LST1_LST3_m_exp_full = new TH1D("h1_dtime_LST1_LST3_m_exp_full","h1_dtime_LST1_m_LST3_norm_full",300,-300.0,300);
  TH1D *h1_dtime_event_acceptance = new TH1D("h1_dtime_event_acceptance","h1_dtime_event_acceptance",300,0.0,300);  
  //
  //
  TH1D *h1_xmax = new TH1D("h1_xmax","h1_xmax", 1000, 0.0,  1000);
  TH1D *h1_hmax = new TH1D("h1_hmax","h1_hmax", 1000, 0.0, 25000);
  //
  TH2D *h2_xmax_vs_hmax = new TH2D("h2_xmax_vs_hmax","h2_xmax_vs_hmax", 400, 0.0, 25000, 400, 0.0,  1000);
  //
  Int_t coincidenceID;
  Int_t coincidenceIDtmp;
  //
  sipmCameraHist *sipm_cam = new sipmCameraHist("sipm_cam","sipm_cam","pixel_mapping.csv",0); 
  //
  TH1D *h1_LST1_sipm_cam_x = new TH1D("h1_LST1_sipm_cam_x","sipm cam x", 260, -1.3, 1.3);
  TH1D *h1_LST2_sipm_cam_x = new TH1D("h1_LST2_sipm_cam_x","sipm cam x", 260, -1.3, 1.3);
  TH1D *h1_LST3_sipm_cam_x = new TH1D("h1_LST3_sipm_cam_x","sipm cam x", 260, -1.3, 1.3);
  TH1D *h1_LST4_sipm_cam_x = new TH1D("h1_LST4_sipm_cam_x","sipm cam x", 260, -1.3, 1.3);
  //
  TH1D *h1_LST1_sipm_cam_y = new TH1D("h1_LST1_sipm_cam_y","sipm cam y", 260, -1.3, 1.3);
  TH1D *h1_LST2_sipm_cam_y = new TH1D("h1_LST2_sipm_cam_y","sipm cam y", 260, -1.3, 1.3);
  TH1D *h1_LST3_sipm_cam_y = new TH1D("h1_LST3_sipm_cam_y","sipm cam y", 260, -1.3, 1.3);
  TH1D *h1_LST4_sipm_cam_y = new TH1D("h1_LST4_sipm_cam_y","sipm cam y", 260, -1.3, 1.3);
  //
  TH2D *h2_LST1_sipm_cam_y_vs_x = new TH2D("h2_LST1_sipm_cam_y_vs_x","LST1 sipm cam y vs x", 260, -1.3, 1.3, 260, -1.3, 1.3);
  TH2D *h2_LST2_sipm_cam_y_vs_x = new TH2D("h2_LST2_sipm_cam_y_vs_x","LST2 sipm cam y vs x", 260, -1.3, 1.3, 260, -1.3, 1.3);
  TH2D *h2_LST3_sipm_cam_y_vs_x = new TH2D("h2_LST3_sipm_cam_y_vs_x","LST3 sipm cam y vs x", 260, -1.3, 1.3, 260, -1.3, 1.3);
  TH2D *h2_LST4_sipm_cam_y_vs_x = new TH2D("h2_LST4_sipm_cam_y_vs_x","LST4 sipm cam y vs x", 260, -1.3, 1.3, 260, -1.3, 1.3);
  //
  Double_t x_mean_LST1, y_mean_LST1;
  Double_t x_mean_LST2, y_mean_LST2;
  Double_t x_mean_LST3, y_mean_LST3;
  Double_t x_mean_LST4, y_mean_LST4;
  Double_t dr_LST1_LST2;
  Double_t dr_LST1_LST3;
  Double_t dr_LST1_LST4;
  Double_t dr_LST2_LST3;
  Double_t dr_LST2_LST4;
  Double_t dr_LST3_LST4;
  Double_t dalpha_LST1_LST2;
  Double_t dalpha_LST1_LST3;
  Double_t dalpha_LST1_LST4;
  Double_t dalpha_LST2_LST3;
  Double_t dalpha_LST2_LST4;
  Double_t dalpha_LST3_LST4;
  TH1D *h1_dr_LST1_LST2 = new TH1D ("h1_dr_LST1_LST2","h1_dr_LST1_LST2",1000, -2.0,2.0);
  TH1D *h1_dalpha_LST1_LST2 = new TH1D ("h1_dalpha_LST1_LST2","h1_dalpha_LST1_LST2",1000, -3.0,3.0);
  TH1D *h1_dalpha_LST1_LST3 = new TH1D ("h1_dalpha_LST1_LST3","h1_dalpha_LST1_LST3",1000, -3.0,3.0);
  TH1D *h1_dalpha_LST1_LST4 = new TH1D ("h1_dalpha_LST1_LST4","h1_dalpha_LST1_LST4",1000, -3.0,3.0);
  TH1D *h1_dalpha_LST2_LST3 = new TH1D ("h1_dalpha_LST2_LST3","h1_dalpha_LST2_LST3",1000, -3.0,3.0);
  TH1D *h1_dalpha_LST2_LST4 = new TH1D ("h1_dalpha_LST2_LST4","h1_dalpha_LST2_LST4",1000, -3.0,3.0);
  TH1D *h1_dalpha_LST3_LST4 = new TH1D ("h1_dalpha_LST3_LST4","h1_dalpha_LST3_LST4",1000, -3.0,3.0);
  //
  TH2D *h2_dalpha_vs_dtLST1_LST3 = new TH2D("h2_dalpha_vs_dtLST1_LST3","h2_dalpha_vs_dtLST1_LST3",400,-200.0,200.0,400, 0.0,3.0);
  TH2D *h2_x_vs_dt_LST1_LST3 = new TH2D("h2_x_vs_dt_LST1_LST3","h2_x_vs_dt_LST1_LST3",400,-200,200,400, -1.3,1.3);
  //
  TH2D *h2_y_vs_x_LST1_LST2 = new TH2D("h2_y_vs_x_LST1_LST2","h2_y_vs_x_LST1_LST2",100, -1.3,1.3,100, -1.3,1.3);
  TH2D *h2_y_vs_x_LST1_LST3 = new TH2D("h2_y_vs_x_LST1_LST3","h2_y_vs_x_LST1_LST3",100, -1.3,1.3,100, -1.3,1.3);
  TH2D *h2_y_vs_x_LST1_LST4 = new TH2D("h2_y_vs_x_LST1_LST4","h2_y_vs_x_LST1_LST4",100, -1.3,1.3,100, -1.3,1.3);
  //
  TH2D *h2_dt_vs_x_LST1_LST2 = new TH2D("h2_dt_vs_x_LST1_LST2","h2_dt_vs_x_LST1_LST2",400, -1.3,1.3, 400, -200.0,200.0);  
  //
  TH2D *h2_y_vs_x_LST2_LST3 = new TH2D("h2_y_vs_x_LST2_LST3","h2_y_vs_x_LST2_LST3",100, -1.3,1.3,100, -1.3,1.3);
  TH2D *h2_y_vs_x_LST2_LST4 = new TH2D("h2_y_vs_x_LST2_LST4","h2_y_vs_x_LST2_LST4",100, -1.3,1.3,100, -1.3,1.3);
  TH2D *h2_y_vs_x_LST3_LST4 = new TH2D("h2_y_vs_x_LST3_LST4","h2_y_vs_x_LST3_LST4",100, -1.3,1.3,100, -1.3,1.3);  
  //
  TH2D *h2_y_vs_x_LST1_LST3_norm_dt = new TH2D("h2_y_vs_x_LST1_LST3_norm_dt","h2_y_vs_x_LST1_LST3_norm_dt",100, -1.3,1.3,100, -1.3,1.3);  
  //TH2D *h2_y_vs_x_LST1_LST3_norm_dt_NORM = new TH2D("h2_y_vs_x_LST1_LST3_norm_dt","h2_y_vs_x_LST1_LST3_norm_dt",100, -1.3,1.3,100, -1.3,1.3);  
  //
  Int_t npe_trg_max = 1000;
  TH1D *h1_npe_trg_max = new TH1D("h1_npe_trg_max","npe trg max",npe_trg_max+1,-0.5,npe_trg_max+0.5);
  TH1D *h1_npe_trg_max_w = new TH1D("h1_npe_trg_max_w","npe trg max w",npe_trg_max+1,-0.5,npe_trg_max+0.5);
  TH1D *h1_npe_trg_max_OR_w = new TH1D("h1_npe_trg_max_OR_w","npe trg max OR w",npe_trg_max+1,-0.5,npe_trg_max+0.5);
  TH1D *h1_npe_trg_max_LST1_w = new TH1D("h1_npe_trg_max_LST1_w","npe trg max LST1 w",npe_trg_max+1,-0.5,npe_trg_max+0.5);
  TH1D *h1_npe_trg_max_LST1_PMT_w = new TH1D("h1_npe_trg_max_LST1_PMT_w","npe trg max LST1 PMT w",npe_trg_max+1,-0.5,npe_trg_max+0.5);
  Double_t tot_Weight = 0.0;
  TH1D *h1_proton_trg_pate_vs_stereo_npe = new TH1D("h1_proton_trg_pate_vs_stereo_npe","h1_proton_trg_pate_vs_stereo_npe",npe_trg_max+1,-0.5,npe_trg_max+0.5);
  TH1D *h1_proton_trg_pate_vs_OR_npe = new TH1D("h1_proton_trg_pate_vs_OR_npe","h1_proton_trg_pate_vs_OR_npe",npe_trg_max+1,-0.5,npe_trg_max+0.5);
  TH1D *h1_proton_trg_pate_vs_LST1_npe = new TH1D("h1_proton_trg_pate_vs_LST1_npe","h1_proton_trg_pate_vs_LST1_npe",npe_trg_max+1,-0.5,npe_trg_max+0.5);
  TH1D *h1_proton_trg_pate_vs_LST1_npe_PMT = new TH1D("h1_proton_trg_pate_vs_LST1_npe_PMT","h1_proton_trg_pate_vs_LST1_npe_PMT",npe_trg_max+1,-0.5,npe_trg_max+0.5);
  //
  TH1D *h1_n_pe_LST1 = new TH1D("h1_n_pe_LST1","h1_n_pe_LST1",npe_trg_max+1,-0.5,npe_trg_max+0.5);
  TH1D *h1_n_pe_LST2 = new TH1D("h1_n_pe_LST2","h1_n_pe_LST2",npe_trg_max+1,-0.5,npe_trg_max+0.5);
  TH1D *h1_n_pe_LST3 = new TH1D("h1_n_pe_LST3","h1_n_pe_LST3",npe_trg_max+1,-0.5,npe_trg_max+0.5);
  TH1D *h1_n_pe_LST4 = new TH1D("h1_n_pe_LST4","h1_n_pe_LST4",npe_trg_max+1,-0.5,npe_trg_max+0.5);
  //
  Double_t PDE_LST1;
  Double_t PDE_LST2;
  Double_t PDE_LST3;
  Double_t PDE_LST4;
  //
  TH1D *h1_PDE_LST1 = new TH1D ("h1_PDE_LST1","h1_PDE_LST1",200, 0.0,1.1);
  TH1D *h1_PDE_LST2 = new TH1D ("h1_PDE_LST2","h1_PDE_LST2",200, 0.0,1.1);  
  TH1D *h1_PDE_LST3 = new TH1D ("h1_PDE_LST3","h1_PDE_LST3",200, 0.0,1.1);
  TH1D *h1_PDE_LST4 = new TH1D ("h1_PDE_LST4","h1_PDE_LST4",200, 0.0,1.1);  
  TH1D *h1_PDE_LST1_PMT = new TH1D ("h1_PDE_LST1_PMT","h1_PDE_LST1_PMT",200, 0.0,1.1);
  TH1D *h1_PDE_LST2_PMT = new TH1D ("h1_PDE_LST2_PMT","h1_PDE_LST2_PMT",200, 0.0,1.1);  
  TH1D *h1_PDE_LST3_PMT = new TH1D ("h1_PDE_LST3_PMT","h1_PDE_LST3_PMT",200, 0.0,1.1);
  TH1D *h1_PDE_LST4_PMT = new TH1D ("h1_PDE_LST4_PMT","h1_PDE_LST4_PMT",200, 0.0,1.1);  
  //
  Double_t part_energy_Tev;
  Double_t part_energy_GeV;
  //
  Long64_t nentries = fChain->GetEntriesFast();
  cout<<"nentries = "<<nentries<<endl;
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    if(jentry%_anaConf.jentry_modulo == 0)
      cout<<jentry<<endl;
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    //if (ientry > 10) break;
    nb = fChain->GetEntry(jentry); nbytes += nb;
    //
    part_energy_Tev = energy;
    part_energy_GeV = energy*1000;
    //
    h1_n_tot_event->Fill(0.0);
    h1_n_tot_event_w->Fill( 0.0, evstHist::get_Weight_ETeV(part_energy_Tev));
    tot_Weight+=evstHist::get_Weight_ETeV(part_energy_Tev);
    //
    h1_xmax->Fill(xmax);
    h1_hmax->Fill(hmax);
    //
    h1_n_pe_LST1->Fill(n_pe_LST1);
    h1_n_pe_LST2->Fill(n_pe_LST2);
    h1_n_pe_LST3->Fill(n_pe_LST3);
    h1_n_pe_LST4->Fill(n_pe_LST4);
    //
    h2_n_pe_LST2_vs_n_pe_LST1->Fill(n_pe_LST1, n_pe_LST2);
    //
    h1_xcore->Fill(xcore);
    h1_ycore->Fill(ycore);
    //
    if(n_pe_LST1>0){
      h1_xcore_LST1->Fill(xcore);
      h1_ycore_LST1->Fill(ycore);
      h1_ev_time_LST1->Fill(ev_time_LST1);
      if(nphotons_LST1>0){
	PDE_LST1 = (Double_t)n_pe_LST1/(Double_t)nphotons_LST1;
	h1_PDE_LST1->Fill(PDE_LST1);
	h1_PDE_LST1_PMT->Fill(PDE_LST1/2.0);
      }
      else{
	cout<<"ERROR --> n_pe_LST1 > 0 && nphotons_LST1<=0"<<endl
	    <<"          n_pe_LST1 = "<<n_pe_LST1<<endl
	    <<"      nphotons_LST1 = "<<nphotons_LST1<<endl;
	assert(0);
      }
    }
    if(n_pe_LST2>0){
      h1_xcore_LST2->Fill(xcore);
      h1_ycore_LST2->Fill(ycore);
      h1_ev_time_LST2->Fill(ev_time_LST2);
      if(nphotons_LST2>0){
	PDE_LST2 = (Double_t)n_pe_LST2/(Double_t)nphotons_LST2;
	h1_PDE_LST2->Fill(PDE_LST2);	
	h1_PDE_LST2_PMT->Fill(PDE_LST2/2.0);
      }
      else{
	cout<<"ERROR --> n_pe_LST2 > 0 && nphotons_LST2 <= 0"<<endl
	    <<"          n_pe_LST2 = "<<n_pe_LST2<<endl
	    <<"      nphotons_LST2 = "<<nphotons_LST2<<endl;
	assert(0);
      }
    }
    if(n_pe_LST3>0){
      h1_xcore_LST3->Fill(xcore);
      h1_ycore_LST3->Fill(ycore);
      h1_ev_time_LST3->Fill(ev_time_LST3);
      if(nphotons_LST3>0){
	PDE_LST3 = (Double_t)n_pe_LST3/(Double_t)nphotons_LST3;
	h1_PDE_LST3->Fill(PDE_LST3);
	h1_PDE_LST3_PMT->Fill(PDE_LST3/2.0);
      }
      else{
	cout<<"ERROR --> n_pe_LST3 > 0 && nphotons_LST3 <= 0"<<endl
	    <<"          n_pe_LST3 = "<<n_pe_LST3<<endl
	    <<"      nphotons_LST3 = "<<nphotons_LST3<<endl;
	assert(0);
      }
    }
    if(n_pe_LST4>0){
      h1_xcore_LST4->Fill(xcore);
      h1_ycore_LST4->Fill(ycore);
      h1_ev_time_LST4->Fill(ev_time_LST4);
      if(nphotons_LST4>0){
	PDE_LST4 = (Double_t)n_pe_LST4/(Double_t)nphotons_LST4;
	h1_PDE_LST4->Fill(PDE_LST4);
	h1_PDE_LST4_PMT->Fill(PDE_LST4/2.0);
      }
      else{
	cout<<"ERROR --> n_pe_LST4 > 0 && nphotons_LST4 <= 0"<<endl
	    <<"          n_pe_LST4 = "<<n_pe_LST4<<endl
	    <<"      nphotons_LST4 = "<<nphotons_LST4<<endl;
	assert(0);
      }
    }
    //
    if(n_pe_LST1>0 && n_pe_LST2>0){
      h1_dtime_LST1_m_LST2->Fill(ev_time_LST1 - ev_time_LST2);
    }
    if(n_pe_LST1>0 && n_pe_LST3>0){
      h1_dtime_LST1_m_LST3->Fill(ev_time_LST1 - ev_time_LST3);
      h1_dtime_LST1_LST3_m_exp->Fill(TMath::Abs(ev_time_LST1 - ev_time_LST3 - -212.2),evstHist::get_Weight_ETeV(part_energy_Tev));
      h1_dtime_LST1_LST3_m_exp_full->Fill(ev_time_LST1 - ev_time_LST3 - -212.2,evstHist::get_Weight_ETeV(part_energy_Tev));
    }
    if(n_pe_LST1>0 && n_pe_LST4>0){
      h1_dtime_LST1_m_LST4->Fill(ev_time_LST1 - ev_time_LST4);
    }
    if(n_pe_LST2>0 && n_pe_LST3>0){
      h1_dtime_LST2_m_LST3->Fill(ev_time_LST2 - ev_time_LST3);
    }
    if(n_pe_LST2>0 && n_pe_LST4>0){
      h1_dtime_LST2_m_LST4->Fill(ev_time_LST2 - ev_time_LST4);
    }
    if(n_pe_LST3>0 && n_pe_LST4>0){
      h1_dtime_LST3_m_LST4->Fill(ev_time_LST3 - ev_time_LST4);
    }
    //
    h1_n_pe_LST_all->Fill(0);
    //
    h2_ycore_vs_xcore_all->Fill( xcore, ycore);
    //
    if(if_single_LST(coincidenceID)){
      h1_n_pe_LST_single->Fill(coincidenceID);
      if(if_double_LST(coincidenceIDtmp) || if_triple_LST(coincidenceIDtmp) || if_four_LST(coincidenceIDtmp))
	cout<<"(if_double_LST(coincidenceIDtmp) || if_triple_LST(coincidenceIDtmp) || if_four_LST(coincidenceIDtmp)"<<endl;
      //
      if(coincidenceID==0)
	h2_ycore_vs_xcore_single_LST1->Fill( xcore, ycore);
    }
    //
    if(if_double_LST(coincidenceID)){
      h1_n_pe_LST_double->Fill(coincidenceID);
      if(if_single_LST(coincidenceIDtmp) || if_triple_LST(coincidenceIDtmp) || if_four_LST(coincidenceIDtmp))
	cout<<"(if_single_LST(coincidenceIDtmp) || if_triple_LST(coincidenceIDtmp) || if_four_LST(coincidenceIDtmp))"<<endl;
      //
      if(coincidenceID==2)
	h2_ycore_vs_xcore_double_LST1LST4->Fill( xcore, ycore);
    }
    //
    if(if_triple_LST(coincidenceID)){
      h1_n_pe_LST_triple->Fill(coincidenceID);
      if(if_single_LST(coincidenceIDtmp) || if_double_LST(coincidenceIDtmp) || if_four_LST(coincidenceIDtmp))
	cout<<"(if_single_LST(coincidenceIDtmp) || if_double_LST(coincidenceIDtmp) || if_four_LST(coincidenceIDtmp))"<<endl;
      if(coincidenceID==0)
	h2_ycore_vs_xcore_triple_LST1LST2LST3->Fill( xcore, ycore);
    }
    //
    if(if_four_LST(coincidenceID)){
      h1_n_pe_LST_four->Fill(coincidenceID);
      if(if_single_LST(coincidenceIDtmp) || if_double_LST(coincidenceIDtmp) || if_triple_LST(coincidenceIDtmp))
	cout<<"(if_single_LST(coincidenceIDtmp) || if_double_LST(coincidenceIDtmp) || if_triple_LST(coincidenceIDtmp))"<<endl;
      if(coincidenceID==0)
	h2_ycore_vs_xcore_four_LST->Fill( xcore, ycore);
    }
    //
    h2_xmax_vs_hmax->Fill( hmax, xmax);
    //
    if(n_pe_LST1>0){
      sipm_cam->Fill_pix_x_y_hist( n_pe_LST1, pe_chID_LST1, h1_LST1_sipm_cam_x, h1_LST1_sipm_cam_y);
      sipm_cam->Fill_pix_hist2D_y_vs_x(n_pe_LST1, pe_chID_LST1, h2_LST1_sipm_cam_y_vs_x);      
      h1_n_pe_LST_all_LST01->Fill(0);
    }
    if(n_pe_LST2>0){
      sipm_cam->Fill_pix_x_y_hist( n_pe_LST2, pe_chID_LST2, h1_LST2_sipm_cam_x, h1_LST2_sipm_cam_y);    
      sipm_cam->Fill_pix_hist2D_y_vs_x(n_pe_LST2, pe_chID_LST2, h2_LST2_sipm_cam_y_vs_x);      
      h1_n_pe_LST_all_LST02->Fill(0);
    }
    if(n_pe_LST3>0){
      sipm_cam->Fill_pix_x_y_hist( n_pe_LST3, pe_chID_LST3, h1_LST3_sipm_cam_x, h1_LST3_sipm_cam_y);
      sipm_cam->Fill_pix_hist2D_y_vs_x(n_pe_LST3, pe_chID_LST3, h2_LST3_sipm_cam_y_vs_x);      
      h1_n_pe_LST_all_LST03->Fill(0);
    }
    if(n_pe_LST4>0){
      sipm_cam->Fill_pix_x_y_hist( n_pe_LST4, pe_chID_LST4, h1_LST4_sipm_cam_x, h1_LST4_sipm_cam_y);    
      sipm_cam->Fill_pix_hist2D_y_vs_x(n_pe_LST4, pe_chID_LST4, h2_LST4_sipm_cam_y_vs_x);      
      h1_n_pe_LST_all_LST04->Fill(0);
    }
    //LST1 and LST2
    if(n_pe_LST1 > 0 && n_pe_LST2 > 0){
      sipm_cam->get_pix_mean( n_pe_LST1, pe_chID_LST1, x_mean_LST1, y_mean_LST1);
      sipm_cam->get_pix_mean( n_pe_LST2, pe_chID_LST2, x_mean_LST2, y_mean_LST2);
      //
      dr_LST1_LST2 = TMath::Sqrt((x_mean_LST1 - x_mean_LST2)*(x_mean_LST1 - x_mean_LST2) +
				 (y_mean_LST1 - y_mean_LST2)*(y_mean_LST1 - y_mean_LST2));
      dalpha_LST1_LST2 = TMath::ATan(dr_LST1_LST2/28.0)*180.0/TMath::Pi();
      h1_dr_LST1_LST2->Fill(dr_LST1_LST2);
      h1_dalpha_LST1_LST2->Fill(dalpha_LST1_LST2, evstHist::get_Weight_ETeV(part_energy_Tev));
      h2_y_vs_x_LST1_LST2->Fill((x_mean_LST1 - x_mean_LST2),(y_mean_LST1 - y_mean_LST2));

      h2_dt_vs_x_LST1_LST2->Fill((x_mean_LST1 - x_mean_LST2),(ev_time_LST1 - ev_time_LST2 - -75.145079));
      
    }
    //LST1 and LST3
    if(n_pe_LST1 > 0 && n_pe_LST3 > 0){
      sipm_cam->get_pix_mean( n_pe_LST1, pe_chID_LST1, x_mean_LST1, y_mean_LST1);
      sipm_cam->get_pix_mean( n_pe_LST3, pe_chID_LST3, x_mean_LST3, y_mean_LST3);
      //
      dr_LST1_LST3 = TMath::Sqrt((x_mean_LST1 - x_mean_LST3)*(x_mean_LST1 - x_mean_LST3) +
				 (y_mean_LST1 - y_mean_LST3)*(y_mean_LST1 - y_mean_LST3));
      dalpha_LST1_LST3 = TMath::ATan(dr_LST1_LST3/28.0)*180.0/TMath::Pi();
      h1_dalpha_LST1_LST3->Fill(dalpha_LST1_LST3, evstHist::get_Weight_ETeV(part_energy_Tev));
      //
      h2_dalpha_vs_dtLST1_LST3->Fill((ev_time_LST1 - ev_time_LST3 - -212.2),dalpha_LST1_LST3);

      //     h2_x_vs_dt_LST1_LST3->Fill((ev_time_LST1 - ev_time_LST3 - -212.2),(x_mean_LST1 - x_mean_LST3));
      h2_x_vs_dt_LST1_LST3->Fill((ev_time_LST1 - ev_time_LST3 - -212.2),(x_mean_LST1 - x_mean_LST3));

      h2_y_vs_x_LST1_LST3->Fill((x_mean_LST1 - x_mean_LST3),(y_mean_LST1 - y_mean_LST3));
      h2_y_vs_x_LST1_LST3_norm_dt->Fill((x_mean_LST1 - x_mean_LST3),(y_mean_LST1 - y_mean_LST3),(ev_time_LST1 - ev_time_LST3 - -212.2));
      
      //
      
      //
    }
    //LST1 and LST4
    if(n_pe_LST1 > 0 && n_pe_LST4 > 0){
      sipm_cam->get_pix_mean( n_pe_LST1, pe_chID_LST1, x_mean_LST1, y_mean_LST1);
      sipm_cam->get_pix_mean( n_pe_LST4, pe_chID_LST4, x_mean_LST4, y_mean_LST4);
      //
      dr_LST1_LST4 = TMath::Sqrt((x_mean_LST1 - x_mean_LST4)*(x_mean_LST1 - x_mean_LST4) +
				 (y_mean_LST1 - y_mean_LST4)*(y_mean_LST1 - y_mean_LST4));
      dalpha_LST1_LST4 = TMath::ATan(dr_LST1_LST4/28.0)*180.0/TMath::Pi();
      h1_dalpha_LST1_LST4->Fill(dalpha_LST1_LST4, evstHist::get_Weight_ETeV(part_energy_Tev));
      h2_y_vs_x_LST1_LST4->Fill((x_mean_LST1 - x_mean_LST4),(y_mean_LST1 - y_mean_LST4));
    }
    //
    for(Int_t jj_npe_trg = 0; jj_npe_trg < npe_trg_max; jj_npe_trg++){
      if(if_stereo_trigger_pne(jj_npe_trg)){
	h1_npe_trg_max->Fill(jj_npe_trg);
	h1_npe_trg_max_w->Fill(jj_npe_trg, evstHist::get_Weight_ETeV(part_energy_Tev));
      }
      if(n_pe_LST1>=jj_npe_trg){
	h1_npe_trg_max_LST1_w->Fill(jj_npe_trg, evstHist::get_Weight_ETeV(part_energy_Tev));
      }
      if( n_pe_LST1 >= jj_npe_trg || n_pe_LST2 >= jj_npe_trg ||
	  n_pe_LST3 >= jj_npe_trg || n_pe_LST4 >= jj_npe_trg ){
	h1_npe_trg_max_OR_w->Fill(jj_npe_trg, evstHist::get_Weight_ETeV(part_energy_Tev));
      }
      if(n_pe_LST1/2.0>=jj_npe_trg){
	h1_npe_trg_max_LST1_PMT_w->Fill(jj_npe_trg, evstHist::get_Weight_ETeV(part_energy_Tev));
      }
    }
    //
    //LST2 and LST3
    if(n_pe_LST2 > 0 && n_pe_LST3 > 0){
      sipm_cam->get_pix_mean( n_pe_LST2, pe_chID_LST2, x_mean_LST2, y_mean_LST2);
      sipm_cam->get_pix_mean( n_pe_LST3, pe_chID_LST3, x_mean_LST3, y_mean_LST3);
      h2_y_vs_x_LST2_LST3->Fill((x_mean_LST2 - x_mean_LST3),(y_mean_LST2 - y_mean_LST3));
    }
    //LST2 and LST4
    if(n_pe_LST2 > 0 && n_pe_LST4 > 0){
      sipm_cam->get_pix_mean( n_pe_LST2, pe_chID_LST2, x_mean_LST2, y_mean_LST2);
      sipm_cam->get_pix_mean( n_pe_LST4, pe_chID_LST4, x_mean_LST4, y_mean_LST4);
      h2_y_vs_x_LST2_LST4->Fill((x_mean_LST2 - x_mean_LST4),(y_mean_LST2 - y_mean_LST4));
    }
    //LST3 and LST4
    if(n_pe_LST3 > 0 && n_pe_LST4 > 0){
      sipm_cam->get_pix_mean( n_pe_LST3, pe_chID_LST3, x_mean_LST3, y_mean_LST3);
      sipm_cam->get_pix_mean( n_pe_LST4, pe_chID_LST4, x_mean_LST4, y_mean_LST4);
      h2_y_vs_x_LST3_LST4->Fill((x_mean_LST3 - x_mean_LST4),(y_mean_LST3 - y_mean_LST4));
    }    
  }
  //
  //----------------------
  //
  cout<<"LST1 -> LST2 "<<getExpectedTimeDelayBetweenTwoLST( LST1_r0, LST2_r0, TMath::Pi(), TMath::Pi()/2.0 - 20.0/180.0*TMath::Pi())<<endl;
  cout<<"LST1 -> LST3 "<<getExpectedTimeDelayBetweenTwoLST( LST1_r0, LST3_r0, TMath::Pi(), TMath::Pi()/2.0 - 20.0/180.0*TMath::Pi())<<endl;
  cout<<"LST1 -> LST4 "<<getExpectedTimeDelayBetweenTwoLST( LST1_r0, LST4_r0, TMath::Pi(), TMath::Pi()/2.0 - 20.0/180.0*TMath::Pi())<<endl;
  cout<<"LST2 -> LST3 "<<getExpectedTimeDelayBetweenTwoLST( LST2_r0, LST3_r0, TMath::Pi(), TMath::Pi()/2.0 - 20.0/180.0*TMath::Pi())<<endl;
  cout<<"LST2 -> LST4 "<<getExpectedTimeDelayBetweenTwoLST( LST2_r0, LST4_r0, TMath::Pi(), TMath::Pi()/2.0 - 20.0/180.0*TMath::Pi())<<endl;
  cout<<"LST3 -> LST4 "<<getExpectedTimeDelayBetweenTwoLST( LST3_r0, LST4_r0, TMath::Pi(), TMath::Pi()/2.0 - 20.0/180.0*TMath::Pi())<<endl;
  //
  //
  //
  anabase::TH1D_divide( h1_dtime_LST1_m_LST2, h1_dtime_LST1_m_LST2_norm, h1_dtime_LST1_m_LST2->GetMaximum());
  anabase::TH1D_divide( h1_dtime_LST1_m_LST3, h1_dtime_LST1_m_LST3_norm, h1_dtime_LST1_m_LST3->GetMaximum());
  anabase::TH1D_divide( h1_dtime_LST1_m_LST4, h1_dtime_LST1_m_LST4_norm, h1_dtime_LST1_m_LST4->GetMaximum());
  anabase::TH1D_divide( h1_dtime_LST2_m_LST3, h1_dtime_LST2_m_LST3_norm, h1_dtime_LST2_m_LST3->GetMaximum());
  anabase::TH1D_divide( h1_dtime_LST2_m_LST4, h1_dtime_LST2_m_LST4_norm, h1_dtime_LST2_m_LST4->GetMaximum());
  anabase::TH1D_divide( h1_dtime_LST3_m_LST4, h1_dtime_LST3_m_LST4_norm, h1_dtime_LST3_m_LST4->GetMaximum());
  //  
  //cout<<LST2_r0.x()<<endl;
  //
  //
  //----------------------
  //
  //
  anabase::TH2D_divide( h2_ycore_vs_xcore_single_LST1, h2_ycore_vs_xcore_all, h2_ycore_vs_xcore_single_LST1_norm);
  anabase::TH2D_divide( h2_ycore_vs_xcore_double_LST1LST4, h2_ycore_vs_xcore_all,  h2_ycore_vs_xcore_double_LST1LST4_norm);
  anabase::TH2D_divide( h2_ycore_vs_xcore_triple_LST1LST2LST3, h2_ycore_vs_xcore_all,  h2_ycore_vs_xcore_triple_LST1LST2LST3_norm);
  anabase::TH2D_divide( h2_ycore_vs_xcore_four_LST, h2_ycore_vs_xcore_all, h2_ycore_vs_xcore_four_LST_norm);
  //
  //
  //----------------------
  //
  //
  get_event_acceptance_vs_dt( h1_dtime_LST1_LST3_m_exp, h1_dtime_event_acceptance);
  //
  //
  //
  cout<<"_n_file_counter = "<<_n_file_counter<<endl;  
  //
  Double_t simulation_sampling_factor = anaTrgB::calculate_simulation_sampling_factor(_anaConf.eMin, _anaConf.eMax);
  Double_t proton_flux_pers_persr_perm2 = anaTrgB::get_proton_flux_pers_persr_perm2(_anaConf.eMin, _anaConf.eMax);
  Double_t proton_sim_solidAngle = anaTrgB::get_solid_angle(_anaConf.thetaSimDeg);
  Double_t proton_generationArea = TMath::Pi()*_anaConf.rsim*_anaConf.rsim;
  Double_t proton_tot_flux = proton_flux_pers_persr_perm2*proton_sim_solidAngle*proton_generationArea;
  Double_t effective_sim_proton = _n_file_counter*_anaConf.nEvPerFileSim*_anaConf.nShowerMultiplicationFactor*simulation_sampling_factor;
  //Double_t effective_sim_proton = (h1_n_tot_event_w->GetBinContent(2)/h1_n_tot_event->GetBinContent(2)*total_corsica_sim_proton);
  //
  Double_t effective_sim_time = effective_sim_proton/proton_tot_flux;
  //
  cout<<"proton_flux_pers_persr_perm2 "<<proton_flux_pers_persr_perm2<<endl
      <<"proton_sim_solidAngle        "<<proton_sim_solidAngle<<endl
      <<"proton_generationArea        "<<proton_generationArea<<endl
      <<"proton_tot_flux              "<<proton_tot_flux<<endl
      <<"simulation_sampling_factor   "<<simulation_sampling_factor<<endl;
  //
  cout<<"h1_n_tot_event->GetBinContent(2)   = "<<h1_n_tot_event->GetBinContent(2)<<endl
      <<"h1_n_tot_event_w->GetBinContent(2) = "<<h1_n_tot_event_w->GetBinContent(2)<<endl
      <<"tot_Weight                         = "<<tot_Weight<<endl
      <<"effective_sim_proton               = "<<effective_sim_proton<<endl
      <<"effective_sim_time                 = "<<effective_sim_time<<endl;
  //
  anabase::TH1D_divide(h1_npe_trg_max_w, h1_proton_trg_pate_vs_stereo_npe, effective_sim_time);
  anabase::TH1D_divide(h1_npe_trg_max_OR_w, h1_proton_trg_pate_vs_OR_npe, effective_sim_time);
  anabase::TH1D_divide(h1_npe_trg_max_LST1_w, h1_proton_trg_pate_vs_LST1_npe, effective_sim_time);
  anabase::TH1D_divide(h1_npe_trg_max_LST1_PMT_w, h1_proton_trg_pate_vs_LST1_npe_PMT, effective_sim_time);
  //
  //
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
  h1_n_pe_LST1->Write();
  h1_n_pe_LST2->Write();
  h1_n_pe_LST3->Write();
  h1_n_pe_LST4->Write();
  //
  h2_n_pe_LST2_vs_n_pe_LST1->Write();
  //
  h1_xcore->Write();
  h1_ycore->Write();
  h1_xcore_LST1->Write();
  h1_ycore_LST1->Write();
  h1_xcore_LST2->Write();
  h1_ycore_LST2->Write();
  h1_xcore_LST3->Write();
  h1_ycore_LST3->Write();
  h1_xcore_LST4->Write();
  h1_ycore_LST4->Write();
  //
  h1_ev_time_LST1->Write();
  h1_ev_time_LST2->Write();
  h1_ev_time_LST3->Write();
  h1_ev_time_LST4->Write();
  //
  h1_dtime_LST1_m_LST2->Write();
  h1_dtime_LST1_m_LST3->Write();
  h1_dtime_LST1_m_LST4->Write();
  h1_dtime_LST2_m_LST3->Write();
  h1_dtime_LST2_m_LST4->Write();
  h1_dtime_LST3_m_LST4->Write();
  //
  h1_dtime_LST1_m_LST2_norm->Write();
  h1_dtime_LST1_m_LST3_norm->Write();
  h1_dtime_LST1_m_LST4_norm->Write();
  h1_dtime_LST2_m_LST3_norm->Write();
  h1_dtime_LST2_m_LST4_norm->Write();
  h1_dtime_LST3_m_LST4_norm->Write();
  //
  h1_dtime_LST1_LST3_m_exp->Write();
  h1_dtime_LST1_LST3_m_exp_full->Write();
  h1_dtime_event_acceptance->Write();
  //
  h1_n_pe_LST_all->Write();
  h1_n_pe_LST_all_LST01->Write();
  h1_n_pe_LST_all_LST02->Write();
  h1_n_pe_LST_all_LST03->Write();
  h1_n_pe_LST_all_LST04->Write();
  //
  h1_n_pe_LST_single->Write();
  h1_n_pe_LST_double->Write();
  h1_n_pe_LST_triple->Write();
  h1_n_pe_LST_four->Write();
  //
  h1_xmax->Write();
  h1_hmax->Write();
  //
  h2_xmax_vs_hmax->Write();
  //
  h1_LST1_sipm_cam_x->Write();
  h1_LST2_sipm_cam_x->Write();
  h1_LST3_sipm_cam_x->Write();
  h1_LST4_sipm_cam_x->Write();
  //
  h1_LST1_sipm_cam_y->Write();
  h1_LST2_sipm_cam_y->Write();
  h1_LST3_sipm_cam_y->Write();
  h1_LST4_sipm_cam_y->Write();
  //
  h2_LST1_sipm_cam_y_vs_x->Write();
  h2_LST2_sipm_cam_y_vs_x->Write();
  h2_LST3_sipm_cam_y_vs_x->Write();
  h2_LST4_sipm_cam_y_vs_x->Write();
  //
  h1_dr_LST1_LST2->Write();
  h1_dalpha_LST1_LST2->Write();
  h1_dalpha_LST1_LST3->Write();
  h1_dalpha_LST1_LST4->Write();
  //
  h2_ycore_vs_xcore_all->Write();
  h2_ycore_vs_xcore_single_LST1->Write();
  h2_ycore_vs_xcore_double_LST1LST4->Write();
  h2_ycore_vs_xcore_triple_LST1LST2LST3->Write();
  h2_ycore_vs_xcore_four_LST->Write();
  //
  h2_ycore_vs_xcore_single_LST1_norm->Write();
  h2_ycore_vs_xcore_double_LST1LST4_norm->Write();
  h2_ycore_vs_xcore_triple_LST1LST2LST3_norm->Write();
  h2_ycore_vs_xcore_four_LST_norm->Write();
  //
  h2_dalpha_vs_dtLST1_LST3->Write();
  h2_x_vs_dt_LST1_LST3->Write();
  //
  h2_y_vs_x_LST1_LST2->Write();
  h2_y_vs_x_LST1_LST3->Write();
  h2_y_vs_x_LST1_LST3_norm_dt->Write();
  h2_y_vs_x_LST1_LST4->Write();
  //
  h2_y_vs_x_LST2_LST3->Write();
  h2_y_vs_x_LST2_LST4->Write();
  h2_y_vs_x_LST3_LST4->Write();
  //  
  //
  h2_dt_vs_x_LST1_LST2->Write();
  //
  h1_PDE_LST1->Write();
  h1_PDE_LST2->Write();
  h1_PDE_LST3->Write();
  h1_PDE_LST4->Write();
  //
  h1_PDE_LST1_PMT->Write();
  h1_PDE_LST2_PMT->Write();
  h1_PDE_LST3_PMT->Write();
  h1_PDE_LST4_PMT->Write();
  //
  h1_npe_trg_max->Write();
  h1_npe_trg_max_w->Write();
  //
  h1_n_tot_event->Write();
  h1_n_tot_event_w->Write();
  h1_npe_trg_max_OR_w->Write();
  h1_npe_trg_max_LST1_w->Write();
  h1_npe_trg_max_LST1_PMT_w->Write();
  //
  h1_proton_trg_pate_vs_stereo_npe->Write();
  h1_proton_trg_pate_vs_LST1_npe->Write();
  h1_proton_trg_pate_vs_LST1_npe_PMT->Write();
  h1_proton_trg_pate_vs_OR_npe->Write();
  //
  rootFile->Close();
}

void anastereo::get_event_acceptance_vs_dt( TH1D *h1_dtime, TH1D *h1_event){
  //Double_t ntot_ev = h1_dtime->GetEntries();
  Double_t ntot_ev = h1_dtime->Integral(1,h1_dtime->GetNbinsX());
  Double_t dt;
  Double_t t;
  Double_t nev = 0;
  for(Int_t i = 1;i<=h1_event->GetNbinsX();i++){
    dt = h1_event->GetBinCenter(i);
    nev = 0;
    for(Int_t j = 1;j<=h1_dtime->GetNbinsX();j++){
      t = h1_dtime->GetBinCenter(j);
      if(t<dt/2)
	nev = nev + h1_dtime->GetBinContent(j);
      else
	break;
    }
    h1_event->SetBinContent(i,nev);
  }
  for(Int_t i = 1;i<=h1_event->GetNbinsX();i++){
    h1_event->SetBinContent(i,h1_event->GetBinContent(i)/ntot_ev);
  }
}
