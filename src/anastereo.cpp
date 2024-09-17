//my
#include "anastereo.hh"
#include "sipmCameraHist.hh"

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

void anastereo::Loop(TString histOut){
  //
  //
  TH1D *h1_n_pe_LST_all    = new TH1D( "h1_n_pe_LST_all",    "h1_n_pe_LST_all",    3, -1.5, 1.5);
  TH1D *h1_n_pe_LST_single = new TH1D( "h1_n_pe_LST_single", "h1_n_pe_LST_single", 6, -1.5, 4.5);
  TH1D *h1_n_pe_LST_double = new TH1D( "h1_n_pe_LST_double", "h1_n_pe_LST_double", 8, -1.5, 6.5);
  TH1D *h1_n_pe_LST_triple = new TH1D( "h1_n_pe_LST_triple", "h1_n_pe_LST_triple", 6, -1.5, 4.5);
  TH1D *h1_n_pe_LST_four   = new TH1D( "h1_n_pe_LST_four",   "h1_n_pe_LST_four",   3, -1.5, 1.5);
  //
  //
  TH1D *h1_n_pe_LST1 = new TH1D("h1_n_pe_LST1","h1_n_pe_LST1",1000,0.0,1000);
  TH1D *h1_n_pe_LST2 = new TH1D("h1_n_pe_LST2","h1_n_pe_LST2",1000,0.0,1000);
  TH1D *h1_n_pe_LST3 = new TH1D("h1_n_pe_LST3","h1_n_pe_LST3",1000,0.0,1000);
  TH1D *h1_n_pe_LST4 = new TH1D("h1_n_pe_LST4","h1_n_pe_LST4",1000,0.0,1000);
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
  TH1D *h1_ev_time_LST1 = new TH1D("h1_ev_time_LST1","h1_ev_time_LST1",1000,-2000.0,2000);
  TH1D *h1_ev_time_LST2 = new TH1D("h1_ev_time_LST2","h1_ev_time_LST2",1000,-2000.0,2000);
  TH1D *h1_ev_time_LST3 = new TH1D("h1_ev_time_LST3","h1_ev_time_LST3",1000,-2000.0,2000);
  TH1D *h1_ev_time_LST4 = new TH1D("h1_ev_time_LST4","h1_ev_time_LST4",1000,-2000.0,2000);
  //
  TH1D *h1_dtime_LST1_m_LST2 = new TH1D("h1_dtime_LST1_m_LST2","h1_dtime_LST1_m_LST2",1000,-400.0,400);
  TH1D *h1_dtime_LST1_m_LST3 = new TH1D("h1_dtime_LST1_m_LST3","h1_dtime_LST1_m_LST3",1000,-400.0,400);
  TH1D *h1_dtime_LST1_m_LST4 = new TH1D("h1_dtime_LST1_m_LST4","h1_dtime_LST1_m_LST4",1000,-400.0,400);
  //
  TH1D *h1_xmax = new TH1D("h1_xmax","h1_xmax", 1000, 0.0,  1000);
  TH1D *h1_hmax = new TH1D("h1_hmax","h1_hmax", 1000, 0.0, 25000);
  //
  TH2D *h2_xmax_vs_hmax = new TH2D("h2_xmax_vs_hmax","h2_xmax_vs_hmax", 400, 0.0, 25000, 400, 0.0,  1000);
  //
  Int_t coincidenceID;
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
  Long64_t nentries = fChain->GetEntriesFast();
  cout<<"nentries = "<<nentries<<endl;
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    if(jentry%1000 == 0)
      cout<<jentry<<endl;
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    //if (ientry > 10) break;
    nb = fChain->GetEntry(jentry); nbytes += nb;
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
    if(n_pe_LST1>50 && n_pe_LST2>50){
      h1_xcore_LST1->Fill(xcore);
      h1_ycore_LST1->Fill(ycore);
      h1_ev_time_LST1->Fill(ev_time_LST1);
    }
    if(n_pe_LST2>0){
      h1_xcore_LST2->Fill(xcore);
      h1_ycore_LST2->Fill(ycore);
      h1_ev_time_LST2->Fill(ev_time_LST2);
    }
    if(n_pe_LST3>0){
      h1_xcore_LST3->Fill(xcore);
      h1_ycore_LST3->Fill(ycore);
      h1_ev_time_LST3->Fill(ev_time_LST3);
    }
    if(n_pe_LST4>0){
      h1_xcore_LST4->Fill(xcore);
      h1_ycore_LST4->Fill(ycore);
      h1_ev_time_LST4->Fill(ev_time_LST4);
    }
    //
    if(n_pe_LST1>50 && n_pe_LST2>50)
      h1_dtime_LST1_m_LST2->Fill(ev_time_LST1 - ev_time_LST2);
    if(n_pe_LST1>50 && n_pe_LST3>50)
      h1_dtime_LST1_m_LST3->Fill(ev_time_LST1 - ev_time_LST3);
    if(n_pe_LST1>50 && n_pe_LST4>50)
      h1_dtime_LST1_m_LST4->Fill(ev_time_LST1 - ev_time_LST4);
    //
    h1_n_pe_LST_all->Fill(0);
    //
    if(if_single_LST(coincidenceID))
      h1_n_pe_LST_single->Fill(coincidenceID);
    //
    if(if_double_LST(coincidenceID))
      h1_n_pe_LST_double->Fill(coincidenceID);
    //
    if(if_triple_LST(coincidenceID))
      h1_n_pe_LST_triple->Fill(coincidenceID);
    //
    if(if_four_LST(coincidenceID))
      h1_n_pe_LST_four->Fill(coincidenceID);
    //
    h2_xmax_vs_hmax->Fill( hmax, xmax);
    //
    if(n_pe_LST1>0){
      sipm_cam->Fill_pix_x_y_hist( n_pe_LST1, pe_chID_LST1, h1_LST1_sipm_cam_x, h1_LST1_sipm_cam_y);
      sipm_cam->Fill_pix_hist2D_y_vs_x(n_pe_LST1, pe_chID_LST1, h2_LST1_sipm_cam_y_vs_x);      
    }
    if(n_pe_LST2>0){
      sipm_cam->Fill_pix_x_y_hist( n_pe_LST2, pe_chID_LST2, h1_LST2_sipm_cam_x, h1_LST2_sipm_cam_y);    
      sipm_cam->Fill_pix_hist2D_y_vs_x(n_pe_LST2, pe_chID_LST2, h2_LST2_sipm_cam_y_vs_x);      
    }
    if(n_pe_LST3>0){
      sipm_cam->Fill_pix_x_y_hist( n_pe_LST3, pe_chID_LST3, h1_LST3_sipm_cam_x, h1_LST3_sipm_cam_y);
      sipm_cam->Fill_pix_hist2D_y_vs_x(n_pe_LST3, pe_chID_LST3, h2_LST3_sipm_cam_y_vs_x);      
    }
    if(n_pe_LST4>0){
      sipm_cam->Fill_pix_x_y_hist( n_pe_LST4, pe_chID_LST4, h1_LST4_sipm_cam_x, h1_LST4_sipm_cam_y);    
      sipm_cam->Fill_pix_hist2D_y_vs_x(n_pe_LST4, pe_chID_LST4, h2_LST4_sipm_cam_y_vs_x);      
    }
    //LST1 and LST2
    if(n_pe_LST1 > 0 && n_pe_LST2 > 0){
      sipm_cam->get_pix_mean( n_pe_LST1, pe_chID_LST1, x_mean_LST1, y_mean_LST1);
      sipm_cam->get_pix_mean( n_pe_LST2, pe_chID_LST2, x_mean_LST2, y_mean_LST2);
      //
      dr_LST1_LST2 = TMath::Sqrt((x_mean_LST1 - x_mean_LST2)*(x_mean_LST1 - x_mean_LST2) +
				 (y_mean_LST1 - y_mean_LST2)*(y_mean_LST1 - y_mean_LST2));
      dalpha_LST1_LST2 = TMath::ATan(dr_LST1_LST2/28.0)*180.0/TMath::Pi();
      h1_dalpha_LST1_LST2->Fill(dalpha_LST1_LST2);
    }
    //LST1 and LST3
    if(n_pe_LST1 > 0 && n_pe_LST3 > 0){
      sipm_cam->get_pix_mean( n_pe_LST1, pe_chID_LST1, x_mean_LST1, y_mean_LST1);
      sipm_cam->get_pix_mean( n_pe_LST3, pe_chID_LST3, x_mean_LST3, y_mean_LST3);
      //
      dr_LST1_LST3 = TMath::Sqrt((x_mean_LST1 - x_mean_LST3)*(x_mean_LST1 - x_mean_LST3) +
				 (y_mean_LST1 - y_mean_LST3)*(y_mean_LST1 - y_mean_LST3));
      dalpha_LST1_LST3 = TMath::ATan(dr_LST1_LST3/28.0)*180.0/TMath::Pi();
      h1_dalpha_LST1_LST3->Fill(dalpha_LST1_LST3);
    }
    //LST1 and LST4
    if(n_pe_LST1 > 0 && n_pe_LST4 > 0 && n_pe_LST2==0 && n_pe_LST3==0){
      sipm_cam->get_pix_mean( n_pe_LST1, pe_chID_LST1, x_mean_LST1, y_mean_LST1);
      sipm_cam->get_pix_mean( n_pe_LST4, pe_chID_LST4, x_mean_LST4, y_mean_LST4);
      //
      dr_LST1_LST4 = TMath::Sqrt((x_mean_LST1 - x_mean_LST4)*(x_mean_LST1 - x_mean_LST4) +
				 (y_mean_LST1 - y_mean_LST4)*(y_mean_LST1 - y_mean_LST4));
      dalpha_LST1_LST4 = TMath::ATan(dr_LST1_LST4/28.0)*180.0/TMath::Pi();
      h1_dalpha_LST1_LST4->Fill(dalpha_LST1_LST4);
    }
  }
  //
  //----------------------
  //
  cout<<"LST1 -> LST2 "<<getExpectedTimeDelayBetweenTwoLST( LST1_r0, LST2_r0, TMath::Pi(), TMath::Pi()/2.0 - 20.0/180.0*TMath::Pi())<<endl;
  cout<<"LST1 -> LST3 "<<getExpectedTimeDelayBetweenTwoLST( LST1_r0, LST3_r0, TMath::Pi(), TMath::Pi()/2.0 - 20.0/180.0*TMath::Pi())<<endl;
  cout<<"LST1 -> LST4 "<<getExpectedTimeDelayBetweenTwoLST( LST1_r0, LST4_r0, TMath::Pi(), TMath::Pi()/2.0 - 20.0/180.0*TMath::Pi())<<endl;
  //cout<<LST2_r0.x()<<endl;
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
  //
  h1_n_pe_LST_all->Write();
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
  rootFile->Close();
}
