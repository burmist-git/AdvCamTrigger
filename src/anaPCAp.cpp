//my
#include "anaPCAp.hh"
#include "sipmCameraHist.hh"
#include "sipmCameraHistCropped.hh"
#include "triggerSim.hh"
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

anaPCAp::anaPCAp(TString fileList, TString anaConfFile) : anashort(fileList), _phi0_shift(77.4)
{
  _anaConf.readFromFile(anaConfFile);
}
							
void anaPCAp::Loop(TString histOut){
  //
  _anaConf.printInfo();
  //
  //TRandom3 *rnd = new TRandom3(123123);
  //
  TH1D *h1_nphotons = new TH1D("h1_nphotons","h1_nphotons",1000,0.0,10000);
  TH1D *h1_n_pe = new TH1D("h1_n_pe","h1_n_pe",10001,-0.5,10000.5);
  TH1D *h1_n_pixels = new TH1D("h1_n_pixels","h1_n_pixels",1000,0.0,10000);
  TH1D *h1_azimuth = new TH1D("h1_azimuth","azimuth",400,-4,4);
  TH1D *h1_altitude = new TH1D("h1_altitude","h1_altitude",400,4,4);
  //
  //
  std::vector<sipmCameraHistCropped*> simp_hist_crop_v;
  TString sipm_hist_crop_name;
  //
  Double_t x0_LST01 = -70.93;
  Double_t y0_LST01 = -52.07;
  Double_t tel_theta = 20.0/180.0*TMath::Pi();
  Double_t tel_phi = 180.0/180.0*TMath::Pi();
  //
  Double_t r_core = 0.0;
  Double_t theta_core = 0.0;
  //
  Double_t x_shift;
  Double_t y_shift;
  //
  Int_t n_ev_cuts = 0;
  //
  //
  Int_t fadc_sum_offset = 15;
  //Float_t fadc_amplitude = 8.25;
  Int_t fadc_MHz = 1024;
  Int_t fadc_offset = 300;
  //Int_t fadc_offset = 0;
  Float_t fadc_sample_in_ns = 1000.0/fadc_MHz;
  Float_t time_offset = fadc_sum_offset*fadc_sample_in_ns;
  //
  sipmCameraHist *sipm_cam = new sipmCameraHist("sipm_cam","sipm_cam","pixel_mapping.csv",0);
  sipmCameraHist *sipm_cam_norot = new sipmCameraHist("sipm_cam_norot","sipm_cam_norot","pixel_mapping.csv",0);
  sipmCameraHist *sipm_cam_shift = new sipmCameraHist("sipm_cam_shift","sipm_cam_shift","pixel_mapping.csv",0);
  sipmCameraHist *sipm_cam_shift_rot = new sipmCameraHist("sipm_cam_shift_rot","sipm_cam_shift_rot","pixel_mapping.csv",0);
  sipmCameraHistCropped *tmp_cam_hist = new sipmCameraHistCropped("tmp","tmp",sipm_cam,"sipmCameraHistCropped_pix.map");
  //
  //////////////////////////
  TVector3 vx_tel(1.0,0.0,0.0);
  TVector3 vy_tel(0.0,1.0,0.0);
  TVector3 vz_tel(0.0,0.0,1.0);
  sipmCameraHist::get_tel_frame( tel_theta, tel_phi, vx_tel, vy_tel, vz_tel);
  //cout<<"vx_tel.x() "<<vx_tel.x()<<endl
  //  <<"vx_tel.y() "<<vx_tel.y()<<endl
  //  <<"vx_tel.z() "<<vx_tel.z()<<endl;    
  //////////////////////////
  //
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  if(_anaConf.nentries_max>0)
    nentries = _anaConf.nentries_max;
  cout<<"nentries = "<<nentries<<endl;
  //
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    //
    if(jentry%_anaConf.jentry_modulo == 0)
      cout<<jentry<<endl;
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    //
    nb = fChain->GetEntry(jentry); nbytes += nb;
    //
    getCore_rel_R_theta( x0_LST01, y0_LST01, xcore, ycore, r_core, theta_core);
    _r_core = r_core;
    _theta_core = theta_core;
    //
    if(cuts()){
      n_ev_cuts++;
      //
      if(n_ev_cuts>_anaConf.n_ev_cuts_max && _anaConf.n_ev_cuts_max>0)
	break;
      h1_nphotons->Fill(nphotons);
      h1_n_pe->Fill(n_pe);
      h1_n_pixels->Fill(n_pixels);
      //
      h1_azimuth->Fill(azimuth);
      h1_altitude->Fill(altitude);
      //
      //
      sipmCameraHist::get_x_y_shift( vx_tel, vy_tel, vz_tel, azimuth, altitude, x_shift, y_shift, _phi0_shift);
      //cout<<"sipmCameraHist::get_theta_p_t_anaFast( azimuth, altitude) = "<<sipmCameraHist::get_theta_p_t_anaFast( azimuth, altitude)*180.0/TMath::Pi()<<endl;
      //
      //
      //sipm_cam->Fill_pe(n_pe, pe_chID);
      sipm_cam->Fill_pe(n_pe, pe_chID, -theta_core);
      sipm_cam_norot->Fill_pe(n_pe, pe_chID);
      //sipm_cam_shift->Fill_pe(n_pe, pe_chID, 0.0, -0.3857, -0.0434);
      //----------------|------------|------------|------------|
      //     phi0       |     8.5    |     90     |    -90     |
      //----------------|------------|------------|------------|
      //x_shift  0.3857 | -0.048154  | -0.404713  |  0.404713  |
      //y_shift  0.0434 |  0.402011  | 0.0117959  | -0.0117959 |
      //----------------|------------|------------|------------|
      sipm_cam_shift->Fill_pe(n_pe, pe_chID, 0.0, -x_shift, -y_shift);
      sipm_cam_shift_rot->Fill_pe(n_pe, pe_chID, -theta_core, -x_shift, -y_shift);
      //sipm_cam_shift->Fill_pe(n_pe, pe_chID, 0.0, 0.0, 0.0);
      //cout<<"x_shift "<<x_shift<<endl
      //  <<"y_shift "<<y_shift<<endl;
      //
      //sipm_cam->Fill_pe(n_pe, pe_chID, -theta_core, h1_theta_pix_rot, h1_theta_deg_pix_rot, h1_r_pix_rot);
      //
      sipm_hist_crop_name = "sipm_cam_crop";
      sipm_hist_crop_name += "_ev";
      sipm_hist_crop_name += n_ev_cuts;
      sipmCameraHistCropped* simp_hist_crop_tmp = new sipmCameraHistCropped(sipm_hist_crop_name.Data(),
									    sipm_hist_crop_name.Data(),
									    sipm_cam,
									    tmp_cam_hist->get_pixel_map());
      //simp_hist_crop_tmp->Fill_pe(n_pe, pe_chID, pe_time, ev_time, time_offset, -theta_core, principal, false, rnd->Gaus(0.0,0.5), rnd->Uniform(-2.5,-0.5));
      //  static void get_x_y_shift( Double_t azimuth, Double_t altitude, Double_t &x_shift, Double_t &y_shift);
      simp_hist_crop_tmp->Fill_pe(n_pe, pe_chID, pe_time, ev_time, time_offset, -theta_core, NULL, false, -x_shift, -y_shift);
      //simp_hist_crop_tmp->Fill_pe(n_pe, pe_chID, pe_time, ev_time, time_offset, 0.0, NULL, false, 0.0, 0.0);
      simp_hist_crop_v.push_back(simp_hist_crop_tmp);
      if(simp_hist_crop_v.size() == 1000){
	tmp_cam_hist->Save_to_csv("data_non_normalized.csv",simp_hist_crop_v);
	for(unsigned int iii = 0; iii<simp_hist_crop_v.size();iii++)
	  delete simp_hist_crop_v.at(iii);
	simp_hist_crop_v.clear();
      }
    }
  }
  //
  tmp_cam_hist->Save_to_csv("data_non_normalized.csv",simp_hist_crop_v);
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
  for(unsigned int ii = 0;ii<simp_hist_crop_v.size();ii++)
    simp_hist_crop_v.at(ii)->Write(); 
  //
  //
  h1_nphotons->Write();
  h1_n_pe->Write();
  h1_n_pixels->Write();
  //
  h1_azimuth->Write();
  h1_altitude->Write();
  //
  sipm_cam->Write();
  sipm_cam_norot->Write();
  sipm_cam_shift->Write();
  sipm_cam_shift_rot->Write();
  //
  TCanvas *c2 = NULL;
  c2 = new TCanvas("c2","c2",10,10,1500,1000);
  tmp_cam_hist->draw_crop_vector( 11, 3, simp_hist_crop_v, c2);
  c2->Write();
  //
  rootFile->Close();
  //
  cout<<"n_ev_cuts = "<<n_ev_cuts<<endl;
}

void anaPCAp::draw_principal(TString histOut){
  //
  load_S_Vh_data("./S.cvs","./Vh.cvs");
  //
  //
  sipmCameraHist *sipm_cam = new sipmCameraHist("sipm_cam","sipm_cam","pixel_mapping.csv",0);
  sipmCameraHistCropped *tmp_cam_hist = new sipmCameraHistCropped("tmp","tmp",sipm_cam,"sipmCameraHistCropped_pix.map");
  //
  std::vector<sipmCameraHistCropped*> sipm_cam_principal_hist_v;
  tmp_cam_hist->Fill_principal( sipm_cam_principal_hist_v, _data_Vh);
  //
  TCanvas *c1 = NULL;
  c1 = new TCanvas("c1","c1",10,10,1500,1000);
  tmp_cam_hist->draw_crop_vector( 11, 3, sipm_cam_principal_hist_v, c1);
  //
  //cout<<"sipm_cam_principal_hist_v.size() = "<<sipm_cam_principal_hist_v.size()<<endl;
  TFile* rootFile = new TFile(histOut.Data(), "RECREATE", " Histograms", 1);
  rootFile->cd();
  if (rootFile->IsZombie()){
    cout<<"  ERROR ---> file "<<histOut.Data()<<" is zombi"<<endl;
    assert(0);
  }
  else
    cout<<"  Output Histos file ---> "<<histOut.Data()<<endl;
  //
  for(unsigned int ii = 0;ii<sipm_cam_principal_hist_v.size();ii++)
    sipm_cam_principal_hist_v.at(ii)->Write();    
  //
  c1->Write();
  rootFile->Close();
}

bool anaPCAp::cuts(){
  //
  if(_anaConf.disable_all_cuts)
    return true;
  //
  if(_anaConf.cuts_set_to_false)
    return false;
  //if((_theta_core*180/TMath::Pi()<10.0) ||
  //  ((_theta_core*180/TMath::Pi()>125.0) && (_theta_core*180/TMath::Pi()<135.0)) ||
  //  ((_theta_core*180/TMath::Pi()>245.0) && (_theta_core*180/TMath::Pi()<255.0)))
  //return true;
  //if(energy>0.1)
  //if(n_pe == 100)
  //return true;
  //if(hmax > 10000 && hmax < 11000)
  //if(n_pe>1000 && n_pe<1100)
  //return true;
  //if(hmax > 10000 && hmax < 11000)
  //if(n_pe>3000 && n_pe<5000)
  //return true;
  //
  Double_t azimuth_min  = (180.0 - 5.0)/180.0*TMath::Pi();
  Double_t azimuth_max  = (180.0 + 5.0)/180.0*TMath::Pi();
  Double_t altitude_min = (90.0 - 20.0 - 2.0)/180.0*TMath::Pi();
  Double_t altitude_max = (90.0 - 20.0 + 2.0)/180.0*TMath::Pi();
  //
  //
  // "GOOD" PCA
  //if(n_pe>1000 && n_pe<2000)
  //return true;
  //
  //
  if(azimuth>azimuth_min && azimuth<azimuth_max)
    if(altitude>altitude_min && altitude<altitude_max)
      if(n_pe>1000 && n_pe<2000)
	return true;
  //
  //
  //
  //if(hmax > 9000 && hmax < 12000)
  //if(n_pe>1000 && n_pe<2000)
  //if(azimuth>azimuth_min && azimuth<azimuth_max)
  //if(altitude>altitude_min && altitude<altitude_max)
  //return true;
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

void anaPCAp::load_S_Vh_data(TString name_S, TString name_Vh){
  //
  ifstream file_cvs_S(name_S);
  ifstream file_cvs_Vh(name_Vh);
  if(file_cvs_S.is_open())
    for(Int_t i = 0; i <_dd_im; i++)
      file_cvs_S>>_data_S[i];
  file_cvs_S.close();
  //
  if(file_cvs_Vh.is_open())
    for(Int_t i = 0; i <_dd_im; i++)
      for(Int_t j = 0; j <_dd_im; j++)
	file_cvs_Vh>>_data_Vh[i][j];
  file_cvs_Vh.close();
}
