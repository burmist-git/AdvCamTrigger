//my
#include "anaPCA.hh"
#include "sipmCameraHist.hh"
#include "sipmCameraHistCropped.hh"
#include "wfCamSim.hh"
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
#include "TPrincipal.h"

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

anaPCA::anaPCA(TString fileList, TString anaConfFile) : anashort(fileList)
{
  _anaConf.readFromFile(anaConfFile);
}
							
void anaPCA::Loop(TString histOut){
  //
  _anaConf.printInfo();
  //assert(0);
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
  //Float_t NGB_rate_in_MHz = 386.0;
  //Float_t fadc_electronic_noise_RMS = 3.94;
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
  TH1D *h1_n_pe = new TH1D("h1_n_pe","h1_n_pe",1001,-0.5,1000.5);
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
  TH1D *h1_theta_pix_rot = new TH1D("h1_theta_pix_rot","h1_theta_pix_rot",4000,-0.1,2*TMath::Pi()+0.1);
  TH1D *h1_r_pix_rot = new TH1D("h1_r_pix_rot","h1_r_pix_rot",4000,-0.1,2.0);
  TH1D *h1_theta_deg_pix_rot = new TH1D("h1_theta_deg_pix_rot","h1_theta_deg_pix_rot",4000, 0.0, 360.0);
  //
  //
  unsigned int n_pe_steps = 50;
  //
  vector<TProfile*> pr_hmax_vs_e_v;
  tProfInit( pr_hmax_vs_e_v, n_pe_steps, "pr_hmax_vs_e_v", "pr_hmax_vs_e_v", 1000, 0, 10);
  //
  std::vector<TH1D*> h1_dx_v;
  h1D1Init(h1_dx_v, n_pe_steps, "h1_dx_v", "h1_dx_v", 100, 0.0, 2.0);
  std::vector<TH1D*> h1_dy_v;
  h1D1Init(h1_dy_v, n_pe_steps, "h1_dy_v", "h1_dy_v", 100, 0.0, 2.0);
  std::vector<TH1D*> h1_dt_v;
  h1D1Init(h1_dt_v, n_pe_steps, "h1_dt_v", "h1_dt_v", 100, 0.0, 75.0);
  std::vector<TH1D*> h1_tmean_v;
  h1D1Init(h1_tmean_v, n_pe_steps, "h1_tmean_v", "h1_tmean_v", 100, 0.0, 75.0);
  std::vector<TH1D*> h1_density_v;
  h1D1Init(h1_density_v, n_pe_steps, "h1_density_v", "h1_density_v", 1000, 0.0, 5000.0);
  std::vector<TH2D*> h2_dy_vs_dx_v;
  h2D2Init(h2_dy_vs_dx_v, n_pe_steps, "h2_dy_vs_dx_v", "h2_dy_vs_dx_v",
	   100, 0.0, 2.0, 100, 0.0, 2.0);
  std::vector<TH1D*> h1_x_std_v;
  h1D1Init(h1_x_std_v, n_pe_steps, "h1_x_std_v", "h1_x_std_v", 1000, 0, 1);
  std::vector<TH1D*> h1_y_std_v;
  h1D1Init(h1_y_std_v, n_pe_steps, "h1_y_std_v", "h1_y_std_v", 1000, 0, 1);
  //
  std::vector<sipmCameraHistCropped*> simp_hist_crop_v;
  TString sipm_hist_crop_name;
  //
  Double_t x0_LST01 = -70.93;
  Double_t y0_LST01 = -52.07;
  //
  Double_t r_core = 0.0;
  Double_t theta_core = 0.0;
  //
  Double_t x_pix_mean, y_pix_mean;
  Double_t x_pix_min, x_pix_max;
  Double_t y_pix_min, y_pix_max;
  Double_t dx_pix, dy_pix;
  Double_t x_pix_std, y_pix_std;
  //
  Double_t t_pix_min;
  Double_t t_pix_max;
  Double_t t_pix_mean;
  Double_t t_pix_std;
  Int_t dt_pix;
  Double_t V_phase;
  //
  TPrincipal *principal = NULL;
  //
  sipmCameraHistCropped *tmp_cam_hist = new sipmCameraHistCropped("tmp","tmp",sipm_cam,"sipmCameraHistCropped_pix.map");  
  //principal = new TPrincipal(tmp_cam_hist->get_n_pixels(),"princ");
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
    getCore_rel_R_theta( x0_LST01, y0_LST01, xcore, ycore, r_core, theta_core);
    _r_core = r_core;
    _theta_core = theta_core;
    if(cuts()){
      n_ev_cuts++;
      if(n_ev_cuts>_anaConf.n_ev_cuts_max && _anaConf.n_ev_cuts_max>0)
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
	Int_t jj_bad = 0;
	if(pe_chID[jj]>=0 && pe_chID[jj]<(Int_t)nn_PMT_channels){
	  sipm_cam->get_pixel_vec().at(pe_chID[jj]).pix_phi;
	  h1_pix_r->Fill(sipm_cam->get_pixel_vec().at(pe_chID[jj]).pix_r);
	  if(sipm_cam->get_pixel_vec().at(pe_chID[jj]).pix_r>0.2)
	    h2_theta_core_vs_pix_theta->Fill(theta_core, sipm_cam->get_pixel_vec().at(pe_chID[jj]).pix_phi);
	}
	else{
	  if(jj_bad == 0){
	    cout<<" ERROR  --> pe_chID[jj]<0 || pe_chID[jj]>=nn_PMT_channels"<<endl
	  	<<" jentry -->  "<<jentry<<endl
		<<"          jj "<<jj<<endl
	  	<<" pe_chID[jj] "<<pe_chID[jj]<<endl;
	  }
	  jj_bad++;
	}
      }
      if(n_pe<(Int_t)n_pe_steps){
	sipm_cam->get_pix_density_info( n_pe, pe_chID,
					x_pix_mean, y_pix_mean,
					x_pix_min, x_pix_max,
					y_pix_min, y_pix_max,
					dx_pix, dy_pix,
					x_pix_std, y_pix_std);
	sipm_cam->get_pix_time_info(n_pe, pe_time,
				    ev_time,
				    time_offset,
				    t_pix_min, t_pix_max,
				    t_pix_mean, t_pix_std,
				    dt_pix);
	pr_hmax_vs_e_v.at(n_pe)->Fill(energy,hmax);
	h1_dx_v.at(n_pe)->Fill(dx_pix);
	h1_dy_v.at(n_pe)->Fill(dy_pix);
	h1_dt_v.at(n_pe)->Fill(dt_pix);
	h1_tmean_v.at(n_pe)->Fill(t_pix_mean);
	h2_dy_vs_dx_v.at(n_pe)->Fill(dx_pix,dy_pix);
	h1_x_std_v.at(n_pe)->Fill(x_pix_std);
	h1_y_std_v.at(n_pe)->Fill(y_pix_std);
	V_phase = dx_pix*dy_pix*dt_pix;
	if(V_phase>0)
	  h1_density_v.at(n_pe)->Fill(n_pe/V_phase);
	//if(n_pe == 1 && dx_pix != 0){
	//cout<<"x_pix_min "<<x_pix_min<<endl
	//    <<"x_pix_max "<<x_pix_max<<endl;
	//}
      }
      //
      //if(cuts()){
      h2_ycore_vs_xcore->Fill(xcore,ycore);
      h2_ycore_vs_xcore_w->Fill(xcore,ycore,n_pe);
      //
      ///*
      if(_anaConf.wf_trg_sim){
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
	//
	//sipm_cam->Fill_pe_center(n_pe, pe_chID);
	//sipm_cam->Fill_pe(n_pe, pe_chID);
      }
      sipm_cam->Fill_pe(n_pe, pe_chID);
      //sipm_cam->Fill_pe(n_pe, pe_chID, -theta_core);
      //sipm_cam->Fill_pe(n_pe, pe_chID, -theta_core, h1_theta_pix_rot, h1_theta_deg_pix_rot, h1_r_pix_rot);
      //
      sipm_hist_crop_name = "sipm_cam_crop";
      sipm_hist_crop_name += "_ev";
      sipm_hist_crop_name += n_ev_cuts;
      sipmCameraHistCropped* simp_hist_crop_tmp = new sipmCameraHistCropped(sipm_hist_crop_name.Data(),
									    sipm_hist_crop_name.Data(),
									    sipm_cam,
									    tmp_cam_hist->get_pixel_map());
      simp_hist_crop_tmp->Fill_pe(n_pe, pe_chID, pe_time, ev_time, time_offset, -theta_core, principal);
      simp_hist_crop_v.push_back(simp_hist_crop_tmp);
      if(simp_hist_crop_v.size() == 1000){
	tmp_cam_hist->Save_to_csv("data_non_normalized.csv",simp_hist_crop_v);
	for(unsigned int iii = 0; iii<simp_hist_crop_v.size();iii++)
	  delete simp_hist_crop_v.at(iii);
	simp_hist_crop_v.clear();
      }
      //delete simp_hist_crop_tmp;
      //
      //sipm_cam->Fill(0,0);
      //}
      //std::vector<std::vector<unsigned int>> trg_vec;
      //*/
      //printEv();    
    }
  }
  //
  Int_t *pix_id_test_flower;
  Float_t *pe_time_test_flower;
  Int_t npixels_n_test_flower = 7;
  pix_id_test_flower = new Int_t[npixels_n_test_flower];
  pe_time_test_flower = new Float_t[npixels_n_test_flower];
  sipm_cam->simulateFlover_ideal_resp( 0, npixels_n_test_flower, pix_id_test_flower, pe_time_test_flower);
  sipm_cam->get_pix_density_info( npixels_n_test_flower, pix_id_test_flower,
				  x_pix_mean, y_pix_mean,
				  x_pix_min, x_pix_max,
				  y_pix_min, y_pix_max,
				  dx_pix, dy_pix,
				  x_pix_std, y_pix_std, 1);
  sipm_cam->get_pix_time_info(npixels_n_test_flower, pe_time_test_flower,
			      0,
			      0,
			      t_pix_min, t_pix_max,
			      t_pix_mean, t_pix_std,
			      dt_pix, 1);
  V_phase = dx_pix*dy_pix*dt_pix;
  if(V_phase>0)
    cout<<"npixels_n_test_flower/V_phase "<<npixels_n_test_flower/V_phase<<endl;  
  //
  //----------------------
  TH2D_divide( h2_ycore_vs_xcore_w, h2_ycore_vs_xcore, h2_ycore_vs_xcore_norm);
  //----------------------
  //
  // Do the actual analysis
  if(principal != NULL){
    principal->MakePrincipals();
    principal->Print();
  }
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
  tmp_cam_hist->Save_to_csv("data_non_normalized.csv",simp_hist_crop_v);
  //
  load_S_Vh_data("./S.cvs","./Vh.cvs");
  //
  std::vector<sipmCameraHistCropped*> sipm_cam_principal_hist_v;
  tmp_cam_hist->Fill_principal( sipm_cam_principal_hist_v, _data_Vh);
  //
  for(unsigned int ii = 0;ii<simp_hist_crop_v.size();ii++)
    simp_hist_crop_v.at(ii)->Write(); 
  //
  //
  TCanvas *c1 = NULL;
  TCanvas *c2 = NULL;
  c1 = new TCanvas("c1","c1",10,10,1500,1000);
  tmp_cam_hist->draw_crop_vector( 11, 3, sipm_cam_principal_hist_v, c1);
  //c2 = new TCanvas("c2","c2",10,10,1500,1000);
  //tmp_cam_hist->draw_crop_vector( 11, 3, simp_hist_crop_v, c2);
  //
  //
  //cout<<"sipm_cam_principal_hist_v.size() = "<<sipm_cam_principal_hist_v.size()<<endl;
  for(unsigned int ii = 0;ii<sipm_cam_principal_hist_v.size();ii++)
    sipm_cam_principal_hist_v.at(ii)->Write();  
  //
  //
  for(unsigned int ii = 0;ii<n_pe_steps;ii++){
    pr_hmax_vs_e_v.at(ii)->Write();
    h1_dx_v.at(ii)->Write();
    h1_dy_v.at(ii)->Write();
    h1_dt_v.at(ii)->Write();
    h1_tmean_v.at(ii)->Write();
    h1_density_v.at(ii)->Write();
    h2_dy_vs_dx_v.at(ii)->Write();
    h1_x_std_v.at(ii)->Write();
    h1_y_std_v.at(ii)->Write();
  }
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
  h1_theta_pix_rot->Write();
  h1_theta_deg_pix_rot->Write();
  h1_r_pix_rot->Write();
  //
  sipm_cam->Write();
  //
  if(c1 != NULL)
    c1->Write();
  if(c2 != NULL)
    c2->Write();
  //
  rootFile->Close();
  //
  //
  cout<<"n_ev_cuts = "<<n_ev_cuts<<endl;
}

bool anaPCA::cuts(){
  //
  if(_anaConf.disable_all_cuts)
    return true;
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
  if(hmax > 10000 && hmax < 11000)
    if(n_pe>1000 && n_pe<2000)
      return true;
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

void anaPCA::load_S_Vh_data(TString name_S, TString name_Vh){
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
