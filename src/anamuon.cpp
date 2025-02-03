//my
#include "anamuon.hh"
#include "ana.hh"
#include "sipmCameraHist.hh"
#include "wfCamSim.hh"
#include "triggerSim.hh"
#include "dbscan.hh"
#include "muonRingFitter.hh"

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
#include <TRandom3.h>
#include <TDirectory.h>
#include <TMultiGraph.h>
#include <TText.h>
#include <TMinuit.h>

//C, C++
#include <iostream>
#include <stdlib.h>
#include <assert.h>
#include <fstream>
#include <iomanip>
#include <vector>
#include <time.h>
#include <bits/stdc++.h>
#include <sys/stat.h>

using namespace std;

void anamuon::Loop(TString histOut){
  //
  TRandom3 *rnd = new TRandom3(11231);
  //
  sipmCameraHist *sipm_cam = new sipmCameraHist("sipm_cam","sipm_cam","pixel_mapping.csv",0);
  TH1D* h1_chx = new TH1D("h1_chx","h1_chx",1000,-1.2,1.2);
  TH1D* h1_chy = new TH1D("h1_chy","h1_chy",1000,-1.2,1.2);
  TH2D* h2_chy_vs_chx = new TH2D("h2_chy_vs_chx","h2_chy_vs_chx",400,-1.2,1.2,400,-1.2,1.2);
  //
  TH1D* h1_pixID_tot = new TH1D("h1_pixID_tot","h1_pixID_tot",(Int_t)sipm_cam->_n_pixels,-0.5,(sipm_cam->_n_pixels-0.5));  
  //
  TH1D *h1_n_pe = new TH1D("h1_n_pe","h1_n_pe",10000,0.0,100000);
  TH1D *h1_nphotons = new TH1D("h1_nphotons","h1_nphotons",10000,0.0,100000);
  //
  TH1D *h1_h_first_int = new TH1D("h1_h_first_int","h1_h_first_int",1000,0.0,6000);
  TH1D *h1_hmax = new TH1D("h1_hmax","h1_hmax",1000,0.0,6000);
  TH2D *h2_hmax_vs_h_first_int = new TH2D("h2_hmax_vs_h_first_int","h2_hmax_vs_h_first_int",400,0.0,6000,400,0.0,6000);
  //
  TH1D *h1_xmax = new TH1D("h1_xmax","h1_xmax",1000,0.0,1000);
  TH1D *h1_emax = new TH1D("h1_emax","h1_emax",1000,-1000.0,1000);
  TH1D *h1_cmax = new TH1D("h1_cmax","h1_cmax",1000,0.0,1000);
  //
  TH2D *h2_n_pe_vs_nphotons = new TH2D("h2_n_pe_vs_nphotons","h2_n_pe_vs_nphotons",400,0.0,10000,400,0.0,6000);
  //  
  TH1D *h1_energy = new TH1D("h1_energy","h1_energy",10000,0.0,10);
  TH1D *h1_xcore = new TH1D("h1_xcore","h1_xcore",400,-15.0,15);
  TH1D *h1_ycore = new TH1D("h1_ycore","h1_ycore",400,-15.0,15);
  TH1D *h1_theta = new TH1D("h1_theta","h1_theta",1000,-1.0,10.0);  
  //
  TH1D *h1_azimuth_deg = new TH1D("h1_azimuth_deg","h1_azimuth_deg",400,180.0-10.0,180+10.0);
  TH1D *h1_altitude_deg = new TH1D("h1_altitude_deg","h1_altitude_deg",400,80.0-2.0,80.0+2.0);
  //
  TH2D *h2_n_pe_vs_energy = new TH2D("h2_n_pe_vs_energy","h2_n_pe_vs_energy", 2000, 0.0, 0.2, 100, 0.0, 5000.0);
  TH2D *h2_ycore_vs_xcore = new TH2D("h2_ycore_vs_xcore","h2_ycore_vs_xcore", 400, -12.0, 12.0, 400, -12.0, 12.0);
  //
  TH1D *h1_muonx0 = new TH1D("h1_muonx0","h1_muonx0",400, -2.0,2.0);
  TH1D *h1_muony0 = new TH1D("h1_muony0","h1_muony0",400, -2.0,2.0);
  TH1D *h1_muon_r0 = new TH1D("h1_muon_r0","h1_muon_r0",400, 0.0,2.0);
  TH1D *h1_muonR = new TH1D("h1_muonR","h1_muonR",400, 0.0,2.0);
  //
  TH2D *h2_muon_r0_vs_thetaDeg = new TH2D("h2_muon_r0_vs_thetaDeg","h2_muon_r0_vs_thetaDeg", 400, 0.0, 1.0, 400, 0.0, 1.0);
  TH2D *h2_muonR_vs_muon_E = new TH2D("h2_muonR_vs_muon_E","h2_muonR_vs_muon_E",400, 0.0, 0.2, 400, 0.0, 2.0);
  //
  Double_t thetaDeg;
  Double_t coreR;
  Double_t azimuth_deg;
  Double_t altitude_deg;
  //
  TGraph *gr_frame = new TGraph();
  gr_frame->SetNameTitle("gr_frame","gr_frame");
  gr_frame->SetPoint(0, -1.2, -1.2);
  gr_frame->SetPoint(1, -1.2,  1.2);
  gr_frame->SetPoint(2,  1.2,  1.2);
  gr_frame->SetPoint(3,  1.2, -1.2);
  gr_frame->SetPoint(4, -1.2, -1.2);  
  //
  vector <TGraph*> gr_v;
  vector <TGraph*> gr_tmp_dbscan_clean_v;
  vector <TCanvas*> c1_v;
  //
  TGraph *gr_n_pe_vs_jentry = new TGraph();
  gr_n_pe_vs_jentry->SetNameTitle("gr_n_pe_vs_jentry","gr_n_pe_vs_jentry");
  //
  //
  Double_t muonx0;
  Double_t muony0;
  Double_t muon_r0;
  Double_t muonR;
  //
  Double_t x0out_fit, y0out_fit, Rout_fit;
  Double_t x0out_fit_err, y0out_fit_err, Rout_fit_err;
  //
  Long64_t nentries = fChain->GetEntriesFast();
  cout<<"nentries = "<<nentries<<endl;
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    if(jentry%100 == 0)
      cout<<jentry<<endl;
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    //if (ientry > 10) break;
    nb = fChain->GetEntry(jentry); nbytes += nb;
    //
    thetaDeg = get_theta_p_t()*180.0/TMath::Pi();
    coreR = TMath::Sqrt(xcore*xcore + ycore*ycore);
    azimuth_deg = azimuth*180.0/TMath::Pi();
    altitude_deg =  altitude*180.0/TMath::Pi();
    //
    //    
    //if(h_first_int<10){
    //if(coreR<5){
    //if(cmax<791){
    //if(n_pe<2200 && energy<0.016){
    if(energy>0.0130 && energy<0.014){
      //if(n_pe>1500 && n_pe<3000){
    //if(n_pe>3000){
    //if(thetaDeg>0.35 && thetaDeg<0.45){
    //if(azimuth_deg>179.0 && azimuth_deg<181.0){
    //if(energy>0.013){
    //    
    for(Int_t kk = 1;kk<=(Int_t)sipm_cam->_n_pixels;kk++)
      h1_pixID_tot->SetBinContent(kk,0.0);
    for(Int_t kk = 0;kk<n_pe;kk++)
      h1_pixID_tot->Fill(pe_chID[kk]);
    //
    h1_azimuth_deg->Fill(azimuth_deg);
    h1_altitude_deg->Fill(altitude_deg);
    //
    gr_n_pe_vs_jentry->SetPoint(gr_n_pe_vs_jentry->GetN(),(Int_t)jentry,n_pe);
    h1_n_pe->Fill(n_pe);
    h1_nphotons->Fill(nphotons);
    //
    h1_h_first_int->Fill(h_first_int);
    h1_hmax->Fill(hmax);
    h2_hmax_vs_h_first_int->Fill( h_first_int, hmax);
    //
    h1_xmax->Fill(xmax);
    h1_emax->Fill(emax);
    h1_cmax->Fill(cmax);
    //
    h2_n_pe_vs_nphotons->Fill( nphotons, n_pe);
    h1_energy->Fill(energy);
    h2_n_pe_vs_energy->Fill(energy,n_pe);
    h1_xcore->Fill(xcore);
    h1_ycore->Fill(ycore);
    h1_theta->Fill(thetaDeg);
    //
    h2_ycore_vs_xcore->Fill( xcore, ycore);
    //
    sipm_cam->Fill_pix_x_y_hist( n_pe, pe_chID, h1_chx, h1_chy);
    sipm_cam->Fill_pix_hist2D_y_vs_x( n_pe, pe_chID, h2_chy_vs_chx);
    //
    //----------------------
    //
    TString h1nametitle = "h1_pne_per_pix_";
    h1nametitle += (Int_t)jentry;
    TH1D *h1_pne_per_pix = new TH1D(h1nametitle.Data(),h1nametitle.Data(),101,-0.5,100.5);
    for(Int_t kk = 1;kk<=(Int_t)sipm_cam->_n_pixels;kk++)
      h1_pne_per_pix->Fill(h1_pixID_tot->GetBinContent(kk));
    //    
    //----------------------
    Int_t n_pe_tmp = 0;
    Int_t *pe_chID_tmp = new Int_t[(Int_t)sipm_cam->_n_pixels];
    for(Int_t kk = 1;kk<=(Int_t)sipm_cam->_n_pixels;kk++){
      if(h1_pixID_tot->GetBinContent(kk)>1){
	pe_chID_tmp[n_pe_tmp] = (Int_t)h1_pixID_tot->GetBinCenter(kk);
	n_pe_tmp++;
      }
    }
    //
    TGraph *gr_tmp = new TGraph();
    TString nametitle="gr_";
    nametitle += (Int_t)jentry;
    gr_tmp->SetNameTitle(nametitle,nametitle);
    //sipm_cam->Fill_pix_TGraph_y_vs_x( n_pe, pe_chID, gr_tmp);
    sipm_cam->Fill_pix_TGraph_y_vs_x( n_pe_tmp, pe_chID_tmp, gr_tmp);
    //
    if(gr_tmp->GetN()>20){
      //---------------------------------------------------------
      //
      TGraph *gr_tmp_dbscan_clean = new TGraph();
      TString nametitle="gr_tmp_dbscan_clean_";
      nametitle += (Int_t)jentry;
      gr_tmp_dbscan_clean->SetNameTitle(nametitle,nametitle);
      dbscan_cleaning(gr_tmp,gr_tmp_dbscan_clean, 10, 0.1);
      //sipm_cam->Fill_pix_TGraph_y_vs_x( n_pe, pe_chID, gr_tmp);
      sipm_cam->Fill_pix_TGraph_y_vs_x( n_pe_tmp, pe_chID_tmp, gr_tmp);
      //
      //---------------------------------------------------------      
      gr_v.push_back(gr_tmp);      
      gr_tmp_dbscan_clean_v.push_back(gr_tmp_dbscan_clean);
      //get_average_approximate_ring_parameters(gr_tmp, rnd, 100, muonx0, muony0, muonR);  
      get_average_approximate_ring_parameters(gr_tmp_dbscan_clean, rnd, 100, muonx0, muony0, muonR);
      //
      //fit
      //cout<<"muonRingFitter *mrf "<<endl
      //<<"gr_tmp->GetN()      "<<gr_tmp->GetN()<<endl;
      //
      muonRingFitter *mrf = new muonRingFitter();
      //mrf->set_gr(gr_tmp_dbscan_clean);
      mrf->set_gr(gr_tmp);

      mrf->fit_ring(muonx0, muony0, muonR,
		    x0out_fit, y0out_fit, Rout_fit,
		    x0out_fit_err, y0out_fit_err, Rout_fit_err);
      delete mrf;
      //
      //
      h1_muonx0->Fill(muonx0);
      h1_muony0->Fill(muony0);
      muon_r0 = TMath::Sqrt(muonx0*muonx0 + muony0*muony0);
      h1_muon_r0->Fill(muon_r0);
      h1_muonR->Fill(muonR);
      //
      h2_muon_r0_vs_thetaDeg->Fill(thetaDeg,muon_r0);
      h2_muonR_vs_muon_E->Fill(energy,muonR);
      //
      TString c1canvaname="c1_";
      c1canvaname += (Int_t)jentry;
      TCanvas *c1 = new TCanvas(c1canvaname.Data(),c1canvaname.Data(),10,10,1500,1000);
      fit_muon_ring(gr_tmp, gr_frame, gr_tmp_dbscan_clean, c1, (Int_t)jentry,
		    muonx0, muony0, muonR,
		    x0out_fit, y0out_fit, Rout_fit,
		    (Int_t)jentry,
		    energy,
		    xcore,
		    ycore,
		    ev_time,
		    nphotons,
		    n_pe,
		    n_pixels,
		    azimuth,
		    altitude,
		    h_first_int,
		    hmax,
		    thetaDeg,
		    h1_pne_per_pix,
		    "muon",
		    "NONE");
      //
      c1_v.push_back(c1);
    }//if(gr_tmp->GetN()>20){

    }
    //----------------------
  }
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
  TDirectory *trueRing = rootFile->mkdir("trueRing");
  trueRing->cd();
  for(unsigned int i = 0;i<gr_v.size();i++){
    gr_v.at(i)->Write();
    c1_v.at(i)->Write();
  }
  rootFile->cd();
  //
  gr_frame->Write();
  //
  h1_n_pe->Write();
  h1_nphotons->Write();
  h2_n_pe_vs_nphotons->Write();
  h1_energy->Write();
  gr_n_pe_vs_jentry->Write();
  h2_n_pe_vs_energy->Write();
  h1_xcore->Write();
  h1_ycore->Write();
  h1_theta->Write();
  h1_azimuth_deg->Write();
  h1_altitude_deg->Write();
  h2_ycore_vs_xcore->Write();
  //
  h1_chx->Write();
  h1_chy->Write();
  h2_chy_vs_chx->Write();
  //
  h1_h_first_int->Write();
  h1_hmax->Write();
  h2_hmax_vs_h_first_int->Write();
  //
  h1_xmax->Write();
  h1_emax->Write();
  h1_cmax->Write();
  //
  h1_muonx0->Write();
  h1_muony0->Write();
  h1_muon_r0->Write();
  h1_muonR->Write();
  //
  h2_muon_r0_vs_thetaDeg->Write();
  h2_muonR_vs_muon_E->Write();
  //
  h1_pixID_tot->Write();
  //
  rootFile->Close();
}

Double_t anamuon::get_canonical_phi_from_V2phi(Double_t phiv2){
  Double_t phiv2_deg = phiv2*180.0/TMath::Pi();
  if(phiv2_deg > 360.0 || phiv2_deg<0.0){
    cout<<" ERROR --> phiv2_deg > 360.0 || phiv2_deg < 0.0 "<<endl
        <<"                                phiv2_deg = "<<phiv2_deg<<endl;
  }
  if( phiv2_deg >= 0.0 && phiv2_deg <= 180.0){
    return phiv2;
  }  
  else{
    return phiv2 - TMath::TwoPi();
  }
}

Double_t anamuon::get_theta_p_t(){
  TVector3 v_det(1.0*TMath::Sin(10.0/180.0*TMath::Pi()),
		 0,
		 1.0*TMath::Cos(10.0/180.0*TMath::Pi()));
  TVector3 v_prot;
  v_prot.SetMagThetaPhi(1.0,TMath::Pi()/2.0-altitude,TMath::Pi() - azimuth);
  TVector3 v_prot_inv(v_prot.x(),v_prot.y(),v_prot.z());
  return TMath::ACos(v_prot_inv.Dot(v_det)/v_prot_inv.Mag()/v_det.Mag());
}

void anamuon::gen_azimuth_line(TGraph *gr, Int_t np, Double_t x0, Double_t y0, Double_t azimuth_deg){
  TVector2 rc(x0,y0);
  Double_t Rmax = 5.0;
  for(Int_t i = 0;i<np;i++){
    TVector2 p;
    p.SetMagPhi(Rmax/(np-1)*i,azimuth_deg/180.0*TMath::Pi());
    TVector2 pt = rc + p;
    gr->SetPoint( i, pt.X(), pt.Y());
  }
}

void anamuon::gen_altitude_line(TGraph *gr, Int_t np, Double_t x0, Double_t y0, Double_t azimuth_val, Double_t altitude_val){
  TVector3 rc(x0, y0, 0.0);
  Double_t Zmax = 600.0;
  Double_t Rmax = 700.0;
  for(Int_t i = 0;i<np;i++){
    TVector3 p;
    p.SetMagThetaPhi(Rmax/(np-1)*i, TMath::Pi()/2.0-altitude_val, get_canonical_phi_from_V2phi(azimuth_val));
    TVector3 pt = rc + p;
    if(pt.Z()<Zmax)
      gr->SetPoint( i, pt.X(), pt.Z());
  }  
}

void anamuon::gen_ring(TGraph *gr, Int_t np, Double_t x0, Double_t y0, Double_t R){
  Double_t phi = 0.0;
  TVector2 rc(x0,y0);
  for(Int_t i = 0;i<np;i++){
    TVector2 p;
    p.SetMagPhi(R,2*TMath::Pi()/(np-1)*i);
    TVector2 pt = rc + p;
    gr->SetPoint( i, pt.X(), pt.Y());
  }
}

void anamuon::get_approximate_ring_parameters(TGraph *gr, TRandom3 *rnd, Double_t &x0, Double_t &y0, Double_t &R){
  Double_t p1x0, p1y0;
  Double_t p2x0, p2y0; 
  Double_t maxDist = 0.0;
  Int_t p1_id = (Int_t)rnd->Uniform(0,gr->GetN());
  Int_t p2_id = (Int_t)rnd->Uniform(0,gr->GetN());
  gr->GetPoint(p1_id, p1x0, p1y0);
  TVector2 p1(p1x0, p1y0);
  gr->GetPoint(p2_id, p2x0, p2y0);
  TVector2 p2(p2x0, p2y0);
  for(Int_t i = 0;i<gr->GetN();i++){
    gr->GetPoint(i, p2x0, p2y0);
    p2.Set(p2x0, p2y0);
    TVector2 dp = p1 - p2;
    if(maxDist<dp.Mod()){
      maxDist = dp.Mod();
      p2_id = i;
    }
  }
  //
  gr->GetPoint(p2_id, p2x0, p2y0);
  p2.Set(p2x0, p2y0);
  //  
  TVector2 pm = (p1 + p2)/2.0;
  x0 = pm.X();
  y0 = pm.Y();
  R = maxDist/2.0;
}

void anamuon::get_Rvs_theta_and_theta_dist(TGraph *gr, Double_t x0, Double_t y0, Double_t R, TGraph *gr_R, TH1D *h1_theta_deg){
  Double_t x, y;
  TVector2 pc(x0, y0);
  for(Int_t i = 0;i<gr->GetN();i++){
    gr->GetPoint(i, x, y);
    TVector2 p(x,y);
    TVector2 dp = p - pc;
    gr_R->SetPoint(i,dp.Phi()*180.0/TMath::Pi(),(dp.Mod() - R)/R);
    h1_theta_deg->Fill(dp.Phi()*180.0/TMath::Pi());
  }
}

void anamuon::get_average_approximate_ring_parameters(TGraph *gr, TRandom3 *rnd, Int_t niterations, Double_t &x0, Double_t &y0, Double_t &R){
  Double_t x0app, y0app, Rapp;
  Double_t x0app_average = 0.0;
  Double_t y0app_average = 0.0;
  Double_t Rapp_average = 0.0;
  vector<Double_t> x0app_v;
  vector<Double_t> y0app_v;
  vector<Double_t> Rapp_v;
  for(Int_t i = 0;i<niterations;i++){
    get_approximate_ring_parameters( gr, rnd, x0app, y0app, Rapp);
    x0app_v.push_back(x0app);
    y0app_v.push_back(y0app);
    Rapp_v.push_back(Rapp);
    x0app_average += x0app;
    y0app_average += y0app;
    Rapp_average += Rapp;
  }
  //
  x0app_average /= x0app_v.size();
  y0app_average /= x0app_v.size();
  Rapp_average /= x0app_v.size();
  //
  x0 = x0app_average;
  y0 = y0app_average;
  R = Rapp_average;
}

void anamuon::fit_muon_ring(TGraph *gr, TGraph *gr_frame, TGraph *gr_clean, TCanvas *c1, Int_t jentry_int,
			    Double_t x0app_average, Double_t y0app_average, Double_t Rapp_average,
			    Double_t x0out_fit, Double_t y0out_fit, Double_t Rout_fit,
			    Int_t event_id_val,
			    Double_t energy_val,
			    Double_t xcore_val,
			    Double_t ycore_val,
			    Double_t ev_time_val,
			    Int_t  nphotons_val,
			    Int_t n_pe_val,
			    Int_t n_pixels_val,
			    Double_t azimuth_val,
			    Double_t altitude_val,
			    Double_t h_first_int_val,
			    Double_t hmax_val,
			    Double_t thetaDeg_val,
			    TH1D *h1_pne_per_pix,
			    TString particle_type,
			    TString pdf_out_file){
  //
  TGraph *gr_app_average_r0 = new TGraph();
  gr_app_average_r0->SetNameTitle("gr_app_average_r0","gr_app_average_r0");  
  gr_app_average_r0->SetPoint( 0, x0app_average, y0app_average);
  gr_app_average_r0->SetMarkerStyle(43);
  gr_app_average_r0->SetMarkerColor(kMagenta+3);
  gr_app_average_r0->SetMarkerSize(3.0);
  //
  TGraph *gr_app_average_ring = new TGraph();
  gr_app_average_ring->SetNameTitle("gr_app_average_ring","gr_app_average_ring");  
  gr_app_average_ring->SetMarkerStyle(7);
  gr_app_average_ring->SetMarkerColor(kMagenta+3);
  gr_app_average_ring->SetMarkerSize(3.0);
  gen_ring(gr_app_average_ring, 360, x0app_average, y0app_average, Rapp_average);
  //
  TGraph *gr_fit_ring = new TGraph();
  gr_fit_ring->SetNameTitle("gr_fit_ring","gr_fit_ring");  
  gr_fit_ring->SetMarkerStyle(1);
  gr_fit_ring->SetMarkerColor(kBlue+3);
  gr_fit_ring->SetMarkerSize(1.0);
  gen_ring(gr_fit_ring, 360, x0out_fit, y0out_fit, Rout_fit);
  //
  TGraph *gr_fit_ring_c = new TGraph();
  gr_fit_ring_c->SetNameTitle("gr_fit_ring_c","gr_fit_ring_c");  
  gr_fit_ring_c->SetMarkerStyle(43);
  gr_fit_ring_c->SetMarkerColor(kBlue+3);
  gr_fit_ring_c->SetMarkerSize(1.0);
  gr_fit_ring_c->SetPoint( 0, x0out_fit, y0out_fit);
  //
  TGraph *gr_R_app_average = new TGraph();
  TString name_R_app_average = "gr_R_app_average";
  name_R_app_average += jentry_int;
  gr_R_app_average->SetNameTitle(name_R_app_average.Data(),name_R_app_average.Data());
  TString name_theta_app_average_deg = "h1_theta_app_average_deg";
  name_theta_app_average_deg += jentry_int;
  TH1D *h1_theta_app_average_deg = new TH1D(name_theta_app_average_deg.Data(),name_theta_app_average_deg.Data(), 36, 0.0, 360.0); 
  get_Rvs_theta_and_theta_dist( gr, x0app_average, y0app_average, Rapp_average, gr_R_app_average, h1_theta_app_average_deg);
  //
  //
  TGraph *gr_generation_area = new TGraph();
  gr_generation_area->SetNameTitle("gr_generation_area","gr_generation_area");
  gen_ring(gr_generation_area, 400, 0.0, 0.0, 10.0);  
  //
  TGraph *gr_dish_area = new TGraph();
  gr_dish_area->SetNameTitle("gr_dish_area","gr_dish_area");
  gr_dish_area->SetMarkerColor(kMagenta);
  gen_ring(gr_dish_area, 400, 0.0, 0.0, 23.0/2.0);
  //
  TGraph *gr_core_pos_frame = new TGraph();
  Double_t core_pos_frame_L = 26.0;
  gr_core_pos_frame->SetNameTitle("gr_core_pos_frame","gr_core_pos_frame");
  gr_core_pos_frame->SetPoint(0, -core_pos_frame_L/2.0, -core_pos_frame_L/2.0);
  gr_core_pos_frame->SetPoint(1,  core_pos_frame_L/2.0, -core_pos_frame_L/2.0);
  gr_core_pos_frame->SetPoint(2,  core_pos_frame_L/2.0,  core_pos_frame_L/2.0);
  gr_core_pos_frame->SetPoint(3,  core_pos_frame_L/2.0, -core_pos_frame_L/2.0);
  gr_core_pos_frame->SetPoint(4, -core_pos_frame_L/2.0, -core_pos_frame_L/2.0);  
  gr_core_pos_frame->SetDrawOption("APL");
  TGraph *gr_core_pos = new TGraph();
  gr_core_pos->SetNameTitle("gr_core_pos","gr_core_pos");
  gr_core_pos->SetPoint(0,xcore_val,ycore_val);
  gr_core_pos->SetMarkerColor(kRed);
  gr_core_pos->SetMarkerStyle(34);
  gr_core_pos->SetMarkerSize(2.0);
  TGraph *gr_azimuth_line = new TGraph();
  gr_azimuth_line->SetNameTitle("gr_azimuth_line","gr_azimuth_line");
  gen_azimuth_line(gr_azimuth_line, 40, xcore_val, ycore_val, azimuth_val*180.0/TMath::Pi());
  //
  //
  TGraph *gr_altitude_line_frame = new TGraph();
  gr_altitude_line_frame->SetNameTitle("gr_altitude_line_frame","gr_altitude_line_frame");
  Double_t altitude_line_frame_H_min = -20.0;
  Double_t altitude_line_frame_H_max = 600.0;
  Double_t altitude_line_frame_L_min = -23.0/2.0;
  Double_t altitude_line_frame_L_max =  23.0/2.0;
  gr_altitude_line_frame->SetPoint(0,altitude_line_frame_L_min,altitude_line_frame_H_min);
  gr_altitude_line_frame->SetPoint(1,altitude_line_frame_L_min,altitude_line_frame_H_max);
  gr_altitude_line_frame->SetPoint(2,altitude_line_frame_L_max,altitude_line_frame_H_max);
  gr_altitude_line_frame->SetPoint(3,altitude_line_frame_L_max,altitude_line_frame_H_min);
  gr_altitude_line_frame->SetPoint(4,altitude_line_frame_L_min,altitude_line_frame_H_min);
  gr_altitude_line_frame->SetDrawOption("APL");
  TGraph *gr_opt_ax_line = new TGraph();
  gr_opt_ax_line->SetNameTitle("gr_opt_ax_line","gr_opt_ax_line");
  gen_altitude_line(gr_opt_ax_line, 100, 0.0, 0.0, TMath::Pi(), 80.0/180.0*TMath::Pi());
  TGraph *gr_altitude_line = new TGraph();
  gr_altitude_line->SetNameTitle("gr_altitude_line","gr_altitude_line");
  gen_altitude_line(gr_altitude_line, 400, xcore_val, ycore_val, azimuth_val, altitude_val);
  //
  //TCanvas *c1 = new TCanvas("c1","c1",10,10,1800,600);
  c1->Divide(3,2);
  c1->cd(1);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetStatColor(kWhite);
  gStyle->SetOptStat(kFALSE);
  //
  //
  gr->SetMarkerStyle(8);
  gr->SetMarkerSize(0.5);
  gr->SetMarkerColorAlpha(kBlack, 1.0);
  gr_clean->SetMarkerStyle(8);
  gr_clean->SetMarkerSize(0.5);
  gr_clean->SetMarkerColorAlpha(kRed, 1.0);
  //
  //
  TMultiGraph *mg = new TMultiGraph();
  mg->Add(gr);
  mg->Add(gr_clean);
  mg->Add(gr_frame);
  mg->Add(gr_app_average_r0);
  mg->Add(gr_app_average_ring);
  mg->Add(gr_fit_ring);
  mg->Add(gr_fit_ring_c);
  mg->Draw("AP");
  //
  c1->cd(2);
  if(h1_pne_per_pix == NULL){
    gr_R_app_average->SetTitle("");
    gr_R_app_average->SetMarkerStyle(20);
    gr_R_app_average->Draw("AP");  
  }
  else{
    h1_pne_per_pix->SetMinimum(0.5);
    h1_pne_per_pix->SetStats(kFALSE);
    h1_pne_per_pix->SetTitle("");
    h1_pne_per_pix->SetLineColor(kBlack);
    h1_pne_per_pix->SetLineWidth(2.0);
    h1_pne_per_pix->Draw();
    gPad->SetLogy();
  }
  //
  c1->cd(3);
  h1_theta_app_average_deg->SetStats(kFALSE);
  h1_theta_app_average_deg->SetTitle("");
  h1_theta_app_average_deg->SetLineColor(kBlack);
  h1_theta_app_average_deg->SetLineWidth(2.0);
  h1_theta_app_average_deg->Draw();
  //
  //
  //
  //
  c1->cd(4);
  TMultiGraph *mg02 = new TMultiGraph();
  mg02->Add(gr_core_pos_frame);
  mg02->Add(gr_generation_area);
  mg02->Add(gr_dish_area);
  mg02->Add(gr_core_pos);
  mg02->Add(gr_azimuth_line);
  mg02->Draw("AP");
  //
  //
  //
  //
  c1->cd(5);
  gr_altitude_line->SetMarkerStyle(9);
  gr_altitude_line->SetMarkerSize(2);
  TMultiGraph *mg03 = new TMultiGraph();
  mg03->Add(gr_altitude_line_frame);
  mg03->Add(gr_altitude_line);
  mg03->Add(gr_opt_ax_line);
  mg03->Draw("AP");
  //
  //
  c1->cd(6);
  ////////////////////////////////
  TString event_id_str = "event_id : "; event_id_str += event_id_val;
  TString energy_str   = "energy   : "; energy_str += (Int_t)(energy_val*1000.0); energy_str += " GeV";
  TString xcore_str    = "xcore    : "; xcore_str += (Int_t)xcore_val; xcore_str += " m";
  TString ycore_str    = "ycore    : "; ycore_str += (Int_t)ycore_val; ycore_str += " m";
  TString ev_time_str  = "ev_time  : "; ev_time_str += (Int_t)ev_time_val; ev_time_str += " ns";
  TString nphotons_str = "nphotons : "; nphotons_str += nphotons_val;
  TString n_pe_str     = "n_pe     : "; n_pe_str += n_pe_val;
  TString n_pixels_str = "n_pixels : "; n_pixels_str += n_pixels_val;
  //
  TString thetaDeg_str;
  TString azimuth_str;
  TString altitude_str;
  TString h_first_int_str;
  TString hmax_str;
  //
  thetaDeg_str	  = "thetaDeg    : "; thetaDeg_str += (Int_t)(thetaDeg_val*100); thetaDeg_str += "/100 deg";
  azimuth_str     = "azimuth     : "; azimuth_str += (Int_t)(azimuth_val*180.0/TMath::Pi()*10); azimuth_str += "/10 deg";
  altitude_str    = "altitude    : "; altitude_str += (Int_t)(altitude_val*180.0/TMath::Pi()*10); altitude_str += "/10 deg";
  h_first_int_str = "h_first_int : "; h_first_int_str += (Int_t)h_first_int_val*1000; h_first_int_str += " m";
  hmax_str        = "hmax        : "; hmax_str += (Int_t)hmax_val*1000; hmax_str += " m";
  //
  TString n_graph_points_str = "n_graph_points : "; n_graph_points_str += (Int_t)gr->GetN();
  //
  TText *t = new TText(0.5,0.9,particle_type.Data());
  t->SetTextAlign(22);
  t->SetTextFont(43);
  t->SetTextSize(40);
  t->Draw();
  TText *t2 = new TText(0.5,0.85,event_id_str.Data());
  t2->SetTextAlign(22);
  t2->SetTextFont(43);
  t2->SetTextSize(40);
  t2->Draw("same");
  TText *t3 = new TText(0.5,0.80,energy_str.Data());
  t3->SetTextAlign(22);
  t3->SetTextFont(43);
  t3->SetTextSize(40);
  t3->Draw("same");
  TText *t4 = new TText(0.5,0.75,xcore_str.Data());
  t4->SetTextAlign(22);
  t4->SetTextFont(43);
  t4->SetTextSize(40);
  t4->Draw("same");
  TText *t5 = new TText(0.5,0.70,ycore_str.Data());
  t5->SetTextAlign(22);
  t5->SetTextFont(43);
  t5->SetTextSize(40);
  t5->Draw("same");
  TText *t6 = new TText(0.5,0.65,ev_time_str.Data());
  t6->SetTextAlign(22);
  t6->SetTextFont(43);
  t6->SetTextSize(40);
  t6->Draw("same");
  TText *t7 = new TText(0.5,0.60,nphotons_str.Data());
  t7->SetTextAlign(22);
  t7->SetTextFont(43);
  t7->SetTextSize(40);
  t7->Draw("same");
  TText *t8 = new TText(0.5,0.55,n_pe_str.Data());
  t8->SetTextAlign(22);
  t8->SetTextFont(43);
  t8->SetTextSize(40);
  t8->Draw("same");
  TText *t9 = new TText(0.5,0.50,n_pixels_str.Data());
  t9->SetTextAlign(22);
  t9->SetTextFont(43);
  t9->SetTextSize(40);
  t9->Draw("same");
  //
  //
  TText *t10 = new TText(0.5,0.45,azimuth_str.Data());
  t10->SetTextAlign(22);
  t10->SetTextFont(43);
  t10->SetTextSize(40);
  t10->Draw("same");
  TText *t11 = new TText(0.5,0.40,altitude_str.Data());
  t11->SetTextAlign(22);
  t11->SetTextFont(43);
  t11->SetTextSize(40);
  t11->Draw("same");
  TText *t12 = new TText(0.5,0.35,h_first_int_str.Data());
  t12->SetTextAlign(22);
  t12->SetTextFont(43);
  t12->SetTextSize(40);
  t12->Draw("same");
  TText *t13 = new TText(0.5,0.30,hmax_str.Data());
  t13->SetTextAlign(22);
  t13->SetTextFont(43);
  t13->SetTextSize(40);
  t13->Draw("same");
  TText *t14 = new TText(0.5,0.25,thetaDeg_str.Data());
  t14->SetTextAlign(22);
  t14->SetTextFont(43);
  t14->SetTextSize(40);
  t14->Draw("same");
  TText *t15 = new TText(0.5,0.20,n_graph_points_str.Data());
  t15->SetTextAlign(22);
  t15->SetTextFont(43);
  t15->SetTextSize(40);
  t15->Draw("same");
  //
  if(pdf_out_file != "NONE")
    c1->SaveAs(pdf_out_file.Data());
  //  
  ////////////////////////////////
}

void anamuon::dbscan_cleaning( TGraph *gr, TGraph *gr_db, unsigned int minPts, Double_t eps){
  //
  _dbs->run( minPts, eps, gr);
  //_dbs->get_cluster_stats();
  //_dbs->print_cluster_stats();
  //
  for(unsigned int ii = 0; ii<_dbs->get_points_v().size(); ii++)
    if(_dbs->get_points_v().at(ii).point_type>0)
      gr_db->SetPoint(gr_db->GetN(),_dbs->get_points_v().at(ii).x,_dbs->get_points_v().at(ii).y);
  //
  _dbs->clear();
}
