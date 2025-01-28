//C, C++
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <vector>

#include <time.h>

using namespace std;

void gen_ring(TGraph *gr, Int_t np, Double_t x0, Double_t y0, Double_t R);
void get_approximate_ring_parameters(TGraph *gr, TRandom3 *rnd, Double_t &x0, Double_t &y0, Double_t &R);
void get_Rvs_theta_and_theta_dist(TGraph *gr, Double_t x0, Double_t y0, Double_t R, TGraph *gr_R, TH1D *h1_theta_deg);

Int_t fit_muon_ring(){
  //
  TRandom3 *rnd = new TRandom3(12312312);
  //
  TString fileN01;
  fileN01 = "../scratch/simtel_data/muon/hist/hist_run1_muon.root";
  //
  TFile *f01 = new TFile(fileN01.Data());
  //
  TGraph *gr = (TGraph*)f01->Get("trueRing/gr_8944");
  TGraph *gr_frame = (TGraph*)f01->Get("gr_frame");
  //
  Double_t x0app, y0app, Rapp;
  Double_t x0app_average = 0.0;
  Double_t y0app_average = 0.0;
  Double_t Rapp_average = 0.0;
  vector<Double_t> x0app_v;
  vector<Double_t> y0app_v;
  vector<Double_t> Rapp_v;
  for(Int_t i = 0;i<50;i++){
    get_approximate_ring_parameters( gr, rnd, x0app, y0app, Rapp);
    x0app_v.push_back(x0app);
    y0app_v.push_back(y0app);
    Rapp_v.push_back(Rapp);
    x0app_average += x0app;
    y0app_average += y0app;
    Rapp_average += Rapp;
    //cout<<"x0app "<<x0app<<endl
    //    <<"y0app "<<y0app<<endl
    //    <<"Rapp  "<<Rapp<<endl;
  }
  //
  x0app_average /= x0app_v.size();
  y0app_average /= x0app_v.size();
  Rapp_average /= x0app_v.size();
  //
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
  TGraph *gr_app_r0 = new TGraph();
  gr_app_r0->SetNameTitle("gr_app_r0","gr_app_r0");
  for(unsigned int ii = 0;ii<x0app_v.size();ii++)
    gr_app_r0->SetPoint( ii, x0app_v.at(ii), y0app_v.at(ii));
  //
  gr->SetMarkerStyle(20);
  gr_frame->SetMarkerStyle(25);
  gr_app_r0->SetMarkerStyle(43);
  gr_app_r0->SetMarkerColor(kBlue);
  gr_app_r0->SetMarkerSize(2.0);
  //  
  TGraph *gr_R_app_average = new TGraph();
  gr_R_app_average->SetNameTitle("gr_R_app_average","gr_R_app_average");
  TH1D *h1_theta_app_average_deg = new TH1D("h1_theta_app_average_deg","h1_theta_app_average_deg", 36, 0.0, 360.0); 
  get_Rvs_theta_and_theta_dist( gr, x0app_average, y0app_average, Rapp_average, gr_R_app_average, h1_theta_app_average_deg);
  //
  TCanvas *c1 = new TCanvas("c1","c1",10,10,1800,600);
  c1->Divide(3,1);
  c1->cd(1);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetStatColor(kWhite);
  //gStyle->SetOptStat(kFALSE);
  TMultiGraph *mg = new TMultiGraph();
  mg->Add(gr);
  mg->Add(gr_frame);
  mg->Add(gr_app_r0);  
  mg->Add(gr_app_average_r0);
  mg->Add(gr_app_average_ring);
  mg->Draw("AP");
  //
  c1->cd(2);
  gr_R_app_average->SetMarkerStyle(20);
  gr_R_app_average->Draw("AP");  
  //
  c1->cd(3);
  h1_theta_app_average_deg->SetLineColor(kBlack);
  h1_theta_app_average_deg->SetLineWidth(2.0);
  h1_theta_app_average_deg->Draw();
  //
  return 0;
}

void gen_ring(TGraph *gr, Int_t np, Double_t x0, Double_t y0, Double_t R){
  Double_t phi = 0.0;
  TVector2 rc(x0,y0);
  for(Int_t i = 0;i<np;i++){
    TVector2 p;
    p.SetMagPhi(R,2*TMath::Pi()/(np-1)*i);
    TVector2 pt = rc + p;
    gr->SetPoint( i, pt.X(), pt.Y());
  }
}

void get_approximate_ring_parameters(TGraph *gr, TRandom3 *rnd, Double_t &x0, Double_t &y0, Double_t &R){
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

void get_Rvs_theta_and_theta_dist(TGraph *gr, Double_t x0, Double_t y0, Double_t R, TGraph *gr_R, TH1D *h1_theta_deg){
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
