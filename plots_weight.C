//C, C++
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

#include <time.h>

using namespace std;

void read_file(TString name, TH1D *h1, TH1D *h1_w, TH1D *h1_wo, TGraphErrors *gr);

Double_t power_law(Double_t *x, Double_t *par){
  //
  Double_t A = par[0];
  Double_t p = par[1];
  return A/TMath::Power(x[0],p);
}

Double_t weight_f(Double_t *x, Double_t *par){
  //
  Double_t A = par[0];
  Double_t p = par[1];
  return A/TMath::Power(x[0],p);
}

Int_t plots_weight(){
  //
  Double_t Emin = 0.006;
  Double_t Emax = 100;
  TH1D *h1_weight_01 = new TH1D("h1_weight_01","h1_weight_01",1000000,0.0,101);
  //
  TH1D *h1_E = new TH1D("h1_E","h1_E",500,0.0,100);
  TH1D *h1_Ew = new TH1D("h1_Ew","h1_Ew",500,0.0,100);
  TGraphErrors *gr_w_vs_E = new TGraphErrors();
  gr_w_vs_E->SetNameTitle("gr_w_vs_E","p0/E^p1");
  //
  //
  read_file("./weight.dat", h1_weight_01, h1_Ew, h1_E, gr_w_vs_E);
  //
  TCanvas *c1 = new TCanvas("c1","c1",10,10,600,600);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetStatColor(kWhite);
  //gStyle->SetOptStat(kFALSE);
  //
  TF1 *f_power_law = new TF1("power_law", power_law, Emin, Emax, 2);
  //
  //f_power_law->FixParameter(0, 100);
  //f_power_law->FixParameter(1, 2);
  f_power_law->SetParameter(0, 8000.0*0.011*0.011);
  f_power_law->SetParameter(1, 2);
  //
  TF1 *f_weight = new TF1("weight_f", weight_f, Emin, Emax, 2);
  //
  f_weight->SetParameter(0, 1.0);
  f_weight->SetParameter(1, 0.7); 
  //
  h1_weight_01->SetLineColor(kBlack);
  h1_weight_01->SetLineWidth(3.0);
  //
  //h1_weight_01->Fit(f_power_law, "",  "", Emin, Emax);
  //
  //h1_weight_01->Draw();
  //
  h1_Ew->SetLineColor(kRed);
  h1_Ew->SetLineWidth(3.0);  
  h1_E->SetLineColor(kBlack);
  h1_E->SetLineWidth(3.0);  
  //
  h1_E->Draw();
  h1_Ew->Draw("same");
  //
  TCanvas *c2 = new TCanvas("c2","c2",10,10,600,600);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetStatColor(kWhite);
  //gStyle->SetOptStat(kFALSE);
  gr_w_vs_E->Fit(f_weight, "",  "", Emin, Emax);
  gr_w_vs_E->Draw("AP");  
  //
  return 0;
}

void read_file(TString name, TH1D *h1, TH1D *h1_w, TH1D *h1_wo, TGraphErrors *gr){
  std::ifstream fFile(name.Data());
  Double_t E,w;
  if(fFile.is_open()){
    while(fFile>>E>>w){
      //h1->Fill(E,w);
      h1->Fill(E);
      h1_w->Fill(E,w);
      h1_wo->Fill(E);
      gr->SetPoint(gr->GetN(),E,w);
      gr->SetPointError(gr->GetN()-1,E*0.001,w*0.001);
    }
    fFile.close();
  }
}
