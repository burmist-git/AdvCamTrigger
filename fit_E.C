//C, C++
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

#include <time.h>

Double_t e_pdf(Double_t *x, Double_t *par){
  Double_t A  = par[0];
  Double_t B  = par[1];
  Double_t C  = par[2];
  Double_t val = A/(TMath::Power(x[0],B)) + C;
  return val;
}

Int_t fit_E(){
  //
  TString fileN01;
  //fileN01 = "./hist_fast_gamma_on_nsb_1x_cut_15918ev.root";
  //fileN01 = "./hist_fast_gamma_diffuse_nsb_1x_cut_12674ev.root";
  //fileN01 = "./hist_fast_electron_nsb_1x_cut_5704ev.root";  
  fileN01 = "./hist_fast_proton_nsb_1x_cut_10172_ev.root";
  //
  //fileN01 = "./hist_fast_proton_nsb_1x.root";
  //fileN01 = "./hist_fast_gamma_diffuse_nsb_1x_cut.root";
  //fileN01 = "./hist_fast_gamma_on_nsb_1x_cut_12130ev.root";
  //fileN01 = "./hist_fast_electron_nsb_1x.root";
  TFile *f01 = new TFile(fileN01.Data());
  TH1D *h1_E_proton = (TH1D*)f01->Get("h1_energy");
  //
  h1_E_proton->SetLineColor(kBlack);
  h1_E_proton->SetLineWidth(2.0);
  h1_E_proton->SetTitle("");
  //
  h1_E_proton->Rebin(100);
  //
  TString fileN;
  //
  const Int_t n = 10000;
  Double_t I[n];
  Double_t etot[n];
  //
  //Double_t etot_min = 0.3;  //in TeV
  //Double_t etot_max = 5;    //in TeV
  Double_t etot_min = 4;    //in TeV
  //Double_t etot_max = 50; //in TeV
  Double_t etot_max = 100; //in TeV
  ///////////////////////////
  //Fit
  const Int_t npar = 3;
  Double_t inParameters[npar];
  Double_t e_min  = etot_min;
  Double_t e_max  = etot_max;
  inParameters[1] = 2.0;
  inParameters[2] = 0.0;
  inParameters[0] = 1300.0;
  //
  TF1 *f_e_pdf = new TF1( "e_pdf", e_pdf, e_min, e_max, npar);
  f_e_pdf->SetParameters(inParameters);
  f_e_pdf->SetParName(0, "A");
  f_e_pdf->SetParName(1, "B");
  f_e_pdf->SetParName(2, "C");
  f_e_pdf->FixParameter(2,inParameters[2]);
  h1_E_proton->Fit("e_pdf","","",e_min, e_max);
  ///////////////////////////  
  TCanvas *c1 = new TCanvas("c1",fileN.Data(),10,10,600,600);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetStatColor(kWhite);
  //gStyle->SetOptStat(kFALSE);
  //
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetLogx();
  gPad->SetLogy();
  //
  h1_E_proton->SetTitle("A/E^{B}");
  h1_E_proton->Draw();
  h1_E_proton->GetXaxis()->SetTitle("E, TeV");
  return 0;
}
