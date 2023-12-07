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

Int_t plots_npe_gamma_gamma_diff_ele_diff(){
  //
  TString fileN01;
  TString fileN02;
  TString fileN03;
  TString fileN04;
  TString fileN05;
  TString fileN06;
  TString fileN07;
  TString fileN08;
  //
  fileN01 = "./hist_fast_gamma_on_nsb_1x_42987507ev.root";
  fileN02 = "./hist_fast_gamma_on_nsb_1x_cut_15918ev.root";
  fileN03 = "./hist_fast_gamma_diffuse_nsb_1x_138771820ev.root";
  fileN04 = "./hist_fast_gamma_diffuse_nsb_1x_cut_12674ev.root";
  fileN05 = "./hist_fast_electron_nsb_1x_48347162ev.root";
  fileN06 = "./hist_fast_electron_nsb_1x_cut_5704ev.root";
  fileN07 = "./hist_fast_proton_nsb_1x_112643300ev.root";
  fileN08 = "./hist_fast_proton_nsb_1x_cut_10172_ev.root";
  //
  TFile *f01 = new TFile(fileN01.Data());
  TFile *f02 = new TFile(fileN02.Data());
  TFile *f03 = new TFile(fileN03.Data());
  TFile *f04 = new TFile(fileN04.Data());
  TFile *f05 = new TFile(fileN05.Data());
  TFile *f06 = new TFile(fileN06.Data());
  TFile *f07 = new TFile(fileN07.Data());
  TFile *f08 = new TFile(fileN08.Data());
  //
  TH1D *h1_01 = (TH1D*)f01->Get("h1_n_pe");
  TH1D *h1_02 = (TH1D*)f02->Get("h1_n_pe");
  TH1D *h1_03 = (TH1D*)f03->Get("h1_n_pe");
  TH1D *h1_04 = (TH1D*)f04->Get("h1_n_pe");
  TH1D *h1_05 = (TH1D*)f05->Get("h1_n_pe");
  TH1D *h1_06 = (TH1D*)f06->Get("h1_n_pe");
  TH1D *h1_07 = (TH1D*)f07->Get("h1_n_pe");
  TH1D *h1_08 = (TH1D*)f08->Get("h1_n_pe");
  //
  h1_01->Rebin(10);
  h1_02->Rebin(10);
  //
  h1_03->Rebin(10);
  h1_04->Rebin(10);
  //
  h1_05->Rebin(10);
  h1_06->Rebin(10);
  //
  h1_07->Rebin(10);
  h1_08->Rebin(10);
  //
  TCanvas *c1 = new TCanvas("c1","c1",10,10,600,600);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetStatColor(kWhite);
  //gStyle->SetOptStat(kFALSE); 
  //
  h1_01->SetLineColor(kBlack);
  h1_01->SetMarkerColor(kBlack);
  h1_01->SetLineWidth(3.0);
  h1_02->SetLineColor(kRed+2);
  h1_02->SetMarkerColor(kRed+2);
  h1_02->SetLineWidth(3.0);
  //
  h1_03->SetLineColor(kBlue+2);
  h1_03->SetMarkerColor(kBlue+2);
  h1_03->SetLineWidth(3.0);
  h1_04->SetLineColor(kMagenta+2);
  h1_04->SetMarkerColor(kMagenta+2);
  h1_04->SetLineWidth(3.0);
  //
  h1_05->SetLineColor(kRed);
  h1_05->SetMarkerColor(kRed);
  h1_05->SetLineWidth(3.0);
  h1_06->SetLineColor(kGreen+2);
  h1_06->SetMarkerColor(kGreen+2);
  h1_06->SetLineWidth(3.0);
  //
  h1_07->SetLineColor(kYellow+2);
  h1_07->SetMarkerColor(kYellow+2);
  h1_07->SetLineWidth(3.0);
  h1_08->SetLineColor(kBlue);
  h1_08->SetMarkerColor(kBlue);
  h1_08->SetLineWidth(3.0);
  //
  h1_01->SetTitle("");
  //
  h1_01->SetMaximum(1.0e+9);
  h1_01->SetMinimum(1.0e+0);
  h1_01->Draw();
  h1_02->Draw("sames");
  //
  h1_03->Draw("sames");
  h1_04->Draw("sames");
  //
  h1_05->Draw("sames");
  h1_06->Draw("sames");
  //
  h1_07->Draw("sames");
  h1_08->Draw("sames");
  //
  //h1_01->GetXaxis()->SetTitle("Number of photons on the telescope sphere");
  h1_01->GetXaxis()->SetTitle("Number of photo electrons");
  gPad->SetLogy();
  gPad->SetGridx();
  gPad->SetGridy();
  //
  TLegend *leg = new TLegend(0.6,0.6,0.9,0.9,"","brNDC");
  leg->AddEntry(h1_01, "gamma on axis  all", "apl");
  leg->AddEntry(h1_02, "gamma on axis  cut", "apl");
  leg->AddEntry(h1_03, "gamma diffused all", "apl");
  leg->AddEntry(h1_04, "gamma diffused cut", "apl");
  leg->AddEntry(h1_05, "e- diffused all", "apl");
  leg->AddEntry(h1_06, "e- diffused cut", "apl");
  leg->AddEntry(h1_07, "p diffused all", "apl");
  leg->AddEntry(h1_08, "p diffused cut", "apl");
  leg->Draw();
  //  
  return 0;
}
