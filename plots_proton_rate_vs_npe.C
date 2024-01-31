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

Int_t plots_proton_rate_vs_npe(){

  TString fileN01;
  TString fileN02;
  //
  fileN01 = "./hist_fast_scan_proton_nsb_1x.root";
  //
  TFile *f01 = new TFile(fileN01.Data());
  //
  TH2Poly *h_01 = (TH2Poly*)f01->Get("evH_flux");
  TH2Poly *h_02 = (TH2Poly*)f01->Get("evH_flux_10npe");
  TH2Poly *h_03 = (TH2Poly*)f01->Get("evH_flux_20npe");
  TH2Poly *h_04 = (TH2Poly*)f01->Get("evH_flux_30npe");
  TH2Poly *h_05 = (TH2Poly*)f01->Get("evH_flux_40npe");
  TH2Poly *h_06 = (TH2Poly*)f01->Get("evH_flux_50npe");
  TH2Poly *h_07 = (TH2Poly*)f01->Get("evH_flux_60npe");
  TH2Poly *h_08 = (TH2Poly*)f01->Get("evH_flux_70npe");
  TH2Poly *h_09 = (TH2Poly*)f01->Get("evH_flux_80npe");
  TH2Poly *h_10 = (TH2Poly*)f01->Get("evH_flux_90npe");
  TH2Poly *h_11 = (TH2Poly*)f01->Get("evH_flux_100npe");
  //
  TH2Poly *h_02_w = (TH2Poly*)f01->Get("evH_flux_10npe_w");
  TH2Poly *h_03_w = (TH2Poly*)f01->Get("evH_flux_20npe_w");
  TH2Poly *h_04_w = (TH2Poly*)f01->Get("evH_flux_30npe_w");
  TH2Poly *h_05_w = (TH2Poly*)f01->Get("evH_flux_40npe_w");
  TH2Poly *h_06_w = (TH2Poly*)f01->Get("evH_flux_50npe_w");
  TH2Poly *h_07_w = (TH2Poly*)f01->Get("evH_flux_60npe_w");
  TH2Poly *h_08_w = (TH2Poly*)f01->Get("evH_flux_70npe_w");
  TH2Poly *h_09_w = (TH2Poly*)f01->Get("evH_flux_80npe_w");
  TH2Poly *h_10_w = (TH2Poly*)f01->Get("evH_flux_90npe_w");
  TH2Poly *h_11_w = (TH2Poly*)f01->Get("evH_flux_100npe_w");
  //
  TCanvas *c1 = new TCanvas("c1","c1",10,10,600,600);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetStatColor(kWhite);
  //gStyle->SetOptStat(kFALSE); 
  //
  TGraph *gr = new TGraph();
  gr->SetPoint(0,1,h_01->Integral());
  gr->SetPoint(1,10,h_02->Integral());
  gr->SetPoint(2,20,h_03->Integral());
  gr->SetPoint(3,30,h_04->Integral());
  gr->SetPoint(4,40,h_05->Integral());
  gr->SetPoint(5,50,h_06->Integral());
  gr->SetPoint(6,60,h_07->Integral());
  gr->SetPoint(7,70,h_08->Integral());
  gr->SetPoint(8,80,h_09->Integral());
  gr->SetPoint(9,90,h_10->Integral());
  gr->SetPoint(9,100,h_11->Integral());
  //
  //gr->SetMaximum(1000000);
  //gr->SetMinimum(1000);
  //
  //
  TGraph *gr_w = new TGraph();
  gr_w->SetPoint(0,10,h_02_w->Integral());
  gr_w->SetPoint(1,20,h_03_w->Integral());
  gr_w->SetPoint(2,30,h_04_w->Integral());
  gr_w->SetPoint(3,40,h_05_w->Integral());
  gr_w->SetPoint(4,50,h_06_w->Integral());
  gr_w->SetPoint(5,60,h_07_w->Integral());
  gr_w->SetPoint(6,70,h_08_w->Integral());
  gr_w->SetPoint(7,80,h_09_w->Integral());
  gr_w->SetPoint(8,90,h_10_w->Integral());
  gr_w->SetPoint(9,100,h_11_w->Integral());
  //
  //gr->SetMaximum(1000000);
  //gr->SetMinimum(1000);
  //
  //gr->Draw("APL");

  TMultiGraph *mg = new TMultiGraph();
  mg->Add(gr);
  mg->Add(gr_w);
  mg->Draw("apl");

  
  /*
  h1_01->SetLineColor(kBlack);
  h1_01->SetLineWidth(3.0);
  h1_02->SetLineColor(kRed);
  h1_02->SetLineWidth(3.0);
  //
  h1_01->SetTitle("");
  //
  h1_01->Draw();
  h1_02->Draw("sames");
  //
  h1_01->GetXaxis()->SetTitle("Number of photons on the telescope sphere");
  gPad->SetLogy();
  gPad->SetGridx();
  gPad->SetGridy();
  //
  TLegend *leg = new TLegend(0.6,0.6,0.9,0.9,"","brNDC");
  leg->AddEntry(h1_01, "all", "apl");
  leg->AddEntry(h1_02, "cut", "apl");
  leg->Draw();
  */
  //  
  return 0;
}
