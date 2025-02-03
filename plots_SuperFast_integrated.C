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

Int_t plots_SuperFast_integrated(){
  //
  TString fileN01;
  fileN01 = "./hist_superfast_gamma_on_nsb_1x.root";
  //
  TFile *f1 = new TFile(fileN01.Data());
  //TH1D *h1_p2p2_spectrum = (TH1D*)f1->Get("h1_npe_per_ch_max_norm_soft_integrated");
  //TH1D *h1_p2p7_spectrum = (TH1D*)f1->Get("h1_npe_per_ch_max_norm_integrated");
  TH1D *h1_p2p2_spectrum = (TH1D*)f1->Get("h1_npe_per_ch_max_normsim_soft_integrated");
  TH1D *h1_p2p7_spectrum = (TH1D*)f1->Get("h1_npe_per_ch_max_normsim_integrated");
  //
  TCanvas *c1 = new TCanvas("c1",fileN01.Data(),10,10,600,600);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetStatColor(kWhite);
  gStyle->SetOptStat(kFALSE);
  //
  gPad->SetLogx();
  gPad->SetLogy();
  gPad->SetGridx();
  gPad->SetGridy();
  //
  h1_p2p2_spectrum->SetLineWidth(2.0);
  h1_p2p2_spectrum->SetLineColor(kRed+2);
  h1_p2p2_spectrum->SetMarkerColor(kRed+2);
  h1_p2p7_spectrum->SetLineWidth(2.0);
  h1_p2p7_spectrum->SetLineColor(kBlack);
  h1_p2p7_spectrum->SetMarkerColor(kBlack);
  //
  h1_p2p2_spectrum->SetMaximum(1);
  //
  h1_p2p2_spectrum->SetTitle("");
  h1_p2p2_spectrum->Draw();
  h1_p2p7_spectrum->Draw("sames");
  //
  h1_p2p2_spectrum->GetYaxis()->SetTitle("Cumulative event fraction");
  h1_p2p2_spectrum->GetXaxis()->SetTitle("Maximum number of p.e. per channel");
  //
  TLegend *leg = new TLegend(0.6,0.6,0.9,0.9,"","brNDC");
  leg->AddEntry(h1_p2p2_spectrum, "1/E^2.20 (hard)", "apl");
  leg->AddEntry(h1_p2p7_spectrum, "1/E^2.68 (soft)", "apl");
  leg->Draw();  
  //
  return 0;
}
