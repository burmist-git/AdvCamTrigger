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

Int_t plots_fadc_val(){

  TString fileN01;
  TString fileN02;
  //
  fileN01 = "../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/nsb_1x_268MHz/trgB/0000/hist_trgB_corsika_0000ID.root";
  fileN02 = "../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/proton_nsb_1x/trgB/0000/hist_trgB_corsika_0000ID.root";
  //
  TFile *f01 = new TFile(fileN01.Data());
  TFile *f02 = new TFile(fileN02.Data());
  //
  TH1D *h1_01 = (TH1D*)f01->Get("h1_fadc_val");
  TH1D *h1_02 = (TH1D*)f02->Get("h1_fadc_val");
  //
  TCanvas *c1 = new TCanvas("c1","c1",10,10,600,600);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetStatColor(kWhite);
  //gStyle->SetOptStat(kFALSE); 
  //
  h1_01->SetLineColor(kBlack);
  h1_01->SetLineWidth(3.0);
  h1_02->SetLineColor(kRed+2);
  h1_02->SetLineWidth(3.0);
  //
  h1_01->SetTitle("");
  //
  //h1_01->SetMaximum(1.0e+13);
  //h1_01->SetMinimum(1.0e+0);
  h1_01->Draw();
  h1_02->Draw("sames");
  //
  h1_01->GetXaxis()->SetTitle("Number of photons on the telescope sphere");
  gPad->SetLogy();
  gPad->SetGridx();
  gPad->SetGridy();
  //
  TLegend *leg = new TLegend(0.6,0.6,0.9,0.9,"","brNDC");
  leg->AddEntry(h1_01, "NSB", "apl");
  leg->AddEntry(h1_02, "p", "apl");
  leg->Draw();
  //  
  return 0;
}
