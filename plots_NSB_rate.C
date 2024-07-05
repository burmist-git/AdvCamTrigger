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

Int_t plots_NSB_rate(){

  TString fileN01;
  TString fileN02;
  TString fileN03;
  //
  fileN01 = "../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/nsb_1x_268MHz/trgB/0000/hist_trgB_corsika_0000ID_10ev.root";
  fileN02 = "../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/nsb_1x_268MHz/trgB/0000/hist_trgB_corsika_0000ID_100ev.root";
  fileN03 = "../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/nsb_1x_268MHz/trgB/0000/hist_trgB_corsika_0000ID_1000ev.root";
  //
  TFile *f01 = new TFile(fileN01.Data());
  TFile *f02 = new TFile(fileN02.Data());
  TFile *f03 = new TFile(fileN03.Data());
  //
  TH1D *h1_01 = (TH1D*)f01->Get("h1_digital_sum_rate");
  TH1D *h1_02 = (TH1D*)f02->Get("h1_digital_sum_rate");
  TH1D *h1_03 = (TH1D*)f03->Get("h1_digital_sum_rate");
  //
  //h1_01->Rebin(10);
  //h1_02->Rebin(10);
  //h1_03->Rebin(10);
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
  h1_03->SetLineColor(kBlue+2);
  h1_03->SetLineWidth(3.0);
  //
  h1_01->SetTitle("");
  //
  h1_01->SetMaximum(1.0e+13);
  h1_01->SetMinimum(1.0e+0);
  h1_01->Draw();
  h1_02->Draw("sames");
  //
  h1_03->Draw("sames");
  //
  h1_01->GetXaxis()->SetTitle("Number of photons on the telescope sphere");
  gPad->SetLogy();
  gPad->SetGridx();
  gPad->SetGridy();
  //
  TLegend *leg = new TLegend(0.6,0.6,0.9,0.9,"","brNDC");
  leg->AddEntry(h1_01, "10 eV", "apl");
  leg->AddEntry(h1_02, "100 eV", "apl");
  leg->AddEntry(h1_03, "1000 eV", "apl");
  leg->Draw();
  //  
  return 0;
}
