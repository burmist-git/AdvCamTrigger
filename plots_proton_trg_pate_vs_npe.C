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

Int_t plots_proton_trg_pate_vs_npe(){

  TString fileN01;
  //
  fileN01 = "./hist_proton_st.root";
  //
  TFile *f01 = new TFile(fileN01.Data());
  //
  TH1D *h1_01 = (TH1D*)f01->Get("h1_proton_trg_pate_vs_OR_npe");
  TH1D *h1_02 = (TH1D*)f01->Get("h1_proton_trg_pate_vs_LST1_npe");
  TH1D *h1_03 = (TH1D*)f01->Get("h1_proton_trg_pate_vs_LST1_npe_PMT");
  TH1D *h1_04 = (TH1D*)f01->Get("h1_proton_trg_pate_vs_stereo_npe");
  //
  TCanvas *c1 = new TCanvas("c1","c1",10,10,600,600);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetStatColor(kWhite);
  gStyle->SetOptStat(kFALSE); 
  //
  h1_01->SetLineColor(kBlack);
  h1_01->SetLineWidth(3.0);
  h1_02->SetLineColor(kRed+2);
  h1_02->SetLineWidth(3.0);
  h1_03->SetLineColor(kBlue+2);
  h1_03->SetLineWidth(3.0);
  h1_04->SetLineColor(kMagenta+2);
  h1_04->SetLineWidth(3.0);
  //
  h1_01->SetTitle("");
  //
  h1_01->SetMaximum(200000);
  h1_01->SetMinimum(1000);
  h1_01->Draw();
  h1_02->Draw("sames");
  h1_03->Draw("sames");
  h1_04->Draw("sames");
  //
  h1_01->GetXaxis()->SetTitle("Number of p.e. per telescope");
  h1_01->GetYaxis()->SetTitle("Rate, Hz");
  //
  gPad->SetLogy();
  gPad->SetGridx();
  gPad->SetGridy();
  TLegend *leg = new TLegend(0.6,0.6,0.9,0.9,"","brNDC");
  leg->AddEntry(h1_01, "4 - LSTs in OR", "apl");
  leg->AddEntry(h1_02, "LST1", "apl");
  leg->AddEntry(h1_03, "LST1 PMT", "apl");
  leg->AddEntry(h1_04, "4 - LSTs in stereo (min. 2 telescopes)", "apl");
  leg->Draw();

  //
  /*
  */
  //
  //LST1 -> LST2 -75.1
  //LST1 -> LST3 -212.2
  //LST1 -> LST4 -150.7
  //LST2 -> LST3 -137.1
  //LST2 -> LST4 -75.6
  //LST3 -> LST4  61.5
  //  
  return 0;
}
