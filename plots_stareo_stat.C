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

Int_t plots_stareo_stat(){

  TString fileN01;
  //
  fileN01 = "./hist_proton_st.root";
  //
  TFile *f01 = new TFile(fileN01.Data());
  //
  TH1D *h1_01 = (TH1D*)f01->Get("h1_n_pe_LST_all");
  TH1D *h1_02 = (TH1D*)f01->Get("h1_n_pe_LST_all_LST01");
  TH1D *h1_03 = (TH1D*)f01->Get("h1_n_pe_LST_all_LST02");
  TH1D *h1_04 = (TH1D*)f01->Get("h1_n_pe_LST_all_LST03");
  TH1D *h1_05 = (TH1D*)f01->Get("h1_n_pe_LST_all_LST04");
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
  h1_03->SetLineColor(kBlue+2);
  h1_03->SetLineWidth(3.0);
  h1_04->SetLineColor(kMagenta+2);
  h1_04->SetLineWidth(3.0);
  h1_05->SetLineColor(kGreen+2);
  h1_05->SetLineWidth(3.0);
  //
  h1_01->SetTitle("");
  //
  //h1_01->SetMaximum(1.0e+9);
  //h1_01->SetMinimum(1.0e+0);
  h1_01->Draw();
  h1_02->Draw("sames");
  h1_03->Draw("sames");
  h1_04->Draw("sames");
  h1_05->Draw("sames");
  //
  //h1_01->GetXaxis()->SetTitle("Number of photons on the telescope sphere");
  gPad->SetLogy();
  gPad->SetGridx();
  gPad->SetGridy();
  //
  /*
  TLegend *leg = new TLegend(0.6,0.6,0.9,0.9,"","brNDC");
  leg->AddEntry(h1_01, "gamma on axis  all", "apl");
  leg->AddEntry(h1_02, "gamma on axis  cut", "apl");
  leg->AddEntry(h1_03, "gamma diffused all", "apl");
  leg->AddEntry(h1_04, "gamma diffused cut", "apl");
  leg->Draw();
  */
  //  
  return 0;
}
