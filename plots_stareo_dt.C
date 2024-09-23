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

Int_t plots_stareo_dt(){

  TString fileN01;
  //
  fileN01 = "./hist_proton_st.root";
  //
  TFile *f01 = new TFile(fileN01.Data());
  //
  TH1D *h1_01 = (TH1D*)f01->Get("h1_dtime_LST1_m_LST2_norm");
  TH1D *h1_02 = (TH1D*)f01->Get("h1_dtime_LST1_m_LST3_norm");
  TH1D *h1_03 = (TH1D*)f01->Get("h1_dtime_LST1_m_LST4_norm");
  TH1D *h1_04 = (TH1D*)f01->Get("h1_dtime_LST2_m_LST3_norm");
  TH1D *h1_05 = (TH1D*)f01->Get("h1_dtime_LST2_m_LST4_norm");
  TH1D *h1_06 = (TH1D*)f01->Get("h1_dtime_LST3_m_LST4_norm");
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
  h1_06->SetLineColor(kYellow+2);
  h1_06->SetLineWidth(3.0);
  //
  h1_01->SetTitle("");
  //
  h1_01->SetMaximum(1.13);
  h1_01->Draw();
  h1_02->Draw("sames");
  h1_03->Draw("sames");
  h1_04->Draw("sames");
  h1_05->Draw("sames");
  h1_06->Draw("sames");

  TLine *l1_01 = new TLine( -75.1,0.0, -75.1, 1.1);
  TLine *l1_02 = new TLine(-212.2,0.0,-212.2, 1.1);
  TLine *l1_03 = new TLine(-150.7,0.0,-150.7, 1.1);
  TLine *l1_04 = new TLine(-137.1,0.0,-137.1, 1.1);
  TLine *l1_05 = new TLine( -75.6,0.0, -75.6, 1.1);
  TLine *l1_06 = new TLine(  61.5,0.0,  61.5, 1.1);
  //
  l1_01->SetLineWidth(3.0);
  l1_02->SetLineWidth(3.0);
  l1_03->SetLineWidth(3.0);
  l1_04->SetLineWidth(3.0);
  l1_05->SetLineWidth(3.0);
  l1_06->SetLineWidth(3.0);
  // 
  l1_01->SetLineStyle(kDashed);
  l1_02->SetLineStyle(kDashed);
  l1_03->SetLineStyle(kDashed);
  l1_04->SetLineStyle(kDashed);
  l1_05->SetLineStyle(kDashed);
  l1_06->SetLineStyle(kDashed);
  //
  l1_01->SetLineColor(kBlack);
  l1_02->SetLineColor(kRed+2);
  l1_03->SetLineColor(kBlue+2);
  l1_04->SetLineColor(kMagenta+2);
  l1_05->SetLineColor(kGreen+2);
  l1_06->SetLineColor(kYellow+2);
  //
  l1_01->Draw("sames");
  l1_02->Draw("sames");
  l1_03->Draw("sames");
  l1_04->Draw("sames");
  l1_05->Draw("sames");
  l1_06->Draw("sames");  
  //
  //h1_01->GetXaxis()->SetTitle("Number of photons on the telescope sphere");
  //gPad->SetLogy();
  gPad->SetGridx();
  gPad->SetGridy();
  //
  TLegend *leg = new TLegend(0.6,0.6,0.9,0.9,"","brNDC");
  leg->AddEntry(h1_01, "LST1 - LST2", "apl");
  leg->AddEntry(h1_02, "LST1 - LST3", "apl");
  leg->AddEntry(h1_03, "LST1 - LST4", "apl");
  leg->AddEntry(h1_04, "LST2 - LST3", "apl");
  leg->AddEntry(h1_05, "LST2 - LST4", "apl");
  leg->AddEntry(h1_06, "LST3 - LST4", "apl");
  leg->Draw();
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
