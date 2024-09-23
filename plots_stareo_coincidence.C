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

Int_t plots_stareo_coincidence(){

  TString fileN01;
  //
  fileN01 = "./hist_proton_st.root";
  //
  TFile *f01 = new TFile(fileN01.Data());
  //
  TH1D *h1_01 = (TH1D*)f01->Get("h1_n_pe_LST_all");
  TH1D *h1_02 = (TH1D*)f01->Get("h1_n_pe_LST_single");
  TH1D *h1_03 = (TH1D*)f01->Get("h1_n_pe_LST_double");
  TH1D *h1_04 = (TH1D*)f01->Get("h1_n_pe_LST_triple");
  TH1D *h1_05 = (TH1D*)f01->Get("h1_n_pe_LST_four");
  //
  TH1D *h1_06 = (TH1D*)f01->Get("h1_n_pe_LST_all_LST01");
  //
  //
  Double_t rate50kHz = h1_06->GetBinContent(1+1);
  cout<<"rate50kHz = "<<rate50kHz<<endl;
  //
  Double_t rate50kHz_at20pe = 51830.3;
  //
  h1_01->SetBinContent(1+1,h1_01->GetBinContent(1+1)/rate50kHz*rate50kHz_at20pe);
  //
  h1_02->SetBinContent(1+1,h1_02->GetBinContent(1+1)/rate50kHz*rate50kHz_at20pe);
  h1_02->SetBinContent(1+2,h1_02->GetBinContent(1+2)/rate50kHz*rate50kHz_at20pe);
  h1_02->SetBinContent(1+3,h1_02->GetBinContent(1+3)/rate50kHz*rate50kHz_at20pe);
  h1_02->SetBinContent(1+4,h1_02->GetBinContent(1+4)/rate50kHz*rate50kHz_at20pe);
  //
  h1_03->SetBinContent(1+1,h1_03->GetBinContent(1+1)/rate50kHz*rate50kHz_at20pe);
  h1_03->SetBinContent(1+2,h1_03->GetBinContent(1+2)/rate50kHz*rate50kHz_at20pe);
  h1_03->SetBinContent(1+3,h1_03->GetBinContent(1+3)/rate50kHz*rate50kHz_at20pe);
  h1_03->SetBinContent(1+4,h1_03->GetBinContent(1+4)/rate50kHz*rate50kHz_at20pe);
  h1_03->SetBinContent(1+5,h1_03->GetBinContent(1+5)/rate50kHz*rate50kHz_at20pe);
  h1_03->SetBinContent(1+6,h1_03->GetBinContent(1+6)/rate50kHz*rate50kHz_at20pe);
  //
  h1_04->SetBinContent(1+1,h1_04->GetBinContent(1+1)/rate50kHz*rate50kHz_at20pe);
  h1_04->SetBinContent(1+2,h1_04->GetBinContent(1+2)/rate50kHz*rate50kHz_at20pe);
  h1_04->SetBinContent(1+3,h1_04->GetBinContent(1+3)/rate50kHz*rate50kHz_at20pe);
  h1_04->SetBinContent(1+4,h1_04->GetBinContent(1+4)/rate50kHz*rate50kHz_at20pe);
  //
  h1_05->SetBinContent(1+1,h1_05->GetBinContent(1+1)/rate50kHz*rate50kHz_at20pe);
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
  h1_03->SetMaximum(1.0e+6);
  h1_03->SetMinimum(0.1e+0);
  h1_03->Draw();
  h1_02->Draw("sames");
  h1_01->Draw("sames");
  h1_04->Draw("sames");
  h1_05->Draw("sames");
  //
  //h1_01->GetXaxis()->SetTitle("Number of photons on the telescope sphere");
  gPad->SetLogy();
  gPad->SetGridx();
  gPad->SetGridy();
  //
  TLegend *leg = new TLegend(0.6,0.6,0.9,0.9,"","brNDC");
  leg->AddEntry(h1_01, "<OR> all LSTs @ 20p.e.", "apl");
  leg->AddEntry(h1_02, "Isolated LST @ 20p.e.", "apl");
  leg->AddEntry(h1_03, "Double coincidence @ 20p.e.", "apl");
  leg->AddEntry(h1_04, "Triple coincidence @ 20p.e.", "apl");
  leg->AddEntry(h1_05, "<AND> all LSTs @ 20p.e.", "apl");
  leg->Draw();

  //  
  return 0;
}
