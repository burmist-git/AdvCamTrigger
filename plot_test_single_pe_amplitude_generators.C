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

Int_t plot_test_single_pe_amplitude_generators(){

  TString fileN01;
  //
  fileN01 = "./hist_test_single_pe_amplitude_generator.root";
  //
  TFile *f01 = new TFile(fileN01.Data());
  //
  TH1D *h1_01 = (TH1D*)f01->Get("_h1_wf_ampl_ADC");
  TH1D *h1_02 = (TH1D*)f01->Get("h1_single_pe_amplitude_generator");
  TH1D *h1_03 = (TH1D*)f01->Get("h1_single_pe_amplitude_from_hist_generator");
  TH1D *h1_04 = (TH1D*)f01->Get("h1_single_pe_amplitude_invf_generate");
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
  h1_04->SetLineColor(kGreen+2);
  h1_04->SetLineWidth(3.0);
  //
  h1_01->SetTitle("");
  //
  h1_01->Draw();
  h1_02->Draw("sames");
  h1_03->Draw("sames");
  h1_04->Draw("sames");
  //
  h1_01->GetXaxis()->SetTitle("ADC counts");
  gPad->SetLogy();
  gPad->SetGridx();
  gPad->SetGridy();
  //
  TLegend *leg = new TLegend(0.6,0.6,0.9,0.9,"","brNDC");
  leg->AddEntry(h1_01, "Single p.e. dist", "apl");
  leg->AddEntry(h1_02, "2 Gb array", "apl");
  leg->AddEntry(h1_03, "From histogram", "apl");
  leg->AddEntry(h1_04, "Inverse function methode", "apl");
  leg->Draw();
  //  
  return 0;
}
